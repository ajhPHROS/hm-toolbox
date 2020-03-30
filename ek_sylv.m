function [Xu, Xv, VA, VB] = ek_sylv(A, B, u, v, k, tol, debug, nrmtype)
%EK_SYLV Approximate the solution of a Sylvester equation AX + XB + U*V' = 0.
%
% [XU,XV] = EK_SYLV(A, B, U, V, K) approximates the solution of the
%     Sylvester equation in the factored form XU * XV'. The integer K is
%     the maximum number of iteration to perform, and the default tolerance
%     is 1e-8. The maximum number of steps K can be set to inf; in that
%     case, the iteration is stopped only when the solution is accurate
%     enough. 
%
% [XU, XV] = EK_SYLV(A, B, U, V, K, TOL, DEBUG) also returns the bases VA
%     and VB, and the optional parameters TOL and DEBUG control the
%     stopping criterion and the debugging during the iteration.
%
% The default stopping criterion is to continue until the residual is
% smaller than TOL * norm(XU * XV'). 
%
% The tolerance TOL can also be specified as a function TOL(R, N) that
% takes as input the residual norm and the norm of the solution (R and N,
% respectively), and returns true if the solution can be accepted.

if ~exist('debug', 'var')
    debug = false;
end

if ~exist('tol', 'var')
    tol = 1e-8;
end

if ~exist('nrmtype', 'var')
    nrmtype = 2;
end
if size(u, 2) ~= size(v, 2)
	error('EK_SYLV:: different number of columns in the factorization of the rhs');
end
if size(u, 2) == 0
    Xu = u;
    Xv = v;
    VA = u;
    VB = v;
    return;
end

if ~isstruct(A)
    if issparse(A)
        AA = ek_struct(A, issymmetric(A));
    else
        AA = ek_struct(A, false);
    end
else
    AA = A;
end

if ~isstruct(B)
    if issparse(B)
        BB = ek_struct(B', issymmetric(B));
    else
        BB = ek_struct(B', false);
    end
else
    BB = B';
end

nrmA = AA.nrm;
nrmB = BB.nrm;

% Dimension of the space
sa = size(u, 2);
sb = size(v, 2);

bsa = sa;
bsb = sb;

% tol can be function tol(r, n) that is given the residual and the norm, or
% a scalar. In the latter case, we turn it into a function
if isfloat(tol)
    tol_eps = tol;
    tol = @(r, n) r < tol_eps * n;
end

it=1;
while max(sa-2*bsa, sb-2*bsb) < k
    if exist('VA', 'var') && ( ( size(VA, 2) + 2 * bsa >= size(VA, 1) ) || ...
            ( size(VB, 2) + 2 * bsb >= size(VB, 1) ) )
        
        %warning('HM:EK_SYLV', ...
        %    'Extended Krylov space is equal to the whole space: using dense solver');
        
        na = size(VA, 1);
        nb = size(VB, 1);
        
        A = AA.multiply(1.0, 0.0, eye(na));
        B = BB.multiply(1.0, 0.0, eye(nb));
        
        Y = lyap(A, B', u * v');
        
        [UU,SS,VV] = svd(Y);
        
	switch nrmtype
    		case 2
        	s = diag(SS);
        	rk = sum( arrayfun(@(ss) tol(ss, s(1) / (nrmA + nrmB)), s) == 0);
    		case 'fro'
        	d = sort(diag(SS));
        	s = sqrt(cumsum(d.^2));
        	rk = sum( arrayfun(@(ss) tol(ss, s(end) / (nrmA + nrmB)), s) == 0 );
	end
        
        Xu = UU(:,1:rk) * sqrt(SS(1:rk,1:rk));
        Xv = VV(:,1:rk) * sqrt(SS(1:rk,1:rk));
        
        if debug
            fprintf('%d Dense solver: rank %d, size %d\n', it, rk, max(na,nb));
        end
        
        return;
    end
    
    if ~exist('VA', 'var')
        [VA, KA, HA, param_A] = ek_krylov(AA, u);
        [VB, KB, HB, param_B] = ek_krylov(BB, v);
    else
        [VA, KA, HA, param_A] = ek_krylov(VA, KA, HA, param_A);
        [VB, KB, HB, param_B] = ek_krylov(VB, KB, HB, param_B);
    end
    
    sa = size(VA, 2);
    sb = size(VB, 2);
    
    % Compute the solution and residual of the projected Lyapunov equation
    As = HA / KA(1:end-bsa,:);
    Bs = HB / KB(1:end-bsb,:);
    Cs = zeros(size(As, 1), size(Bs, 1));
    
    if ~exist('Cprojected', 'var')
        Cprojected = ( VA(:,1:size(u,2))' * u ) * ( v' * VB(:,1:size(v,2)) );
    end
    
    Cs(1:size(u,2), 1:size(v,2)) = Cprojected;
    
    [Y, res] = lyap_galerkin(As, Bs, Cs, bsa, bsb, nrmtype);
    
    % You might want to enable this for debugging purposes
    if debug
        fprintf('%d Residue: %e (space dim: %d)\n', it, res / norm(Y), ...
            size(Y, 1));
    end
    
    if tol(res, norm(Y, nrmtype)) % res < norm(Y) * tol
        break
    end
    it=it+1;
end

% fprintf('SYLV : it = %d\n', it);

% fprintf('lyap its = %d, nrmA = %e\n', it, nrmA)
[UU,SS,VV] = svd(Y);

switch nrmtype
    case 2
        s = diag(SS);
        rk = sum( arrayfun(@(ss) tol(ss, s(1) / (nrmA + nrmB)), s) == 0);
    case 'fro'
        d = sort(diag(SS));
        s = sqrt(cumsum(d.^2));
        rk = sum( arrayfun(@(ss) tol(ss, s(end) / (nrmA + nrmB)), s) == 0 );
end

Xu = VA(:,1:size(Y,1)) * UU(:,1:rk) * sqrt(SS(1:rk,1:rk));
Xv = VB(:,1:size(Y,2)) * VV(:,1:rk) * sqrt(SS(1:rk,1:rk));

end

