function B = hss_from_random_sampling(Afun, Afunt, Aeval, m, n)
%HSS_FROM_RANDOM_SAMPLING Build the HSS representation of a matrix
% 			   using mat-vec multiplication with random (block) vectors
%			   and access to (block) diagonal entries.
%
% B = HSS_RANDOM_SAMPLING is based on the algorithm in [1].
%
% [1] Martinsson, Per-Gunnar. "A fast randomized algorithm for computing a
%     hierarchically semiseparable representation of a matrix." SIAM 
%     Journal on Matrix Analysis and Applications 32.4 (2011): 1251-1274

tol = hssoption('threshold');

B = hss_build_hss_tree(m, n, hssoption('block-size'));

failed = true;

k = 15;

Ocol = randn(n, k);
Scol = Afun(Ocol);
Orow = randn(m, k);
Srow = Afunt(Orow);

bs = 5;

nrm = svds(@(x,t) normest_afun(Afun, Afunt, x, t), [n n], 1, 'largest', ...
    struct('tol', 1e-2)); nrm = nrm(1);
% nrm = 1;

while failed
    % fprintf('HSS_RANDOM_FROM_SAMPLING :: columns: %d\n', size(Scol, 2));
    [B, ~, ~, ~, ~, ~, ~, ~, ~, failed] = ...
        hss_from_random_sampling_rec(B, Aeval, Scol, Srow, ...
            Ocol, Orow, 0, 0, tol, nrm);
    
    if failed
        % fprintf('HSS_FROM_RANDOM_SAMPLING :: Enlarging sampling space to %d\n', k + bs);        
        Ocol = [ Ocol, randn(n, bs) ];
        Scol = [ Scol, Afun(Ocol(:, end-bs+1:end))  ];
        Orow = [ Orow, randn(m, bs) ];
        Srow = [ Srow, Afunt(Orow(:,end-bs+1:end)) ];
        k = k + bs;
    end
end

% Remove any temporary data that we might have stored in the leafnodes B12
% anb B21 -- these are not part of the final HSS data structure. 
B = clean_structure(B);

% [ hssrank(B), k ]
% B = hss_proper(B);
% B = hss_compress(B, tol);
% hssrank(B)

end

function [B, Scol, Srow, Ocol, Orow, Jcol, ...
          Jrow, U, V, failed] = ...
              hss_from_random_sampling_rec(B, Aeval, Scol, Srow, ...
              Ocol, Orow, row, col, tol, nrm)
          
	if B.leafnode == 1
        failed = false;
        
		[m, n] = size(B.D);
        
        % Only store the matrix D if it has not been initialized yet; in
        % general it might have been initialized in a previous pass. 
        if isequal(B.D, zeros(m, n))        
            B.D = Aeval(row + 1: row + m, col + 1:col + n);
            B.B12 = [];
            B.B21 = [];
        end
        
        % Compute span of the columns, but only if the old one was not ok
        if isempty(B.B12)
            Scol = Scol - B.D * Ocol;

            [Q, rk] = colspan(Scol, tol, nrm);

            if rk >= size(Scol, 2) - 10
                failed = true;
                Jcol = []; Jrow = []; U = []; V = [];                
                return;
            end

            [Xcol, Jcol] = interpolative(Q', tol);
            B.U = Xcol';
            Scol = Scol(Jcol, :);

            U = Xcol';                
            B.B12 = Jcol;        
            Jcol = row + Jcol;
        else
            U = B.U;
            Jcol = B.B12;
            Scol = Scol(Jcol, :) - B.D(Jcol, :) * Ocol;
            Jcol = row + Jcol;
        end

        % Compute span of the rows
        if isempty(B.B21)
            Srow = Srow - B.D' * Orow;
        
            [Q, rk] = colspan(Srow, tol, nrm);          

            if rk >= size(Srow, 2) - 10
                failed = true;
                Jcol = []; Jrow = []; U = []; V = [];
                return;
            end

            [Xrow, Jrow] = interpolative(Q.', tol);

            B.V = Xrow.';
            Srow = Srow(Jrow, :);

            V = Xrow.';

            B.B21 = Jrow;

            Jrow = col + Jrow;
        else
            V = B.V;
            Jrow = B.B21;
            Srow = Srow(Jrow, :) - B.D(:, Jrow)' * Orow;
            Jrow = col + Jrow;
        end
	else
		[B.A11, Scol1, Srow1, Ocol1, Orow1, Jcol1, Jrow1, U1, V1, failed1]  = hss_from_random_sampling_rec(B.A11, Aeval, Scol(1:B.ml, :), Srow(1:B.nl, :), Ocol(1:B.nl, :), Orow(1:B.ml, :), row, col, tol, nrm);
        
        if ~failed1        
            [B.A22, Scol2, Srow2, Ocol2, Orow2, Jcol2, Jrow2, U2, V2, failed2]  = hss_from_random_sampling_rec(B.A22, Aeval, Scol(B.ml + 1:end, :), Srow(B.nl + 1:end,:), Ocol(B.nl + 1:end, :), Orow(B.ml + 1:end, :), row + B.ml, col + B.nl, tol, nrm);
        end
        
        if (failed1 || failed2)
            failed = true;
            
            Scol = []; Srow = []; Ocol = []; Orow = []; 
            Jcol = []; Jrow = []; U = []; V = [];
            
            return;
        else
            failed = false;
        end
        
        Ocol2 = V2' * Ocol2;
        Ocol1 = V1' * Ocol1;
        Orow2 = U2' * Orow2;
        Orow1 = U1' * Orow1;

        Jcol = [Jcol1, Jcol2]; Jrow = [Jrow1, Jrow2];
        Ocol = [Ocol1; Ocol2]; Orow = [Orow1; Orow2];
        
        if isempty(B.B12)
        
            B.B12 = Aeval(Jcol1, Jrow2);
            B.B21 = Aeval(Jcol2, Jrow1);            

            Scol = [Scol1 - B.B12  * Ocol2;  Scol2 - B.B21  * Ocol1 ]; 
            Srow = [Srow1 - B.B21' * Orow2;  Srow2 - B.B12' * Orow1 ];

            if B.topnode == 0                
                [Q, R, ~] = qr(Scol, 0); rk = sum(abs(diag(R)) > abs(R(1,1)) * tol * size(R, 2));

                Q = Q(:,1:rk);

                [Xcol, Jcolloc] = interpolative(Q', tol);

                B.Rl = Xcol(:, 1:size(Scol1, 1))';
                B.Rr = Xcol(:, size(Scol1, 1)+1:end)';
                Scol = Scol(Jcolloc, :);
                Jcol = Jcol(Jcolloc);
                U = [ B.Rl ; B.Rr ];

                [Q, R, ~] = qr(Srow, 0); rk = sum(abs(diag(R)) > abs(R(1,1)) * tol * size(R, 2));

                Q = Q(:,1:rk);

                [Xrow, Jrowloc] = interpolative(Q', tol);

                B.Wl = Xrow(:, 1:size(Srow1, 1))';
                B.Wr = Xrow(:, size(Srow1, 1)+1:end)';
                Srow = Srow(Jrowloc, :);
                Jrow = Jrow(Jrowloc);
                V = [ B.Wl ; B.Wr ];
                
                B.U = Jcolloc;
                B.V = Jrowloc;
            else
                Scol = []; Srow = []; Ocol = []; Orow = []; 
                Jcol = []; Jrow = []; U = []; V = [];
            end
        else
            Jcolloc = B.U;
            Jrowloc = B.V;
            
            Scol = [Scol1 - B.B12  * Ocol2;  Scol2 - B.B21  * Ocol1 ]; 
            Srow = [Srow1 - B.B21' * Orow2;  Srow2 - B.B12' * Orow1 ];
            
            U = [ B.Rl ; B.Rr ];
            V = [ B.Wl ; B.Wr ];
            
            Scol = Scol(Jcolloc, :);
            Jcol = Jcol(Jcolloc);
            Srow = Srow(Jrowloc, :);
            Jrow = Jrow(Jrowloc);
        end
	end
end

function [Q, rk] = colspan(S, tol, nrm)
    use_qr = true;
    
    if use_qr
        [Q, R, ~] = qr(S, 0);
        rk = sum(abs(diag(R)) > nrm * tol * sqrt(size(S, 2)));
        Q = Q(:,1:rk);
    else
        [Q, S] = svd(S);
        rk = sum(diag(S) > nrm * tol);
        Q = Q(:,1:rk);
    end
end

function y = normest_afun(Afun, Afunt, x, transp)

if strcmp(transp, 'transp')
    y = Afunt(x);
else
    y = Afun(x);
end

end

function B = clean_structure(B)
%CLEAN_STRUCTURE Clean the HSS structure from temporary data. 

if B.leafnode
    B.B12 = [];
    B.B21 = [];
else
    B.A11 = clean_structure(B.A11);
    B.A22 = clean_structure(B.A22);
end

end