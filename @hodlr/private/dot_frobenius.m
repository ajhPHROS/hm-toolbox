function F = dot_frobenius(H1,H2)
%DOT_FROBENIUS Compute the Frobenius inner product tr(H1'*H2).
%
% F = DOT_FROBENIUS(H1,H2) evaluates the Frobenius inner product between H1
% and H2. The computed inner product is exact up to rounding errors.

if xor(~is_leafnode(H1), ~is_leafnode(H2))
    error('DOT_FROBENIUS:: the two hodlr matrices do not have compatible partitioning')
end
if is_leafnode(H1)
    F = sum(H1.F .* H2.F,'all');
else
    if size(H1.U21, 1) ~= size(H2.U21, 1) || size(H1.U12, 1) ~= size(H2.U12, 1)|| ...
            size(H1.V12, 1) ~= size(H2.V12, 1) || size(H2.V21, 1) ~= size(H2.V21, 1)
        error('DOT_FROBENIUS:: the two hodlr matrices do not have compatible partitioning')
    end
    F11 = dot_frobenius(H1.A11, H2.A11);
    F22 = dot_frobenius(H1.A22, H2.A22);
    
    % Contribution from the block (2,1)
    F21 = sum((H1.V21'*H2.V21).*(H1.U21'*H2.U21));

    % Contribution from the block (1,2)
    F12 = sum((H1.V12'*H2.V12).*(H1.U12'*H2.U12));
    
    F = sum([ F11, F12, F21, F22 ]);
end


end

