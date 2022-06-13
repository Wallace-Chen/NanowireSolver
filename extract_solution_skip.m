% pick wavefunctions & energies from solutions of system of eigen-equations
function [psi, E] = extract_solution_skip(N, v, l)
    
    n_col = floor( (length(l)-1)/N ) + 1;
    % skip every N solutions
    indices = zeros(N, n_col);
    indices(1, :) = 1;
    v = v(:, find(indices));
    l = l(find(indices));
    n_ele = min(n_col, N);
    % the number of mesh nodes
    Np = length( v(:,1) ) / N;
    % build mask matrix to speed up the solution extraction
    cells = cellmat(1, n_ele, Np, 1, 1);
    mask = blkdiag( cells{:} );
    
    % extract psi from eigenvectors v, pick the diagonal Np*1 block
    % matrices of v as wavefunctions:
    % psi_11 psi_12 psi_13
    % psi_21 psi_22 psi_23
    % psi_31 psi_32 psi_33
    % then chosen wavefunctions are [psi_11; psi_22; psi_33] with
    % eigen-energies increasing in order.
    psi_tmp = sum( v(1:Np*n_ele, 1:n_ele) .* mask, 2 );
    psi_tmp = reshape(psi_tmp, Np, n_ele);
    
    psi = zeros(Np, max(n_ele, n_col));
    psi(:, 1:n_ele) = psi_tmp;
    for i = N+1 : n_col
        psi( :, i ) = v( end-Np+1:end, i);
    end
    % slice the corresponding eigen-energies
    E = l;
    
end