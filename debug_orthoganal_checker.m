% a function to calculate <psi1|psi2>, and check the orthogonality
function corr=debug_orthoganal_checker(psi_schrodinger)

    global triangle_areas_schrod;
    global t_schrod;
    
    N = size(psi_schrodinger, 2);
    % get the normalization
    [normalization, ~] = normalize_and_sqrt_m_triangular(psi_schrodinger);
    
    psi_normalized = psi_schrodinger .* normalization;
    
    integral = zeros(N);
    for i = 1 : N
        for j = i : N
            psi_ij =  psi_normalized(:, i) .* psi_normalized(:, j);
           
            avg_psi_ij = 1/3 * ...
                ( psi_ij(t_schrod(1, :)) + psi_ij(t_schrod(2, :)) + psi_ij(t_schrod(3, :)) );
            integral(i, j) = sum( avg_psi_ij' .* triangle_areas_schrod, 'all' );
        end
    end
    
    figure
    imagesc(integral)
    colorbar
    
    integral
    corr = integral;
    save("./correlation.mat", "corr");
end