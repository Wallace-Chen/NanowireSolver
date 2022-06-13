
function debug_surface_charge(N, epsilon_F, E_schrodinger, psi_schrodinger, V_poiss, heaviside_n_D)

    %% debug purpose
    global triangle_areas_poiss;
    global p_poiss t_poiss e_poiss;
    global l0;
     
    [normalization, ~] = normalize_and_sqrt_m_triangular( psi_schrodinger );
    psi_sqrt_eps_sum = psi_sqrt_eps_summation( N, normalization, psi_schrodinger, epsilon_F, E_schrodinger );
    
    x = p_poiss(1, :);
    y = p_poiss(2, :);
    nD = n_D_func(x, y, heaviside_n_D);
    nE = n_e_func(x, y, psi_sqrt_eps_sum);
    nSpont = sigma_spont(x, y);
    nPiezo = sigma_piezo(x, y);
    nSurface = sigma_surface(x, y, V_poiss);
    
    p1 = find(e_poiss(5, :) == 1);
    ave1_nSurface = 0.5 * (nSurface(e_poiss(1, p1)) + nSurface(e_poiss(2, p1)));
    p1_l = sqrt( (p_poiss(1, e_poiss(1, p1))-p_poiss(1, e_poiss(2, p1))).^2 + (p_poiss(2, e_poiss(1, p1))-p_poiss(2, e_poiss(2, p1))).^2 );
    p2 = find(e_poiss(5, :) == 2);
    ave2_nSurface = 0.5 * (nSurface(e_poiss(1, p2)) + nSurface(e_poiss(2, p2)));
    p2_l = sqrt( (p_poiss(1, e_poiss(1, p2))-p_poiss(1, e_poiss(2, p2))).^2 + (p_poiss(2, e_poiss(1, p2))-p_poiss(2, e_poiss(2, p2))).^2 );
    p3 = find(e_poiss(5, :) == 3);
    ave3_nSurface = 0.5 * (nSurface(e_poiss(1, p3)) + nSurface(e_poiss(2, p3)));
    p3_l = sqrt( (p_poiss(1, e_poiss(1, p3))-p_poiss(1, e_poiss(2, p3))).^2 + (p_poiss(2, e_poiss(1, p3))-p_poiss(2, e_poiss(2, p3))).^2 );
     
    ave_surface = 1/3 * ...
        (nSurface( t_poiss(1, :) ) + ...
         nSurface( t_poiss(2, :) ) + ...
         nSurface( t_poiss(3, :) ));
    ave_nSpont = 1/3 * ...
        (nSpont( t_poiss(1, :) ) + ...
         nSpont( t_poiss(2, :) ) + ...
         nSpont( t_poiss(3, :) ));
    ave_nPiezo = 1/3 * ...
        (nPiezo( t_poiss(1, :) ) + ...
         nPiezo( t_poiss(2, :) ) + ...
         nPiezo( t_poiss(3, :) ));
     ave_nD = 1/3 * ...
        (nD( t_poiss(1, :) ) + ...
         nD( t_poiss(2, :) ) + ...
         nD( t_poiss(3, :) ));
     ave_nE = 1/3 * ...
        (nE( t_poiss(1, :) ) + ...
         nE( t_poiss(2, :) ) + ...
         nE( t_poiss(3, :) ));
     fprintf('total surface charge: %f\n', p1_l*ave1_nSurface' + p2_l*ave2_nSurface' + p3_l*ave3_nSurface');
     fprintf('surface charge: %f\n', triangle_areas_poiss*ave_surface');
     fprintf('piezo: %f\n', triangle_areas_poiss*ave_nPiezo' );
     fprintf('spontaneous: %f\n', triangle_areas_poiss*ave_nSpont' );
     fprintf('dopants: %f\n', triangle_areas_poiss*ave_nD' );
     fprintf('electron density: %f\n', triangle_areas_poiss*ave_nE' );
     
    figure
    pdesurf(l0 * p_poiss, t_poiss,  V_poiss);
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    
    global DX DY;
    figure
    pdesurf(l0 * p_poiss, t_poiss,  DX * V_poiss);
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    
    figure
    pdesurf(l0 * p_poiss, t_poiss,  DY * V_poiss);
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    
    figure
    pdesurf(l0 * p_poiss, t_poiss,  nSurface');
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    
    figure
    pdesurf(l0 * p_poiss, t_poiss,  nPiezo');
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    
    figure
    pdesurf(l0 * p_poiss, t_poiss,  nSpont');
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    
    figure
    pdesurf(l0 * p_poiss, t_poiss,  nD');
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    
    figure
    pdesurf(l0 * p_poiss, t_poiss,  nE');
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    
    figure
    pdesurf(l0 * p_poiss, t_poiss,  nPiezo'+nSpont'+nD'+nE');
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    %% debug end

end