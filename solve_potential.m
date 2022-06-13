% solve electrostat from poisson equations

function [V_poiss] = solve_potential(N, epsilon_F, E_schrodinger, psi_schrodinger, heaviside_n_D, factor, V_poiss)

    %% global variable declaration
    
    global debug;
    global p_poiss; global e_poiss; global t_poiss;
    global C;
    global degree_of_polygon;
    global psi_sqrt_eps_xc_m;
    global degree_of_polygon;

    a_poiss = 0;
    b_poiss = generate_dirichlet_boundary( 1 );
    c_poiss = 'c_coefficient_poiss(sd)';
    %% debug purpose
    %c_poiss = [ mat2str(-vector_of_eps(1)) ];
    % for j = 2 : number_of_interfaces
    %     c_poiss = [c_poiss, ' + ', mat2str(vector_of_eps(j-1) - vector_of_eps(j)), ...
    %         sprintf(' * heaviside_core_poiss(%s, x, y)', mat2str(vector_of_side_lengths(j)))];
    % end
    %% debug end
    
    psi_sqrt_eps_sum = zeros( size(psi_schrodinger(:, 1)) );
    
    if N > 0
        % find normalization and |psi|^2*sqrt(m) term integrated over (x,y), normalized.
        [normalization, ~] = normalize_and_sqrt_m_triangular( psi_schrodinger );
        % calculate the psi_sqrt_eps term needed to calculate the total
        % electron density
        psi_sqrt_eps_sum = psi_sqrt_eps_summation( N, normalization, psi_schrodinger, epsilon_F, E_schrodinger );
    end
    
    
    %% SCALAR PDE EQUATION
    % f coefficient for the scalar PDE equation, solve the electrostatic
    % potentials, due to electrons, dopants, and polarizations. V_cb do not
    % need to be solved (assumed to remain unchanged through the iterations).
    f_poiss_1 = sprintf('n_D_func(x, y, %s) + n_e_func(x, y, %s)', mat2str(heaviside_n_D), mat2str(psi_sqrt_eps_sum));
    if abs(degree_of_polygon) == 3
        f_poiss_1 = sprintf('n_D_func(x, y, %s) + n_e_func(x, y, %s) + sigma_spont(x, y) + sigma_piezo(x, y)', mat2str(heaviside_n_D), mat2str(psi_sqrt_eps_sum));
        %f_poiss_1 = sprintf('n_D_func(x, y, %s) + n_e_func(x, y, %s) + sigma_spont(x, y) + sigma_piezo(x, y) + sigma_surface(x, y, %s)', mat2str(heaviside_n_D), mat2str(psi_sqrt_eps_sum), mat2str(V_poiss));
    end
    if debug
       disp("debug: initiate_system, showing the coefficients for the scalar poisson equation...");
       c_poiss
       a_poiss
       f_poiss_1
    end
    V_poiss = assempde(b_poiss, p_poiss, e_poiss, t_poiss, c_poiss, a_poiss, f_poiss_1);
    
    V_poiss = factor * V_poiss / C;
    if abs(degree_of_polygon) ~= 3
        V_poiss = V_poiss - min(V_poiss, [], 'all');
    end

end