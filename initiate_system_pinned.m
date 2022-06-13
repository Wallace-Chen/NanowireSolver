% function to initiate the system, solve the indepedent schrodinger
% equation, find the proper number of electrons satisfying the charge
% neutrality equation.

% the fermi level is pinned
function [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N] = initiate_system_pinned()

    global debug;
    
    global p_schrod; global e_schrod; global t_schrod;
    global p_poiss; global e_poiss; global t_poiss;
    global degree_of_polygon;
    global vector_of_eps;
    global number_of_interfaces;
    global vector_of_side_lengths;
    global C;
    global psi_sqrt_eps_xc_m;
    global pinned_fermi_level;
    global damping_factor;
    global save_file;
    
    global include_ex;
    global K_poiss N_poiss;
    global Q G H R M_schrod;
    global DX DY;
    
    disp("Initializing the system...");
    %% setup coefficient matrix
    % boundary matrix
    bmatrix = generate_dirichlet_boundary( 1 );
    % string representing c coefficient matrix
    % c_schrod = c_coefficient();
    c_schrod = 'c_coefficient_schrod(sd)';
    % string representing a coefficient matrix, only include conduction
    % band potential, no polarization, dopant, static or exchange
    % potentials.
    a_schrod = 'V_conduction_band(x,y)';
    % d coefficient is constant for all electrons and mesh points, a scalar
    d_schrod = 1;
    % initial guess for the potential, evaluated at mesh points
    V_elestat = V_conduction_band(p_poiss(1, :), p_poiss(2, :))';
    V_poiss = zeros(size(V_elestat));
    %% 0. for the triangular geometry, we need to first compute the initial potential before solving schrodinger equation to include polarization charges
    if abs(degree_of_polygon) == 3
        f_poiss_1 = sprintf('n_D_func(x, y, 1) + sigma_spont(x, y) + sigma_piezo(x, y)');
        a_poiss = 0;
        b_poiss = generate_dirichlet_boundary( 1 );
        c_poiss = 'c_coefficient_poiss(sd)';
        V_poiss = assempde(b_poiss, p_poiss, e_poiss, t_poiss, c_poiss, a_poiss, f_poiss_1);
        V_poiss = V_poiss / C;
        V_elestat = V_elestat + V_poiss;
        a_schrod = sprintf('V_total_without_xc(x, y, %s)', mat2str(V_poiss));
    end
    
    %% debug purpose
%         x = p_poiss(1, :);
%         y = p_poiss(2, :);
%         debug_sigma = n_D_func(x, y, 1) + sigma_spont(x, y) + sigma_piezo(x, y);
%         
%         figure
%         pdesurf(10 * p_schrod, t_schrod, debug_sigma' );
%         colormap jet;
%         set(gcf, 'Renderer', 'zbuffer');
%         
%         figure
%         pdesurf(10 * p_schrod, t_schrod, V_poiss );
%         colormap jet;
%         set(gcf, 'Renderer', 'zbuffer');
    %% debug end
    
    %% variable declaration
    psi_schrodinger = []; % initialize the matrix psi_schrodinger
    E_schrodinger = []; % initialize the variable E_schrodinger
    epsilon_F = pinned_fermi_level;
    min_V = min(V_elestat, [], 'all');
    r_schrod = [min_V, min_V+10];
    
    %% 1. solve independent schrodinger equation, keep adding electrons until the charge neutrality constraint satisfied.
    while isempty(E_schrodinger) || E_schrodinger(end) <= epsilon_F
        [~, psi_schrodinger_add, E_schrodinger_add] = ...
            evalc('pdeeig(bmatrix, p_schrod, e_schrod, t_schrod, c_schrod, a_schrod, d_schrod, r_schrod)');
        
        % find new electrons with higher energies
        if ~isempty(E_schrodinger_add)
            
            psi_schrodinger = [psi_schrodinger, psi_schrodinger_add];
            E_schrodinger =   [E_schrodinger; E_schrodinger_add];
        end
        
        r_schrod(1) = r_schrod(2);
        r_schrod(2) = r_schrod(2) + 10;
        
    end
    
    % number of energy states below fermi level
    N = sum( epsilon_F >= E_schrodinger );
    %% debug purpose
    %N = 3;
    %epsilon_F = E_schrodinger(N);
    %% debug end
    
    % a row vector that calculates/contains the spatial regions of the
    % potential that are above the Fermi level; these spatial regions have
    % been ionized and contribute a positively-valued dopant density to the
    % Poisson equation
    heaviside_n_D = zeros( size(V_elestat') );
    heaviside_n_D( V_elestat > epsilon_F ) = 1;
    
    %heaviside_n_D = ones( size(V_elestat') );
    %% debug purpose, showing the wave functions, energies and fermi level
    if debug
       disp("Debug: initiate_system, the calculated eigen vectors, eigen-energies and fermi level are: ");
       psi_schrodinger
       E_schrodinger
       epsilon_F
    end
    
    %% 2. pre-compute some matrices, used later
    % Laplacian operator for possion equation
    b_poiss = generate_dirichlet_boundary( 1 );
    c_poiss = 'c_coefficient_poiss(sd)';
    thePde = pde.PDEModel(1);
    
    % poisson Boundary contributions
    [ ~, ~, H_poiss, R_poiss] = assemb(b_poiss, p_poiss, e_poiss);
    [N_poiss, O_poiss] = pdenullorth(H_poiss);

    % get laplacian operator of poisson equation
    [K_poiss, ~, ~] = thePde.assema(p_poiss, t_poiss, c_poiss, 0, 0);
    % find  K_poiss operator in the null space of H_poiss.
    K_poiss = N_poiss' * K_poiss * N_poiss;
    % If the K matrix is very close to symmetric, make it exactly symmetric
    % so \ will use the symmetric factorization routine.
    K_poiss = (K_poiss + K_poiss') / 2;
    
    % boundary matrices for schrodinger equation
    [Q, G, H, R] = assemb(bmatrix, p_schrod, e_schrod);
    % M matrice from 'd' coefficient
    [~, M_schrod] = thePde.assema(p_schrod, t_schrod, 0, d_schrod, 0);
    
    % gradient operators
    [DX, DY] = getgrad(p_poiss, t_poiss);
    [~, M] = thePde.assema(p_poiss, t_poiss, 0, 1, 0);
    DX = inv(M) * DX;
    DY = inv(M) * DY;
    
    %% 3. compute potential, including elestatic, and exchange potential
    V_poiss_old = V_poiss;
    [V_poiss] = solve_potential(N, epsilon_F, E_schrodinger, psi_schrodinger, heaviside_n_D, damping_factor, V_poiss);

    if abs(degree_of_polygon) == 3
        V_poiss = V_poiss_old + (V_poiss-V_poiss_old) * damping_factor;
    end
    
    %% debug purpose
%         figure
%         pdesurf(10 * p_schrod, t_schrod, V_poiss );
%         colormap jet;
%         set(gcf, 'Renderer', 'zbuffer');
%         
%         figure
%         pdesurf(10 * p_schrod, t_schrod, DX * V_poiss );
%         colormap jet;
%         set(gcf, 'Renderer', 'zbuffer');
%         
%         figure
%         pdesurf(10 * p_schrod, t_schrod, DY * V_poiss );
%         colormap jet;
%         set(gcf, 'Renderer', 'zbuffer');
%         
%         figure
%         pdesurf(10 * p_schrod, t_schrod, sigma_surface(p_schrod(1, :), p_schrod(2, :), V_poiss)' );
%         colormap jet;
%         set(gcf, 'Renderer', 'zbuffer');
    %%
    
    save(save_file, '-append');
    disp("Initialization finished.");
end