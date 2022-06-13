% function to solve poisson and schrodinger equations iteratively until the
% potential gets converged.
% @param V_elestat: initial guess, contains only eletrcostatic potential.

% the fermi level is pinned
function [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N] = ...
    iterate_run_spectral_pinned(psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N)

    global debug;
    global p_poiss e_poiss t_poiss;
    global p_schrod e_schrod t_schrod;
    global include_ex;
    global draw_ex;
    global M_schrod;
    global V_ex_old;
    global pinned_fermi_level;
    
    global save_file;
    
    % build local potential
    V_elestat = V_poiss' + V_conduction_band(p_poiss(1, :), p_poiss(2, :));
    if include_ex && draw_ex
        V_ex_real = inv(M_schrod) * V_ex_old;
        V_elestat = V_elestat + diag(full(V_ex_real))';
    end
    V_elestat = V_elestat';
    
    step = 10; % search range
    %% section 1. solve schrodinger equation given the updated potentials
    fprintf('  1. Solving schrodinger equation...\n');
    
    % coefficient preparations for schordinger equations
    % c_schrod = c_coefficient();
    c_schrod = 'c_coefficient_schrod(sd)';
    % d coefficient is constant for all electrons and mesh points, a scalar
    %d_schrod = 1;
    % build the a coefficient matrix
    a_schrod = sprintf('V_total_without_xc(x, y, %s)', mat2str(V_poiss));
    %a_schrod = 'V_conduction_band(x,y)';
    % find the minimum of the whole potential matrix
    % mini = min(a_schrod, [], 'all');
    % fprintf('The minimum potential is: %f\n', mini);
    min_V =  min(V_elestat, [], 'all');
    
    r_schrod = [min_V, min_V+step];
    %r_schrod = [0 1]
    %b_schrod = generate_dirichlet_boundary( 1 );
    
    % before solving eigen equation, we assemble coefficient matrices first
    % which are the same within the iteration
    %b_poiss = generate_dirichlet_boundary( 1 );
    %c_poiss = 'c_coefficient_poiss(sd)';
    [K, M, B] = assem_pde(c_schrod, a_schrod, psi_schrodinger, E_schrodinger, epsilon_F, N);
    
    psi_schrodinger = []; % initialize the matrix psi_schrodinger
    E_schrodinger = []; % initialize the variable E_schrodinger
    epsilon_F = pinned_fermi_level;
    n_iter = 1;
    
    while isempty(E_schrodinger) || E_schrodinger(end) <= epsilon_F
        
        % solve for N coupled schrodinger equations
        
        [~, v_add, l_add, ires] = evalc("sptarn(K, M, r_schrod(1), r_schrod(2), 1)");
%         [~, v_add, l_add] = ...
%             evalc('pdeeig(b_schrod, p_schrod, e_schrod, t_schrod, c_schrod, a_schrod, d_schrod, r_schrod)');
        if ires < 0
            error(message('main:iterate_run_spectral:MoreEigenvaluesMayExist!'));
        end
        
        if ~isempty(v_add)
            v_add = B * v_add;
        end
        
        %% debug purpose
%            fprintf('    number of newly-found electrons is: %d\n', length(l_add));
%            fprintf('    number of found electrons is: %d\n', length(l));
        %%
            
        % add newly-found electrons
        if ~isempty(l_add)
            E_schrodinger = [E_schrodinger; l_add];
            psi_schrodinger = [psi_schrodinger, v_add];
        end
        % increase the search range
        r_schrod(1) = r_schrod(2);
        r_schrod(2) = r_schrod(2) + step;
        
        n_iter = n_iter + 1;
    end
    
    % number of energy states below fermi level
    N = sum( epsilon_F >= E_schrodinger );
    
    % a row vector that calculates/contains the spatial regions of the
    % potential that are above the Fermi level; these spatial regions have
    % been ionized and contribute a positively-valued dopant density to the
    % Poisson equation
    heaviside_n_D = zeros( size(V_elestat') );
    heaviside_n_D( V_elestat > epsilon_F ) = 1;
    
    %heaviside_n_D = ones( size(V_elestat') );
%     global degree_of_polygon;
%     global outer_shell_side_length;
%     global l0;
%     global C;
%     figure
%     x = linspace(0, 0, 200);
%     if degree_of_polygon == 0
%         x_min = -outer_shell_side_length;
%         x_max = outer_shell_side_length;
%     elseif abs(degree_of_polygon) == 3
%         x_min = -sqrt(3)/6*outer_shell_side_length;
%         x_max = sqrt(3)/3*outer_shell_side_length;
% %     elseif degree_of_polygon == 3
% %         x_max = sqrt(3)/6*outer_shell_side_length;
% %         x_min = -sqrt(3)/3*outer_shell_side_length;
%     elseif degree_of_polygon == 6
%        x_min = -sqrt(3)/2*outer_shell_side_length;
%        x_max = sqrt(3)/2*outer_shell_side_length;
%     end
%     y = linspace(x_min, x_max, 200);
%     V_total_with_xc_val = V_total_without_xc(x, y, V_poiss);
%     if include_ex && draw_ex
%         V_ex_real = diag(full(inv(M_schrod) * V_ex_old));
%         V_total_with_xc_val = V_total_with_xc_val - griddata_NaN(p_poiss(1,:),p_poiss(2,:),V_ex_real',x,y);
%     end
%     plot( l0*y, C*V_total_with_xc_val, 'LineWidth', 3 );
%     hold on;
%     plot( l0*[x_min, x_max], C*[epsilon_F epsilon_F], '--r', 'LineWidth', 3 );
%     temp_axis = axis;
%     %axis( [l0*x_min, l0*x_max 0 temp_axis(4)] );
%     axis( [l0*x_min, l0*x_max temp_axis(3) temp_axis(4)] );
    
    %% section 2. compute potential, including elestatic, and exchange potential
    fprintf('  2. Solving poisson equations...\n');
    
    [V_poiss] = solve_potential(N, epsilon_F, E_schrodinger, psi_schrodinger, heaviside_n_D, 1, V_poiss);
    
    save(save_file, '-append');
    
end