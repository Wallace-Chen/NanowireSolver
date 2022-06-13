%--------------------------------------------------------------------------
% main function to calculate eigen wavefunctions of the first fewv lowest eigen
% energies for nanowires with various geometries, including the non-local
% exchange potentials between electrons.
%--------------------------------------------------------------------------

%function main()
function main(i)

%     clear all;
%     clear global;
    tic;
    disp("Welcome to the PROGRAM solving the schrodinger equations for the given nanowires geometry.");
    disp("Loading variables and generating geometry...");
    %% global variable declaration
    % for debug purpose in development phase, set to false in production
    % mode 
    global debug;
    debug = false;
    
    global avg_side_length_schrod;
    global avg_side_length_poiss;
    global p_poiss; global e_poiss; global t_poiss;
    global p_schrod; global e_schrod; global t_schrod;
    global number_mesh_points_schrod;
    global number_mesh_points_poiss;
    global triangle_indices_of_each_region_schrod;
    global triangle_indices_of_each_region_poiss;
    global triangle_areas_schrod;
    global triangle_areas_poiss;
    global damping_factor;
    global include_ex;
    global draw_ex;
    global V_ex_old;
    global save_file;
    global fd;
    
    global ckp_n;
    global ckp_average_change;
    global ckp_epsilon_F;
    global ckp_E_schrodinger;
    global ckp_E_schrodinger_debug;
    global ckp_psi_schrodinger;
    
    %% metric variables declaration
    ckp_n = []; % array to hold number of eletrons below fermi level for each iterations
    ckp_average_change = []; % array to hold average change of potentials (exclude exchange potential) for each iterations
    ckp_epsilon_F = []; % array to hold fermi levels for each iterations
    ckp_E_schrodinger = {}; % cell array to hold electrons' energies for each iterations
    ckp_E_schrodinger_debug = {}; % cell array to hold all electrons' energies for each iterations
    ckp_psi_schrodinger = {}; % cell array to hold electrons' wavefunctions for each iterations
    %% Preparation, load constants and variables

    % load constants from the file
    run('./constant.m'); 

    % load user-defined variables for the configurations considered in the
    % problem 
    %run('./input_parameters.m');
    run(sprintf('./inputFiles/input_parameters_%d.m', i));
    out_folder = sprintf('%d_%s', i, out_folder);
    
    % create output folder
    fd = './data/';
    mkdir(fd, out_folder);
    
    % the name of save file
    name_of_save_file = sprintf('%d_charge_neutrality', degree_of_polygon);
    for i = 1 : number_of_interfaces
        name_of_save_file = [name_of_save_file, '_', num2str(vector_of_side_lengths(i))];
    end
    if include_ex
        save_file = ['./data/', out_folder, '/', name_of_save_file, '_ex', '.mat'];
    else
        save_file = ['./data/', out_folder, '/', name_of_save_file, '_noex', '.mat'];
    end
    
    save(save_file);
    %% preprocessing of variables
    
    if abs(degree_of_polygon) ~= 3 && degree_of_polygon ~= 6 && degree_of_polygon ~= 0
       disp("Error: main function, the geometry not supported!");
       return;
    end
    
    if include_ex
       fprintf('Exchange potential will be considered in the computation.\n');
       if draw_ex
           fprintf('Exchange potential will be plot and included for fermi level, depletion region calculation.\n');
       end
    else
       fprintf('Note: Exchange potential will NOT be considered in the computation.\n');
    end
    
    
    pinned_fermi_level = pinned_fermi_level / C;
    
    if abs(degree_of_polygon) == 3
        % set the potential energy of the outer shell edge to be zero
        vector_of_V0 = (vector_of_V0 - vector_of_V0(1)) / C;
    else
        % set the minimum of the bandedge energy to be zero and re-scale the
        % bandedge energy in dimensionless units, for the charge-neutrality
        % case
        vector_of_V0 = (vector_of_V0 - min(vector_of_V0)) / C;
    end
    
    % generate the mesh points and data to be used by the pdeeig program in the
    % MATLAB PDE Toolbox for solving the Schrodinger equation
    [p_schrod, e_schrod, t_schrod, triangle_indices_of_each_region_schrod, triangle_areas_schrod] = ...
        generate_mesh();

    % use the same mesh points and data to be used by the assempde program in
    % the MATLAB PDE Toolbox for solving the Poisson equation (one can modify
    % this if necessary)
    p_poiss = p_schrod; e_poiss = e_schrod; t_poiss = t_schrod;
    number_mesh_points_schrod = length(p_schrod(1, :));
    number_mesh_points_poiss = length(p_poiss(1, :));
    triangle_indices_of_each_region_poiss = triangle_indices_of_each_region_schrod;
    triangle_areas_poiss = triangle_areas_schrod;

    % compute the average side length of the triangular mesh
    avg_side_length_schrod = average_side_length(p_schrod, t_schrod);
    avg_side_length_poiss = average_side_length(p_poiss, t_poiss);
    V_ex_old = sparse(size(p_poiss,2), size(p_poiss,2));
    %% initializing the system
    damping_factor = 0.1;
    
%     if abs(degree_of_polygon) == 3
%         [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N] = initiate_system_pinned();
%     else
%         [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N] = initiate_system();
%     end
    [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N] = initiate_system_new();
    average_change = 1 / C;
    N_old = 0;
    n_iter = 1;
    
    ckp_n = [N];
    ckp_average_change = [0];
    ckp_epsilon_F = [epsilon_F];
    ckp_E_schrodinger{end+1} = E_schrodinger(1:N);
    ckp_E_schrodinger_debug{end+1} = E_schrodinger;
    ckp_psi_schrodinger{end+1} = psi_schrodinger;
    %% solve problem iteratively
    good_iters = 0;
    % require at least 10 interations, apply convergence criteria
    while n_iter <= 10 || good_iters < 3

        fprintf('\n ----- the # of iteration: %d ----- \n', n_iter);
        fprintf('  Input: %d electrons, the fermi enery is: %f, the average change: %f the computed energies are: \n', N, epsilon_F, average_change);
        E_schrodinger
        
        N_old = N; V_poiss_old = V_poiss;
        
%         if abs(degree_of_polygon) == 3
%             [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N] = ...
%                 iterate_run_spectral_pinned(psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N);
%         else
%             [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N] = ...
%                 iterate_run_spectral(psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N);
%         end
        [psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N] = ...
                iterate_run_new(psi_schrodinger, E_schrodinger, epsilon_F, V_poiss, N);
        
        
        average_change = mean( abs(V_poiss-V_poiss_old) );
        
        % store checkpoints
        ckp_n = [ckp_n, N];
        ckp_average_change = [ckp_average_change, average_change];
        ckp_epsilon_F = [ckp_epsilon_F, epsilon_F];
        ckp_E_schrodinger{end+1} = E_schrodinger(1:N);
        ckp_E_schrodinger_debug{end+1} = E_schrodinger;
        ckp_psi_schrodinger{end+1} = psi_schrodinger;
        
        % compute damping factor
        if average_change < 70
            damping_factor = min( [damping_factor + 0.01 0.1] );
        else
            damping_factor = 0.05;
        end
        
        fprintf(['Potential change is: ',num2str(average_change*C), '\n']);
        fprintf(['Potential change percent is: ', num2str(average_change/mean(abs(V_poiss))), '\n']);
        fprintf(['Damping factor is: ', num2str(damping_factor), '\n']);
        
        if average_change / mean(abs(V_poiss)) > 0.1 || average_change > 0.01/C || N ~= N_old
            good_iters = 0;
        else
            good_iters = good_iters + 1;
        end
        
        if n_iter <= 10 || good_iters < 3
            % update potentials
            V_poiss = V_poiss_old + damping_factor * (V_poiss - V_poiss_old);
        end
        
        n_iter = n_iter + 1;
        
        %% debug purpose, force stop
        if n_iter > 50
            break
        end
        
    end
    toc;
    %%
    fprintf('Iteration finished. \n');
    fprintf('%d electrons, the fermi enery is: %f, the average change: %f \n', N, epsilon_F, average_change);
    %% save variables in .mat file
    fprintf('saving variables to the file...\n');
    save(save_file, '-append');
    fprintf('\n Variables saved to the file %s. \n', save_file);
    %% plot
    fprintf('ploting figures...\n');
    if include_ex
        plot_figures(N, epsilon_F, psi_schrodinger, E_schrodinger, V_poiss, '_ex.jpg');
        fprintf('plotting checkpoints figures...\n');
        plot_checkpoints('_ex');
    else
        plot_figures(N, epsilon_F, psi_schrodinger, E_schrodinger, V_poiss, '_noex.jpg');
        fprintf('plotting checkpoints figures...\n');
        plot_checkpoints('_noex');
    end
    fprintf('all done.\n');
end