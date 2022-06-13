% function to solve eigen equation with non-local exchange potentials
% Input:
% b_schrod, c_schrod, a_schrod, d_schrod, r are standard inputs of Matlab's
% vanilla pdeeig. 
% Note only matrix type supported, PDEModel not implemented here.
% Note a_schrod should include all potentials except the exchange terms.
% b_poiss is the boundary matrix for poisson equation.
% psi, e, epsilon_F are wavefunctions, energies and fermi level from initial guess or
% previous iteration, which are used to construct the exchange terms. N is
% the number of occupied energy states.
% Output
% K, M are coefficient matrices of the generalized eigen equation:
% KU = \lambda MU
% B is the null space operator of H.
% Attention: the script implicitly assumes the Dirichlet boundary conditions
% for poisson equations. Techniquely, Q, G matrix are taken to be zeros.
% No assumptions made for schrodinger equation.
function [K, M, B] = assem_pde(c_schrod, a_schrod, psi, e, epsilon_F, N)

    global p_poiss e_poiss t_poiss;
    global p_schrod e_schrod t_schrod;
    global C;
    global include_ex;
    global V_ex_old;
    global K_poiss N_poiss;
    global Q G H R M_schrod;
    
    global save_file;
    global damping_factor;
    
    np = size(p_poiss, 2);
    V_ex = sparse(np, np);
    thePde = pde.PDEModel(1);
    if include_ex
%         % poisson Boundary contributions
%          [ ~, ~, H_poiss, R_poiss] = assemb(b_poiss, p_poiss, e_poiss);
%          [N_poiss, O_poiss] = pdenullorth(H_poiss);
% 
%         % get laplacian operator of poisson equation
%         [K_poiss, ~, ~] = thePde.assema(p_poiss, t_poiss, c_poiss, 0, 0);
%         % find  K_poiss operator in the null space of H_poiss.
%         K_poiss = N_poiss' * K_poiss * N_poiss;
%         % If the K matrix is very close to symmetric, make it exactly symmetric
%         % so \ will use the symmetric factorization routine.
%         K_poiss = (K_poiss + K_poiss') / 2;

        % build the exchange potentials
        [normalization, ~] = normalize_and_sqrt_m_triangular(psi);
        psi_normalized = psi .* normalization;
        for i = 1 : N
            psi_i = psi_normalized(:, i);
            [~, psi_star_expd, ~] = thePde.assema(p_poiss, t_poiss, 0, sprintf('genVexPsi(x, y, %s, %d)', mat2str(psi_i), 1), 0);
            [~, psi_expd, ~] = thePde.assema(p_poiss, t_poiss, 0, sprintf('genVexPsi(x, y, %s, %d)', mat2str(psi_i), 0), 0);
            % map to null space
            psi_star_expd = N_poiss' * psi_star_expd;
            psi_expd = N_poiss' * psi_expd * N_poiss;
            % compute exchange term, and restore the full space
            F = N_poiss * (psi_expd * (K_poiss \ psi_star_expd));
            V_ex = V_ex + F * sqrt( epsilon_F -  e(i)); 
        end
        
        V_ex = V_ex / C; % do the unit conversion
        
        V_ex = V_ex_old + (V_ex - V_ex_old) * damping_factor;
        V_ex_old = V_ex;
    end
    
    % schrodinger boundary conditions
%     [Q, G, H, R] = assemb(b_schrod, p_schrod, e_schrod);
    [K_, M_ele, F] = thePde.assema(p_schrod, t_schrod, c_schrod, a_schrod, 0);
    M_tot = M_ele + V_ex; % sign is considered already in psi_star_expd
    % cannot set zero for M_tot matrix, only zero your potential via a_schrod
    % coefficient. Since we are now in the finite element representation,
    % not in the real space, an all one's matrix here does NOT correspond to
    % a constant potential in real space.
    %M_tot = M_tot - min(M_tot, [], 'all');
    [K, ~, B, ~] = assempde(K_, M_tot, F, Q, G, H, R);
%     [~, M] = thePde.assema(p_schrod, t_schrod, 0, d_schrod, 0);
    M = B' * M_schrod * B;
    
    %% debug end
    
    % finally, solve the generalized eigen equation by Alnordi algorithm
    % [v,l,ires] = sptarn(K, M, r(1), r(2), 0);
    % v=B*v;
    %% debug purpose
    %[~, M_debug_poiss, ~] = thePde.assema(p_schrod, t_schrod, 0, 'helpler(x,y)', 0);
    %% debug end
    
    save(save_file, '-append');
end