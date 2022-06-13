function debug_exchange_term(psi)

    global p_poiss t_poiss;
    global lhs rhs;
    N = size(psi, 2);
    np = size(psi, 1);
    change_ratio = zeros(N, N);
    thePde = pde.PDEModel(1);
    num = 1;
    for i = 1 : N
        f = psi(:, i);
        f = ones(np, 1);
        [~, lhs0, ~] = thePde.assema(p_poiss, t_poiss, 0, sprintf('genVexPsi(x, y, %s, %d)', mat2str(f), 0), 0);
        for j = 1 : N
            if num > 1
                return
            end
            fprintf('--- working on (%d, %d) ---\n', i, j);
            v = psi(:, j);
            %v = 3*ones(np, 1);
            lhs = lhs0 * v;
            rhs = f .* v;
            [~, ~, rhs] = thePde.assema(p_poiss, t_poiss, 0, 0, sprintf('genVexPsi(x, y, %s, %d)', mat2str(rhs), 0));
            
            abs(lhs - rhs) ./ abs(rhs);
            change_ratio(i, j) = sum(abs(lhs - rhs)) / sum(abs(rhs));
            
            num = num + 1;
        end
    end
    
    change_ratio

end