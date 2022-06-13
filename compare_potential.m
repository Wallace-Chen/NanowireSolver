function compare_potential(f1, f2)

    plot_exchange = true;

    global p_schrod t_schrod;
    global p_poiss t_poiss;
    global is_pinned;
    global vector_of_side_lengths;
    global vector_of_V0;
    global number_of_interfaces;
    global avg_side_length_schrod;
    global degree_of_polygon;
    
    global vector_of_masses;
    global triangle_indices_of_each_region_schrod;
    global triangle_areas_schrod;
    
    commons = load(f1, 'p_schrod', 't_schrod', 'p_poiss', 't_poiss', 'vector_of_side_lengths', 'vector_of_V0', 'number_of_interfaces', 'degree_of_polygon', 'avg_side_length_schrod', 'outer_shell_side_length', 'C', 'l0', 'epsilon_F', 'triangle_areas_schrod', 'vector_of_masses', 'triangle_indices_of_each_region_schrod');
    p_schrod = commons.p_schrod;
    t_schrod = commons.t_schrod;
    p_poiss = commons.p_poiss;
    t_poiss = commons.t_poiss;
    vector_of_side_lengths = commons.vector_of_side_lengths;
    vector_of_V0 = commons.vector_of_V0;
    number_of_interfaces = commons.number_of_interfaces;
    degree_of_polygon = commons.degree_of_polygon;
    avg_side_length_schrod = commons.avg_side_length_schrod;
    vector_of_masses = commons.vector_of_masses;
    triangle_indices_of_each_region_schrod = commons.triangle_indices_of_each_region_schrod;
    triangle_areas_schrod = commons.triangle_areas_schrod;
    outer_shell_side_length = commons.outer_shell_side_length;
    C = commons.C;
    l0 = commons.l0;
    epsilon_F = commons.epsilon_F;
    
    S1 = load(f1, 'V_poiss', 'M_schrod', 'V_ex_old', 'N', 'psi_schrodinger', 'epsilon_F', 'E_schrodinger');
    S2 = load(f2, 'V_poiss', 'M_schrod', 'V_ex_old', 'N', 'psi_schrodinger', 'epsilon_F', 'E_schrodinger');
    
    label_size = 15;
    title_size = 23;
    labelFormat = {'fontsize', label_size};
    titleFormat = {'fontsize', title_size};
    
    % compute electron density
    [~, psi_sqrt_m1] = normalize_and_sqrt_m_triangular(S1.psi_schrodinger);
    [~, psi_sqrt_m2] = normalize_and_sqrt_m_triangular(S2.psi_schrodinger);
    
    ele_den1 = 1 / pi * psi_sqrt_m1(1:S1.N) * sqrt( S1.epsilon_F-S1.E_schrodinger(1:S1.N) );
    ele_den2 = 1 / pi * psi_sqrt_m2(1:S2.N) * sqrt( S2.epsilon_F-S2.E_schrodinger(1:S2.N) );
    
    fprintf('fermi level 1: %f\n', S1.epsilon_F*C);
    fprintf('fermi level 2: %f\n', S2.epsilon_F*C);
    fprintf('total eletrons 1: %f\n', ele_den1);
    fprintf('total eletrons 2: %f\n', ele_den2);
    
    figure
    x = linspace(0, 0, 200);
    if degree_of_polygon == 0
        x_min = -outer_shell_side_length;
        x_max = outer_shell_side_length;
    elseif abs(degree_of_polygon) == 3
        x_min = -sqrt(3)/6*outer_shell_side_length;
        x_max = sqrt(3)/3*outer_shell_side_length;
    elseif degree_of_polygon == 6
       x_min = -sqrt(3)/2*outer_shell_side_length;
       x_max = sqrt(3)/2*outer_shell_side_length;
    end
    y = linspace(x_min, x_max, 200);
    V_total_val1 = V_total_without_xc(x, y, S1.V_poiss);
    V_total_val2 = V_total_without_xc(x, y, S2.V_poiss);
    if plot_exchange
        V_ex1 = diag(full(inv(S1.M_schrod) * S1.V_ex_old));
        V_total_val1 = V_total_val1 + griddata_NaN(p_poiss(1,:),p_poiss(2,:),V_ex1',x,y);
        V_ex2 = diag(full(inv(S2.M_schrod) * S2.V_ex_old));
        V_total_val2 = V_total_val2 + griddata_NaN(p_poiss(1,:),p_poiss(2,:),V_ex2',x,y);
    end
    plot( l0*y, C*V_total_val1, 'LineWidth', 3 );
    hold on;
    plot( l0*y, C*V_total_val2, 'LineWidth', 3 );
    hold on;
    plot( l0*[x_min, x_max], C*[epsilon_F epsilon_F], '--r', 'LineWidth', 3 );
    temp_axis = axis;
    axis( [l0*x_min, l0*x_max temp_axis(3) temp_axis(4)] );
    %axis( [l0*x_min, l0*x_max 0 0.1] );
    xlabel('y (nm)', labelFormat{:});
    ylabel('Energy (eV)', labelFormat{:});
    title('Potential Bend Diagram', titleFormat{:});
    hold off
    legend('No EX','With EX')
    
end