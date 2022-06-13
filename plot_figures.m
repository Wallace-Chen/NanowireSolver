% function to plot checkpoints variables
function plot_figures(N, epsilon_F, psi_schrodinger, E_schrodinger, V_poiss, format)

    %% global variables
    global l0;
    global p_schrod p_poiss;
    global t_schrod;
    global C;
    
    global include_ex;
    global draw_ex;
    global M_schrod;
    global V_ex_old;
    
    global degree_of_polygon;
    global outer_shell_side_length;
    
    global fd out_folder;

    %% parameter check
    if abs(degree_of_polygon) ~= 3 && degree_of_polygon ~= 6 && degree_of_polygon ~= 0
       disp("Error: generateMesh, the geometry not supported!");
       return;
    end
    
    label_size = 15;
    title_size = 23;
    labelFormat = {'fontsize', label_size};
    titleFormat = {'fontsize', title_size};
    
    outf = [fd out_folder];
    %% plot all occupied wave-functions
    [normalization, ~] = normalize_and_sqrt_m_triangular( psi_schrodinger );
    
    for i = 1 : N
        figure
        % pdesurf(l0 * p_schrod, t_schrod, (normalization(i)*psi_schrodinger(:,i)).^2 );
        pdesurf(l0 * p_schrod, t_schrod, psi_schrodinger(:,i) );
        view(0, 90);
        colormap jet;
        set(gcf, 'Renderer', 'zbuffer');
        xlabel('x (nm)', labelFormat{:});
        ylabel('y (nm)', labelFormat{:});
        title(['E_{',num2str(i),'} = ',num2str(C*E_schrodinger(i)),' eV'], titleFormat{:});
        axis equal;
        axis off;
        
        saveas(gcf, [outf, '/figures/wavefunction_', num2str(i), format]);
    end
    
    for i = 1 : N
        figure
        pdesurf(l0 * p_schrod, t_schrod, (normalization(i)*psi_schrodinger(:,i)).^2 );
        view(0, 90);
        colormap jet;
        set(gcf, 'Renderer', 'zbuffer');
        xlabel('x (nm)', labelFormat{:});
        ylabel('y (nm)', labelFormat{:});
        title(['E_{',num2str(i),'} = ',num2str(C*E_schrodinger(i)),' eV'], titleFormat{:});
        axis equal;
        axis off;
        
        saveas(gcf, [outf, '/figures/wavefunction_', num2str(i), '_squared', format]);
    end
    
    %% plot the band-bending diagram and Fermi level
    figure
    x = linspace(0, 0, 200);
    if degree_of_polygon == 0
        x_min = -outer_shell_side_length;
        x_max = outer_shell_side_length;
    elseif abs(degree_of_polygon) == 3
        x_min = -sqrt(3)/6*outer_shell_side_length;
        x_max = sqrt(3)/3*outer_shell_side_length;
%     elseif degree_of_polygon == 3
%         x_max = sqrt(3)/6*outer_shell_side_length;
%         x_min = -sqrt(3)/3*outer_shell_side_length;
    elseif degree_of_polygon == 6
       x_min = -sqrt(3)/2*outer_shell_side_length;
       x_max = sqrt(3)/2*outer_shell_side_length;
    end
    y = linspace(x_min, x_max, 200);
    V_total_with_xc_val = V_total_without_xc(x, y, V_poiss);
    if include_ex && draw_ex
        V_ex_real = diag(full(inv(M_schrod) * V_ex_old));
        V_total_with_xc_val = V_total_with_xc_val + griddata_NaN(p_poiss(1,:),p_poiss(2,:),V_ex_real',x,y);
    end
    plot( l0*y, C*V_total_with_xc_val, 'LineWidth', 3 );
    hold on;
    plot( l0*[x_min, x_max], C*[epsilon_F epsilon_F], '--r', 'LineWidth', 3 );
    temp_axis = axis;
    %axis( [l0*x_min, l0*x_max 0 temp_axis(4)] );
    axis( [l0*x_min, l0*x_max temp_axis(3) temp_axis(4)] );
    xlabel('y (nm)', labelFormat{:});
    ylabel('Energy (eV)', labelFormat{:});
    title('Potential Bend Diagram', titleFormat{:});
    hold off
    
    fprintf('Fermi level is: %f eV');C*epsilon_F
    saveas(gcf, [outf, '/figures/potential_bend_diagram', format]);
    
    %% plot the total electron density
    figure
    psi_sqrt_eps_sum = psi_sqrt_eps_summation( N, normalization, psi_schrodinger, epsilon_F, E_schrodinger );
    pdesurf(l0*p_schrod, t_schrod, n_e_prefactor(p_schrod(1,:),p_schrod(2,:))' .* psi_sqrt_eps_sum);
    colormap jet
    lightangle(-37.5,30);
    set(gcf,'Renderer','zbuffer');
    grid on
    xlabel('x (nm)');
    ylabel('y (nm)');
    title('Total Electron Density', titleFormat{:});
    
    saveas(gcf, [outf, '/figures/total_electron_density', '.fig']);
    
    %% compute normalized distance
    compute_spatial_localization(psi_schrodinger(:, 1:N), epsilon_F, E_schrodinger(1:N));
    
    
end