% function to compute
function [normalized_dis, density_dis] = compute_spatial_localization(psi_schrodinger, epsilon_F, E_schrodinger)
    
    global triangle_areas_schrod;
    global p_schrod;
    global t_schrod;
    global vector_of_side_lengths;
    global degree_of_polygon;
    global l0;
    %% input parameter check
    if (abs(degree_of_polygon) ~= 3) && (degree_of_polygon ~= 6) && (degree_of_polygon ~= 0)
       disp("Error: compute_spatial_localization, the geometry not supported!");
       return;
    end
    
    label_size = 15;
    title_size = 18;
    labelFormat = {'fontsize', label_size};
    titleFormat = {'fontsize', title_size};
    %%
    dis_corner = vector_of_side_lengths(end);
    if abs(degree_of_polygon) == 3
        dis_corner = sqrt(3) * dis_corner / 3;
    end
    
    num = size( psi_schrodinger, 2 );
    normalized_dis = zeros(1, num);
    
    % normalize the wavefunction first
    [normalization, ~] = normalize_and_sqrt_m_triangular( psi_schrodinger );
    
    % psi_schrodinger = psi_schrodinger .* normalization;
    psi_squared = psi_schrodinger .^ 2;
    
    % compute the distance r for each triangles
    x_A = p_schrod(1, t_schrod(1, :));
    x_B = p_schrod(1, t_schrod(2, :));
    x_C = p_schrod(1, t_schrod(3, :));
    y_A = p_schrod(2, t_schrod(1, :));
    y_B = p_schrod(2, t_schrod(2, :));
    y_C = p_schrod(2, t_schrod(3, :));
    x = (x_A + x_B + x_C) / 3;
    y = (y_A + y_B + y_C) / 3;
    triangle_distance = sqrt( x.^2 + y.^2 );
    
    %% compute averaged distance for each wavefunction
    for i = 1:num
        avg_psi_squared = 1/3 * normalization(i) .^ 2 * ...
            ( psi_squared( t_schrod(1, :), i ) + psi_squared( t_schrod(2, :), i ) + psi_squared( t_schrod(3, :), i ) );
        normalized_dis(i) = (triangle_areas_schrod .* triangle_distance) * avg_psi_squared / dis_corner;
    end
    
    fprintf('Computed normalized distances for each wavefunctions are:');
    normalized_dis
    
    %% compute averaged distance for electron density
    psi_sqrt_eps_sum = psi_sqrt_eps_summation( num, normalization, psi_schrodinger, epsilon_F, E_schrodinger );
    density = n_e_prefactor(p_schrod(1,:),p_schrod(2,:))' .* psi_sqrt_eps_sum;
    avg_density = 1/3 * ...
            ( density( t_schrod(1, :) ) + density( t_schrod(2, :) ) + density( t_schrod(3, :) ) );
    density_dis = triangle_areas_schrod .* triangle_distance * avg_density;
    density_dis = density_dis / ( triangle_areas_schrod * avg_density );
    density_dis = density_dis / dis_corner;
    
    fprintf('Computed normalized distance for total electron density: %f\n', density_dis);
    
    %% plot electron density
    figure;
    pdesurf(l0*p_schrod, t_schrod, density);
    view(0, 90);
    colormap jet;
    set(gcf, 'Renderer', 'zbuffer');
    xlabel('x (nm)', labelFormat{:});
    ylabel('y (nm)', labelFormat{:});
    title(sprintf('Total electron density, normalized distance: %4f\n', density_dis), titleFormat{:});
    axis equal;
    axis off;
        
    saveas(gcf, ['./data/plots/figures/electron_density.jpg']);

    %% compute area of x% of the height, for the total electron density
    % only available for the hexagon geometry
    % x: 90, 80, 70
    fprintf('total area is %f\n', sum(triangle_areas_schrod));
    ratios = [90 80 70];
    if degree_of_polygon == 6
        for i = 1 : length(ratios)
           thred = ratios(i);
           area = 0;
           for idx = 1 : 6
              heaviside_corner = select_corner(idx, false, false); 
              density_corner = avg_density' .* heaviside_corner;
              area = area + sum( (density_corner >= thred*max(density_corner)/100) .* triangle_areas_schrod );
              
               figure
               pdesurf(l0*p_schrod, t_schrod, density .* select_corner(idx, true, true)');
%                figure
%                pdesurf(l0*p_schrod, t_schrod, density .* select_corner(idx, true, false)');
              
           end
           fprintf('The area at %d%% maximum is: %f\n', thred, area/6);
        end
    end
    
    %% compute average radius and standard deviation for each peak
    % only available for the hexagon atm
    if degree_of_polygon == 6
        dis_average = 0;
        dis_std = 0;
        for idx = 1 : 6
            heaviside_corner_points = select_corner(idx, true, true); 
            density_corner_points = density' .* heaviside_corner_points;
            max_i = find(density_corner_points >= max(density_corner_points));
            x_max = p_schrod(1, max_i(1));
            y_max = p_schrod(2, max_i(1));
            
            heaviside_corner_area = select_corner(idx, false, true);
            density_corner_area = avg_density .* heaviside_corner_area';
            peak_distance = sqrt( (x-x_max).^2 + (y-y_max).^2 );
            % average distance to the max of the peak
            peak_dis_average = triangle_areas_schrod .* peak_distance * density_corner_area;
            peak_dis_average = peak_dis_average / ( triangle_areas_schrod * density_corner_area );
            % stand deviation
            peak_dis_std = triangle_areas_schrod .* (peak_distance.^2) * density_corner_area;
            peak_dis_std = peak_dis_std / ( triangle_areas_schrod * density_corner_area );
            
            dis_average = dis_average + peak_dis_average;
            dis_std = dis_std + sqrt(peak_dis_std);
        end
        fprintf('Averaged distance to the peak is: %f\n', dis_average/6);
        fprintf('Std of distance to the peak is: %f\n', dis_std/6);
    end
    
end


% corner selector
% idx: 1-6, which corner to choose
% ispoint: true/false, whether select points (true) or the area (false)
% narrow: true/false, whether choose more strict corner area
function heaviside_corner = select_corner(idx, ispoint, narrow)
    global p_schrod;
    global t_schrod;
    global vector_of_side_lengths;
    
    a = vector_of_side_lengths(1);
    if ispoint
        x = p_schrod(1, :);
        y = p_schrod(2, :);
    else
        x_A = p_schrod(1, t_schrod(1, :));
        x_B = p_schrod(1, t_schrod(2, :));
        x_C = p_schrod(1, t_schrod(3, :));
        y_A = p_schrod(2, t_schrod(1, :));
        y_B = p_schrod(2, t_schrod(2, :));
        y_C = p_schrod(2, t_schrod(3, :));
        x = (x_A + x_B + x_C) / 3;
        y = (y_A + y_B + y_C) / 3;
    end
    
    heaviside_corner = zeros( size(x) );
    
    switch idx
        case 1
            heaviside_corner = (x >= 0) .* (y > x/sqrt(3));
        case 2
            heaviside_corner = (y <= x/sqrt(3)) .* (y > -x/sqrt(3));
        case 3
            heaviside_corner = (y <= -x/sqrt(3)) .* (x > 0);
        case 4
            heaviside_corner = (x <= 0) .* (y < x/sqrt(3));
        case 5
            heaviside_corner = (y >= x/sqrt(3)) .* (y < -x/sqrt(3));
        case 6
            heaviside_corner = (y >= -x/sqrt(3)) .* (x < 0);
    end
    
    if narrow
        switch idx
            case 1
                heaviside_corner_narrow = (y >= -x/sqrt(3)+a/sqrt(3)) .* (y <= -x/sqrt(3)+2*a/sqrt(3));
            case 2
                heaviside_corner_narrow = (x >= a/2) .* (x <= a);
            case 3
                heaviside_corner_narrow = (y <= x/sqrt(3)-a/sqrt(3)) .* (y >= x/sqrt(3)-2*a/sqrt(3));
            case 4
                heaviside_corner_narrow = (y <= -x/sqrt(3)-a/sqrt(3)) .* ( y >= -x/sqrt(3)-2*a/sqrt(3));
            case 5
                heaviside_corner_narrow = (x <= -a/2) .* (x >= -a);
            case 6
                heaviside_corner_narrow = (y >= x/sqrt(3)+a/sqrt(3)) .* (y <= x/sqrt(3)+2*a/sqrt(3));
        end
        heaviside_corner = heaviside_corner .* heaviside_corner_narrow;
    end
end