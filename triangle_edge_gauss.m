function val = triangle_edge_gauss(x, y, pols, is_piezo)

    global vector_of_side_lengths;
    global avg_side_length_poiss;
    global degree_of_polygon;
    
    vector_of_side_lengths_tmp = vector_of_side_lengths;
    a = 2 * avg_side_length_poiss;
    
    val = zeros(size(x));
    for i =  2 : length(vector_of_side_lengths_tmp)
        if is_piezo
            amp_bottom = eval(sprintf('(%s(0)-%s(0))', pols(i), pols(i-1)));
            amp_side = eval(sprintf('(%s(120*pi/180)-%s(120*pi/180))', pols(i), pols(i-1)));
        else
            amp_bottom = pols(i) - pols(i-1);
            amp_side = (pols(i) - pols(i-1)) * cosd(120);
        end
        
        if degree_of_polygon < 0 % a N-faced
            amp_bottom = -amp_bottom;
            amp_side = -amp_side;
        end
        
        l = vector_of_side_lengths_tmp(i);
        y_bottom =  -l * sqrt(3) / 6;
        
        val_bottom = zeros(size(x));
        val_left = zeros(size(x));
        val_right = zeros(size(x));
        
        indices = find( x >= -l/2 & x <= l/2 );
        val_bottom(indices) = amp_bottom / (a*sqrt(pi)) * exp( -(y(indices)-y_bottom).^2 / a^2 );
        
        % rotate 120 degree counter-clock wise
        rotated_points = rotate_matrix(120) * [x; y];
        indices = find( rotated_points(1,:) >= -l/2 & rotated_points(1,:) <= l/2 );
        val_left(indices) = amp_side / (a*sqrt(pi)) * exp( -(rotated_points(2,indices)-y_bottom).^2 / a^2 );
        
        % rotate 120 degree, clock wise      
        rotated_points = rotate_matrix(-120) * [x; y];
        indices = find( rotated_points(1,:) >= -l/2 & rotated_points(1,:) <= l/2 );
        val_right(indices) = amp_side / (a*sqrt(pi)) * exp( -(rotated_points(2,indices)-y_bottom).^2 / a^2 );
        
        val = val + (val_bottom + val_left + val_right);
    end
end

% a rotation matrix, by theta degree counter-clock wise
function matrix = rotate_matrix(theta)
    matrix = [ cosd(theta) -sind(theta); sind(theta) cosd(theta)];
end