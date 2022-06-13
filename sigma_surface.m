% a function to calculate the surface charge on the outer shell given a
% potential energy (experienced by an electron), by applying the boundary
% conditions assuming the charge neutrality and equi-potential outer shell.

function sigma = sigma_surface(x, y, V)
    global DX DY;
    global degree_of_polygon;
    global p_poiss e_poiss;
    global vector_of_eps;
    global C;
    
    sigma = zeros(size(x));
    
    % only for triangular geometry atm
    if abs(degree_of_polygon) ~= 3
        return
    end
    % convert back to eV / 10nm
    V = V * C;
    
    sigma = zeros( size(V) );
    
    EX = DX * V;
    EY = DY * V;
    
    % extract boundary points
    p1 = find(e_poiss(5, :) == 1);
    p1 = union( e_poiss(1, p1),  e_poiss(2, p1) );
    p2 = find(e_poiss(5, :) == 2);
    p2 = union( e_poiss(1, p2),  e_poiss(2, p2) );
    p3 = find(e_poiss(5, :) == 3);
    p3 = union( e_poiss(1, p3),  e_poiss(2, p3) );
    
    % surface charge
    sigma(p1) = EX(p1)*sqrt(3)/2 - EY(p1)/2;
    sigma(p2) = EY(p2);
    sigma(p3) = -EX(p3)*sqrt(3)/2 - EY(p3)/2;
    
    sigma = vector_of_eps(1) * griddata_NaN( p_poiss(1, :), p_poiss(2, :), sigma, x, y );
    

end