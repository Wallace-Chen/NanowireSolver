% calculate area charge density due to spontaneous polarization at point
% (x,y), only for hexagonal geometry
function sigma = sigma_piezo(x, y)
    global degree_of_polygon;
    global e_cgs l0 eV;
    global pol_piezo;
    
    sigma = zeros( size(x) );
    if abs(degree_of_polygon) ~= 3
        return
    end
    
    factor = 4 * pi * e_cgs^2 / l0 / eV;
    factor = factor / eV * 1e-18 * l0^2;
    
    sigma = factor * triangle_edge_gauss(x, y, pol_piezo, true);
%     sigma_bottom = -piezo_analytic(0) * factor;
%     sigma_side = piezo_analytic(60 * pi / 180) * factor;
%     
%     if degree_of_polygon == 3 % a Ga-face triangle
%         sigma = triangle_edge_gauss(x, y, sigma_bottom, sigma_side);
%     elseif degree_of_polygon == -3 % a N-face triangle
%         sigma = triangle_edge_gauss(x, y, -sigma_bottom, -sigma_side);
%     end
    
end