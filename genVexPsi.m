% function to interpolate a function defined on schrodinger's meshgrid (p_schrod, t_schrod)
% fac indicate if we need to include the prefactor.
function f = genVexPsi(x, y, f, fac)

    global p_schrod t_schrod;
    global l0;
    global e_cgs;
    global eV;
    
    [x_sym, y_sym] = symmetrize_coordinates(x, y);
    %f = interpolant_NaN( p_schrod, t_schrod, f, x_sym, y_sym )';
    f = griddata_NaN(p_schrod(1,:), p_schrod(2,:), f', x, y);
    if fac > 0
        factor = 4 * pi * e_cgs^2 / l0 / eV;
        f = factor / 2 * n_e_prefactor(x, y) .* f;
    end

end