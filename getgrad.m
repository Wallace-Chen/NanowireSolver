% a function to calculate two gradient matrices given a finite element mesh.
% Note: these gradient matrices are in their weak forms, not in the real
% space, like K, M, F matrices inside the assempde function.
function [DX, DY] = getgrad(p, t)

    np = size(p, 2);

    % Corner point indices
    a1 = t(1,:);
    a2 = t(2,:);
    a3 = t(3,:);
    
    % Triangle sides
    r23x = p(1,a3) - p(1,a2);
    r23y = p(2,a3) - p(2,a2);
    r31x = p(1,a1) - p(1,a3);
    r31y = p(2,a1) - p(2,a3);
    r12x = p(1,a2) - p(1,a1);
    r12y = p(2,a2) - p(2,a1);
    
    % Area
    ar = abs( r31x .* r23y - r31y .* r23x ) / 2;
    
    % gradients of basis function
    g1x = 0.5 * r23y ./ ar;
    g1y = -0.5 * r23x ./ ar;
    g2x = 0.5 * r31y ./ ar;
    g2y = -0.5 * r31x ./ ar;
    g3x = 0.5 * r12y ./ ar;
    g3y = -0.5 * r12x ./ ar;
    
    % inner product between gradient directions and correct directions
    inner_prod1 = [r23y; -r23x] .* [r31x; r31y];
    inner_prod2 = [r31y; -r31x] .* [r12x; r12y];
    inner_prod3 = [r12y; -r12x] .* [r23x; r23y];
    
    % correct the sign of gradients
    g1x = g1x .* sign( sum(inner_prod1, 1) );
    g1y = g1y .* sign( sum(inner_prod1, 1) );
    g2x = g2x .* sign( sum(inner_prod2, 1) );
    g2y = g2y .* sign( sum(inner_prod2, 1) );
    g3x = g3x .* sign( sum(inner_prod3, 1) );
    g3y = g3y .* sign( sum(inner_prod3, 1) );
    
    % compute elements of gradients
    f1x = g1x .* ar / 3;
    f1y = g1y .* ar / 3;
    f2x = g2x .* ar / 3;
    f2y = g2y .* ar / 3;
    f3x = g3x .* ar / 3;
    f3y = g3y .* ar / 3;
    
    % assemble the matrices
    DX = sparse(a1, a1, f1x, np, np);
    DX = DX + sparse(a2, a1, f1x, np, np);
    DX = DX + sparse(a3, a1, f1x, np, np);
    DX = DX + sparse(a1, a2, f2x, np, np);
    DX = DX + sparse(a2, a2, f2x, np, np);
    DX = DX + sparse(a3, a2, f2x, np, np);
    DX = DX + sparse(a1, a3, f3x, np, np);
    DX = DX + sparse(a2, a3, f3x, np, np);
    DX = DX + sparse(a3, a3, f3x, np, np);
    
    DY = sparse(a1, a1, f1y, np, np);
    DY = DY + sparse(a2, a1, f1y, np, np);
    DY = DY + sparse(a3, a1, f1y, np, np);
    DY = DY + sparse(a1, a2, f2y, np, np);
    DY = DY + sparse(a2, a2, f2y, np, np);
    DY = DY + sparse(a3, a2, f2y, np, np);
    DY = DY + sparse(a1, a3, f3y, np, np);
    DY = DY + sparse(a2, a3, f3y, np, np);
    DY = DY + sparse(a3, a3, f3y, np, np);
    

end