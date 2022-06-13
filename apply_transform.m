function mat = apply_transform(mat)
    p = load('./order.mat');
    p = p.p;
    n = length(p);
    
    for i = 1:n
        [~,ind]=sort(p);
        to = ind(i);
        mat = transform(mat, i, to);
        p([i to]) = p([to i]);
    end

end

function mat = transform(mat, i, j)
    mat([i j], :) = mat([j i], :);
    mat(:, [i j]) = mat(:, [j i]);
end