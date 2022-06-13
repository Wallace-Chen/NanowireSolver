function extract_order(matrix)

    p = diag(matrix);
%    [p,ind]=sort(p);
    save('./order.mat', 'p');
end