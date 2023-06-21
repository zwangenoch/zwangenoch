function table_v = verify(table, length, margin)

% this function is to verify the table input
% this function only do minor revision based on previous 
% decisions

if size(table, 1) <= 2
    table_v = table; 
else
    table_v = zeros(size(table));
    table_v(1, :) = table(1, :);
    
    n = 2;
    
    for i = 2:(size(table, 1)-1)
        if abs(table(i, 2) - table_v(n-1, 2)) >= (length - margin)...
                && abs(table(i, 2) - table_v(n-1, 2)) <= (length + margin)
            table_v(n, 2) = table(i, 2);
            table_v(n, 3) = 1;
            n = n + 1;
        elseif table(i, 2) == table_v(n-1, 2) ...
                || abs(table(i, 2) - table_v(n-1, 2)) <= margin
            table_v(n, 3) = min([table_v(n-1, 3), table(i, 3)]);
        elseif abs(table(i, 2) - table_v(n-1, 2)) >= (2*(length - margin))
            table_v (n, 2) = table_v(n-1, 2) + length;
            table_v (n, 3) = 3;     
            n = n + 1;
        else
            table_v(n, 2) = table(i, 2);
            table_v(n, 3) = table(i, 3);
            n = n + 1;
        end
    end
end

idx_end = find ((table_v(:, 2) == 0), 1);
table_v = table_v (1:idx_end-1, :);
table_v(:, 1) = [(1:size(table_v, 1))];

end