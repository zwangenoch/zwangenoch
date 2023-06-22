function table_o=output(table_v, collar_length, length)
    table_o = zeros(size(table_v, 1), 8);
    table_o(:, 1) = (1:size(table_v, 1));
    
    for i = (1:size(table_v, 1))
        % find collar width
        [~, idx] = min(abs(table_v(i, 2) - collar_length(:, 1)));
        
        table_o(i, 5) = collar_length(idx, 2);
        table_o(i, 7) = table_v(i, 2) - table_o(i, 5)/2;
        table_o(i, 6) = table_v(i, 2) + table_o(i, 5)/2;
        table_o(i, 8) = table_v(i, 3);
        table_o(i, 3) = table_o(i, 6);
        
        % fill in previous bodybottom and bodylength
        if i > 1
            table_o(i-1, 2) = table_o(i, 7);
            table_o(i-1, 4) = table_o(i-1, 2) - table_o(i-1, 3);
        end
        % fill in for last depth
        if i == size(table_v, 1)
            table_o(i, 2) = table_o(i, 3) + length;
            table_o(i, 4) = length;    
        end
    end
end