function [collar_table_format, collar_length] = comprs_count(data_preproc, depth, len_spec, clr_len_def)

% input simplified signal of all decision CHs
data_input = data_preproc;
% len_spec = [39.2, 38, 39, 38.5, 40]; % from inner to outer
d_data = zeros(size(data_input));
dd_data = zeros(size(data_input));

table_depth = -999 * ones(size(data_input));
table_width = -999 * ones(size(data_input));

collar_factor = [1, 9/10, 8/10, 6/10, 5/10];
num_pre = 0;
for i_p = 1:size(data_input, 2)
    collar_num = floor(depth(end) / len_spec(i_p));
    
    table = zeros(size(data_input, 1), 3);
    
    % calculate derivation
    for i = 3:(size(data_input, 1)-2)
        d_data(i) = 1/12*(data_input(i-2, i_p) - 8*data_input(i-1, i_p) + 8*data_input(i+1, i_p) - data_input(i+2, i_p));
    end
    
    % calculate derivation
    for i = 3:(size(d_data, 1)-2)
        
        dd_data(i) = 1/12*(d_data(i-2) - 8*d_data(i-1) + 8*d_data(i+1) - d_data(i+2));
        
    end
    
    last_idx = 1;
    for i = 3:(size(d_data, 1)-2)
        if d_data(i-1) >=0 && d_data(i-2) >0 &&...
                d_data(i+1) <=0 && d_data(i+2) <0 && dd_data(i) < 0
            if i-1 == last_idx || i-2 == last_idx
                table(last_idx, 1:3) = [0, 0, 0];
            end
            left = -999;
            right = -999;
            table(i, 1) = depth(i);
            table(i, 2) = data_input(i, i_p);
            % find left side of the peak
            for j = i-2:-1:1
                if d_data(j)<0 && d_data(j-1)<0
                    idx_left = j;
                    left = depth(idx_left);
                    break;
                end
            end
            for k = i+2:size(d_data, 1)
                if d_data(k)<0 && d_data(k-1)<0
                    idx_right = k;
                    right = depth(idx_right);
                    break;
                end
            end
            if left ~= -999 && right~= -999
                table(i, 3) = abs(left-right);
                last_idx = i;
            else
                table(i, 1:3) = [0, 0, 0];
            end
        end
    end
    
    table_order = sortrows(table, -2);
    n = num_pre + floor(collar_num * collar_factor(i_p));
    num_pre = n;
    table_confined = sortrows(table_order(1:n, :), 1);
    table_confined = table_confined(table_confined(:, 1) ~= 0, :);

    table_depth(1:size(table_confined, 1), i_p) = table_confined(:, 1);
    table_width(1:size(table_confined, 1), i_p) = table_confined(:, 3);
end

collar_table = table_depth;
collar_table_format = zeros(size(data_input, 1), size(data_input, 2)+1);
width_table_format = zeros(size(data_input, 1), size(data_input, 2)+1);

% use longest collar table depth as depth and mark the second column as one
collar_table_format(:, 1) = collar_table(:, end);
column_temp = collar_table_format(:, 1);
column_temp = column_temp(column_temp > 0);
collar_table_format(1:size(column_temp, 1), end) = size(data_input, 2) * ones(size(column_temp));
width_table_format(1:size(column_temp, 1), end) = table_width(1:size(column_temp, 1), end);

flag = 0;
for m = 1:(size(data_input, 2)-1)
    % range of overlapping judgement
    local_width = table_width(table_width(:, m)>0, m);
    half_range = mean(local_width, 'all')/2;
    for n = 1:size(collar_table, 1)    
        if collar_table(n, m) > 0
            for i = 1:size(column_temp, 1)
                if flag == 1
                    flag = 0;
                    break;
                elseif collar_table(n, m) >= column_temp(i, 1) - half_range ...
                        && collar_table(n, m) <= column_temp(i, 1) + half_range
                    collar_table_format (i, m + 1) = m;
                    width_table_format (i, m + 1) = table_width(n, m);
                    % flag to end current collar
                    flag = 1;
                elseif collar_table(n, m) > column_temp(i, 1) + half_range
                    continue;
                elseif i>1 && collar_table(n, m) < column_temp(i, 1) - half_range ...
                        && collar_table(n, m) > column_temp(i-1, 1) + half_range
                    % compose a row containing new depth, 1 and 0, and insert
                    % to collar_table_format
                    row_temp = zeros(1, size(data_input, 2) + 1);
                    row_temp(1) = collar_table(n, m);
                    row_temp(m+1) = m;
                    index_temp = find(collar_table_format(:,1) > collar_table(n, m), 1);
                    collar_table_format = [collar_table_format(1:index_temp-1, :); ...
                        row_temp; collar_table_format(index_temp:end, :)];
                    row_temp_w = zeros(1, size(data_input, 2) + 1);
                    row_temp_w(m+1) = table_width(n, m);
                    width_table_format = [width_table_format(1:index_temp-1, :); ...
                        row_temp_w; width_table_format(index_temp:end, :)];
                    flag = 1;
                elseif i==1 && collar_table(n, m) < column_temp(i, 1) - half_range
                    row_temp = zeros(1, size(data_input, 2) + 1);
                    row_temp(1) = collar_table(n, m);
                    row_temp(m+1) = m;
                    index_temp = find(collar_table_format(:,1) > collar_table(n, m), 1);
                    collar_table_format = [row_temp; collar_table_format(1:end, :)];
                    row_temp_w = zeros(1, size(data_input, 2) + 1);
                    row_temp_w(m+1) = table_width(n, m);
                    width_table_format = [row_temp_w; width_table_format(index_temp:end, :)];
                    flag = 1;
                end
            end
        end
        width_table_format(:, 1) = collar_table_format(:, 1);
        column_temp = collar_table_format(:, 1);
        column_temp = column_temp(column_temp > 0);
    end
end

idx = find(collar_table_format(:, 1) < 0, 1);
collar_table_format = collar_table_format(1:idx-1, :);
width_table_format = width_table_format(1:idx-1, :);

% eliminate clamp caused error
if size(collar_table_format, 2) >= 5
    for n = 1:size(collar_table_format, 1)
        for m = 2:size(collar_table_format, 2)-1
            if m <= 3 && collar_table_format(n, m) > 0 ...
                    && collar_table_format(n, end) == 0 ...
                    && collar_table_format(n, end-1) == 0
                collar_table_format(n, m) = 0;
            elseif collar_table_format(n, m) > 0 ...
                    && collar_table_format(n, m+1) == 0
                collar_table_format(n, m+1) = m;
            end
        end
    end
else 
    for n = 1:size(collar_table_format, 1)
        for m = 2:size(collar_table_format, 2)-1
            if collar_table_format(n, m) > 0 ...
                    && collar_table_format(n, m+1) == 0
                collar_table_format(n, m+1) = m;
            end
        end
    end
end

collar_length = width_table_format(:, [1, end]);
collar_length(:, end) = collar_length(:, end)/2;
collar_length (collar_length == 0) = clr_len_def;

end

