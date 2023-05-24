clear;
close all; clc;

% input simplified signal of all decision CHs
simp_data_input = 1;
depth = 1;

data = simp_data_input;
% collar table with high confidence
collar_table_high_conf = -999 * ones(size(simp_data_input, 1), size(simp_data_input, 2));
% collar table with low confidence need to be determined later
collar_table_low_conf = -999 * ones(size(simp_data_input, 1), size(simp_data_input, 2));
collar_size = zeros(1, size(simp_data_input, 2));

%% First level - find collar location and amplitude based on non-zeros data
for CH = 1:size(simp_data_input, 2)
    idx_collar = 1; % track collar number
    % briefly eliminate lower energies, more will be done in data cleaning
    data(:, CH) = data(:, CH) - max(data(:, CH) / (8 * CH));
    top_depth = -999;
    bottom_depth = -999;
    top_array = -999 * ones(size(simp_data_input, 1), 1);
    bottom_array = -999 * ones(size(simp_data_input, 1), 1);
    deriv_1st = -999 * ones(size(simp_data_input, 1), 1);
    
    % first level data selection - find peaks
    for idx_depth = 1:size(simp_data_input, 1)
        if idx_depth < size(simp_data_input, 1) - 4 ...
                && idx_depth > 1 ...
                && data(idx_depth, CH) > 0 ...
                && data(idx_depth - 1, CH) < 0 ...
                && data(idx_depth + 1, CH) > 0 ...
                && data(idx_depth + 2, CH) > 0 ...
                && data(idx_depth + 3, CH) > 0 ...
                && bottom_depth == -999 % pair top and bottom
            top_depth = depth(idx_depth);
            top_array(idx_collar) = top_depth;
        elseif idx_depth > 4 ...
                && idx_depth < size(simp_data_input, 1) - 1 ...
                && data(idx_depth, CH) > 0 ...
                && data(idx_depth + 1, CH) < 0 ...
                && data(idx_depth - 1, CH) > 0 ...
                && data(idx_depth - 2, CH) > 0 ...
                && data(idx_depth - 3, CH) > 0 ...
                && top_depth ~= -999 % first should be top
            bottom_depth = depth(idx_depth);
            bottom_array(idx_collar) = bottom_depth;
            %collar_table(idx_collar, CH) = (top_depth + bottom_depth)/2;            
            idx_collar = idx_collar + 1;
            bottom_depth = -999;
        end    
    end
    
%% second level data selection - find overlaps
    idx_preselect = 1; % track subrange number
    idx_collar_high_conf = 1; % track collar number with high confidence
    idx_collar_low_conf = 1; % track collar number with low confidence
    idx_inner = 1; % index tracker within a subrange
    % calculate first order derivative of raw data with a 3 number gap
%     m = 1;
%     for n = 1:(size(simp_data_input, 1) - 3)
%         deriv_1st(m) = simp_data_input(n + 3) - simp_data_input(n);
%         m = m + 1;
%     end
    deriv_1st = diff(simp_data_input(:, CH));
    for idx_depth = 4:(size(deriv_1st, 1) - 3)
        if depth(idx_depth - 1) >= top_array(idx_preselect) ...
                && depth(idx_depth + 1) <= bottom_array(idx_preselect)
            x = depth(idx_depth); % for debug
            a=1; % for debug
            if deriv_1st(idx_depth) >= 0 ...
                    && deriv_1st(idx_depth - 1) >= 0 ...
                    && deriv_1st(idx_depth + 1) <= 0 ...
                    && deriv_1st(idx_depth - 1) + deriv_1st(idx_depth - 2) > 0 ...                    
                    && deriv_1st(idx_depth + 1) + deriv_1st(idx_depth + 2) < 0                 
                collar_table_high_conf(idx_collar_high_conf, CH) = depth(idx_depth);
                idx_collar_high_conf = idx_collar_high_conf + 1;
%                 && deriv_1st(idx_depth - 1) + deriv_1st(idx_depth - 2) + deriv_1st(idx_depth - 3) > 0 ...
%                     && deriv_1st(idx_depth - 2) + deriv_1st(idx_depth - 3) > 0 ...
%                 && deriv_1st(idx_depth + 1) + deriv_1st(idx_depth + 2) + deriv_1st(idx_depth + 3) < 0 ...
%                     && deriv_1st(idx_depth + 2) + deriv_1st(idx_depth + 3) < 0   
%             elseif deriv_1st(idx_depth - 1) >= 0 ...
%                     && deriv_1st(idx_depth + 1) <= 0 ...
%                     && deriv_1st(idx_depth - 1) >= deriv_1st(idx_depth)...
%                     && deriv_1st(idx_depth) >= deriv_1st(idx_depth + 1)...
%                 collar_table_low_conf(idx_collar_low_conf, CH) = depth(idx_depth);
%                 idx_collar_low_conf = idx_collar_low_conf + 1;                     
            end
        elseif depth(idx_depth) > bottom_array(idx_preselect)...
                && depth(idx_depth - 1) <= bottom_array(idx_preselect)
            idx_preselect = idx_preselect + 1;
        end     
    end
    
    % calculate average collar size (depth range of collar bump)
    size_diff_collar = bottom_array - top_array;
    collar_size(CH) = median(size_diff_collar(size_diff_collar>0));
       
    
end

%% collar table formatting
% target format:
%--------------------------------
% Depth   P5   P4   P3   P2   P1
% 1.1     1    0    0    0    0
% 2.3     1    1    1    0    0
% 4.5     1    1    1    1    1
%...
% 100.3   1    1    0    0    0
%--------------------------------

collar_table = collar_table_high_conf;
collar_table_format = zeros(size(simp_data_input, 1), size(simp_data_input, 2)+1);
% range of overlapping judgement
half_range = mean(collar_size)/2;

% use longest collar table depth as depth and mark the second column as one
collar_table_format(:, 1) = collar_table(:, end);
column_temp = collar_table_format(:, 1);
column_temp = column_temp(column_temp > 0);
collar_table_format(1:size(column_temp, 1), 2) = ones(size(column_temp));

for m = [(size(simp_data_input, 2)-1) : -1 :1]
    for n = 1:size(collar_table, 1)
        if collar_table(n, m) > 0
            for i = 1:size(column_temp, 1)
                if collar_table(n, m) >= column_temp(i, 1) - half_range ...
                        && collar_table(n, m) <= column_temp(i, 1) + half_range
                    collar_table_format(i, (size(simp_data_input, 2) - m + 2)) = 1;
%                     if i > 1
%                         column_temp(i, 1) = 0;
%                     else 
%                         column_temp(i, 1) = 0;
%                     end
                elseif collar_table(n, m) > column_temp(i, 1) + half_range
                    continue;
                elseif i>1 && collar_table(n, m) < column_temp(i, 1) - half_range ...
                        && collar_table(n, m) > column_temp(i-1, 1) + half_range
                    % compose a row containing new depth, 1 and 0, and insert
                    % to collar_table_format
                    row_temp = zeros(1, size(simp_data_input, 2) + 1);
                    row_temp(1) = collar_table(n, m);
                    row_temp(2:(size(simp_data_input, 2) - m + 1)) = ...
                        row_temp(2:(size(simp_data_input, 2) - m + 1)) + 1;
                    index_temp = find(collar_table_format(:,1) > column_temp(i-1, 1), 1);
                    collar_table_format = [collar_table_format(1:index_temp-1, :); ...
                        row_temp; collar_table_format(index_temp:end, :)];
                elseif i==1 && collar_table(n, m) < column_temp(i, 1) - half_range
                    row_temp = zeros(1, size(simp_data_input, 2) + 1);
                    row_temp(1) = collar_table(n, m);
                    row_temp(2:(size(simp_data_input, 2) - m + 1)) = ...
                        row_temp(2:(size(simp_data_input, 2) - m + 1)) + 1;
                    collar_table_format = [row_temp; collar_table_format(1:end, :)];
                end
            end
        end
        column_temp = collar_table_format(:, 1);
        column_temp = column_temp(column_temp > 0);
    end
end

idx = find(collar_table_format(:, 1) < 0, 1);
collar_table_format = collar_table_format(1:idx-1, :);

% eliminate clamp caused error
if size(collar_table_format, 2) >= 3
    for m = size(collar_table_format, 2) : -1 : 3
        for n = 1:size(collar_table_format, 1)
            if collar_table_format(n, m) == 1 && collar_table_format(n, m-1) == 0
                if size(collar_table_format, 2)>= 5 && m >= 5 ...
                    && collar_table_format(n, m-2) == 0 && collar_table_format(n, m-3) == 0
                        collar_table_format(n, m) = 0;
                else
                    collar_table_format(n, m-1) = 1;  
                end
            end
        end
    end
end
    
% time factor to all columns
for k = 2:size(collar_table_format, 2)
    collar_table_format(:, k) = collar_table_format(:, k) *...
        (size(collar_table_format, 2) - k + 1);
end
% reverse column order
temp = flip(collar_table_format(:, 2:end), 2);
collar_table_format(:, 2:end) = temp;

aaa=1;

