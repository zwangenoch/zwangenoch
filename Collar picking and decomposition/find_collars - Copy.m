clear;
close all; clc;

% input simplified signal of all decision CHs
simp_data_input = 1;
depth = 1;

%% First approach - find collar location and amplitude based on non-zeros data
data = simp_data_input;
% collar table with high confidence
collar_table_high_conf = -999 * ones(size(simp_data_input, 1), size(simp_data_input, 2));
% collar table with low confidence need to be determined later
collar_table_low_conf = -999 * ones(size(simp_data_input, 1), size(simp_data_input, 2));

for CH = 1:size(simp_data_input, 2)
    idx_collar = 1; % track collar number
    % briefly eliminate lower energies, more will be done in data cleaning
    data(:, CH) = data(:, CH) - max(data(:, CH) / (5 * CH));
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
    
    % second level data selection - find overlaps
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
            if deriv_1st(idx_depth) > 0 ...
                    && deriv_1st(idx_depth - 1) > 0 ...
                    && deriv_1st(idx_depth + 1) < 0 ...
                    && deriv_1st(idx_depth - 1) + deriv_1st(idx_depth - 2) > 0 ...
                    && deriv_1st(idx_depth - 1) + deriv_1st(idx_depth - 2) + deriv_1st(idx_depth - 3) > 0 ...
                    && deriv_1st(idx_depth - 2) + deriv_1st(idx_depth - 3) > 0 ...
                    && deriv_1st(idx_depth + 1) + deriv_1st(idx_depth + 2) < 0 ...
                    && deriv_1st(idx_depth + 1) + deriv_1st(idx_depth + 2) + deriv_1st(idx_depth + 3) < 0 ...
                    && deriv_1st(idx_depth + 2) + deriv_1st(idx_depth + 3) < 0   
                collar_table_high_conf(idx_collar_high_conf, CH) = depth(idx_depth);
                idx_collar_high_conf = idx_collar_high_conf + 1;
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
       
    
end


aaa=1;

