function [table_P1, table_P2, table_P3, table_P4, table_P5]...
    = decisionmaking(collar_info, depth, len_spec, margin, anchor_depth)

table_P1 = 0;
table_P2 = 0;
table_P3 = 0;
table_P4 = 0;
table_P5 = 0;

%% initialization
% input specified pipe length-length between joints, and margin, in feet
% len_spec = [39.2, 38, 39, 38.5, 40]; % from inner to outer
% margin = 3;
%anchor_depth = [76.199, 105.699, 90.116529, 385.282015, 830.3636];
% anchor_depth = [77.599, 105.699, 90.116529, 385.282015, 830.3636];

% browse according to pipe number
for p_No = 1:size(collar_info, 2)
    % assign length
    length = len_spec(p_No);
    
    [~, anchor_idx] = min(abs(depth - anchor_depth(p_No)));
%     anchor_idx = find(round(depth) == round(anchor_depth(p_No)), 1, 'first');
    % generate collar table filled with 0s
    table = zeros(1200, 3);
    table(:, 1) = (1:1200)';
    table(:, 2) = -999 * ones(1200, 1);
    table(1, 2) = depth(anchor_idx);
    table(1, 3) = 1;
    i_table = 2;
    % generate reverse collar table filled with 0s from anchor to top
    table_r = zeros(1200, 3);
    table_r(:, 1) = (1:1200)';
    table_r(:, 2) = -999 * ones(1200, 1);
    table_r(1, 2) = depth(anchor_idx);
    table_r(1, 3) = 1;
    i_table_r = 2;
        
    % Pipe 1-2 and not outer most pipe: browse all possible collars along depth   
    if p_No < 3 && p_No < size(collar_info, 2)
        % make decision from anchor to bottom
        i = anchor_idx;        
        while i <= size(collar_info, 1) && i>= 1
            flag_repeat = 0;
            if collar_info(i, p_No) == p_No
                if abs(depth(i) - table(i_table - 1, 2)) <= (length+1*margin) && ...
                        abs(depth(i) - table(i_table - 1, 2)) >= (length-1*margin)
                    if collar_info(i, end)~=0 % not clamp
                        table(i_table, 2) = depth(i);
                        table(i_table, 3) = 1;
                        collar_info(i, :) = collar_info(i, :) - 1;                        
                        % if next is high probablity to be green, the previous
                        % that match length rule is green as well
                        if table(i_table-1, 3) == 3 
                            table(i_table-1, 3) = 2;
                        elseif table(i_table-1, 3) == 2
                            table(i_table-1, 3) = 1;
                        end
                        i_table = i_table + 1;
                    else
                        i = i + 1;
                        continue % is clamp
                    end
                    
                else
                    if abs(depth(i) - table(i_table - 1, 2)) < (length-1*margin)
                        if abs(depth(i) - table(i_table - 1, 2)) <= 2* margin
                            if collar_info(i, end)~=0 % not clamp and replace last collar
                                i_restore = find(round(depth) == round(table(i_table - 1, 2)));
                                collar_info(i_restore, :) = collar_info(i_restore, :) + 1;
                                table(i_table - 1, 2) = depth(i);
                                table(i_table - 1, 3) = 2;
                                collar_info(i, :) = collar_info(i, :) - 1;
                            else
                                i = i + 1;
                                continue % is clamp
                            end
                        elseif abs(depth(i) - table(i_table - 1, 2)) > 2* margin    
                            if table(i_table - 1, 3) ~= 1
                                if i_table > 2 && calc_mode(depth(i),table(i_table - 2, 2), length)...
                                        <= calc_mode(table(i_table - 1, 2),table(i_table - 2, 2), length)%
                                    if collar_info(i, end)~=0 % not clamp and replace last collar
                                        i_restore = find(round(depth) == round(table(i_table - 1, 2)));
                                        collar_info(i_restore, :) = collar_info(i_restore, :) + 1;
                                        table(i_table - 1, 2) = depth(i);
                                        table(i_table - 1, 3) = 2;
                                        collar_info(i, :) = collar_info(i, :) - 1;
                                    else
                                        i = i + 1;
                                        continue % is clamp
                                    end
                                elseif i_table <= 2
                                    if collar_info(i, end)~=0 % not clamp and replace last collar
                                        table(i_table, 2) = depth(i);
                                        table(i_table, 3) = 2;
                                        collar_info(i, :) = collar_info(i, :) - 1;
                                        i_table = i_table + 1;
                                    else
                                        i = i + 1;
                                        continue % is clamp
                                    end
                                else
                                    i = i + 1;
                                    continue
                                end
                            elseif table(i_table - 1, 3) == 1
                                if collar_info(i, end)~=0 % not clamp 
                                    table(i_table, 2) = depth(i);
                                    table(i_table, 3) = 2;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                    i_table = i_table + 1;
                                else
                                    i = i + 1;
                                    continue % is clamp
                                end                                
                            end
                        end
                    elseif abs(depth(i) - table(i_table - 1, 2)) > (length+1*margin)
                        if abs(depth(i) - table(i_table - 1, 2)) <= (1.2 * length + margin)
                            if collar_info(i, end)~=0 % not clamp
                                table(i_table, 2) = depth(i);
                                table(i_table, 3) = 2;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table = i_table + 1;
                            else
                                i = i + 1;
                                continue % is clamp
                            end 
                        elseif abs(depth(i) - table(i_table - 1, 2)) > (1.2 * length + margin)
                            drct = 1;
                            [num_inrange, depth_1st, depth_2nd] = ...
                                count_num(depth, table, i, p_No, collar_info, length, margin,...
                                table(i_table - 1, 2), drct);
                            if num_inrange > 0
                                % find index of depth_1st and depth_2nd
                                idx_depth_1st = find(round(depth) == round(depth_1st));
                                if depth_2nd ~= -99999
                                    idx_depth_2nd = find(round(depth) == round(depth_2nd));
                                else
                                    idx_depth_2nd = NaN;
                                end
                                if num_inrange == 1
                                    if collar_info(idx_depth_1st, end)~=0 % not clamp
                                        if collar_info(i, end)~=0 % not clamp
                                            table(i_table, 2) = depth_1st;
                                            table(i_table, 3) = 3;                     
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table = i_table + 1;
                                            flag_repeat = 1;
                                        else 
                                            table(i_table, 2) = depth_1st;
                                            table(i_table, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table = i_table + 1;
                                        end
                                    elseif collar_info(idx_depth_1st, end) ==0 % clamp
                                        if collar_info(i, end)~=0 % not clamp
                                            table(i_table, 2) = depth(i);
                                            table(i_table, 3) = 3;
                                            collar_info(i, :) = collar_info(i, :) - 1;
                                            i_table = i_table + 1;
                                        else
                                            i = i + 1;
                                            continue % is clamp
                                        end
                                    end             
                                elseif num_inrange > 1
                                    if collar_info(idx_depth_1st, end)~=0 % not clamp
                                        if collar_info(i, end)~=0 % not clamp
                                            table(i_table, 2) = depth_1st;
                                            table(i_table, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table = i_table + 1;
                                            flag_repeat = 1;
                                        else 
                                            table(i_table, 2) = depth_1st;
                                            table(i_table, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table = i_table + 1;
                                        end
                                    elseif collar_info(idx_depth_1st, end) ==0 % clamp
                                        if collar_info(i, end)~=0 % not clamp
                                            table(i_table, 2) = depth_2nd;
                                            table(i_table, 3) = 3;
                                            collar_info(idx_depth_2nd, :) = collar_info(idx_depth_2nd, :) - 1;
                                            i_table = i_table + 1;
                                            flag_repeat = 1;
                                        else 
                                            table(i_table, 2) = depth_2nd;
                                            table(i_table, 3) = 3;
                                            collar_info(idx_depth_2nd, :) = collar_info(idx_depth_2nd, :) - 1;
                                            i_table = i_table + 1;
                                        end
                                    end       
                                end    
                            elseif num_inrange <= 0
                                if collar_info(i, end) ~=0 % not clamp
                                    table(i_table, 2) = table(i_table - 1, 2) + length;
                                    table(i_table, 3) = 3;
                                    i_table = i_table + 1;
                                    flag_repeat = 1;
                                else
                                    table(i_table, 2) = table(i_table - 1, 2) + length;
                                    table(i_table, 3) = 3;
                                    i_table = i_table + 1;
                                end
                            end   
                        end
                    end
                end
            elseif collar_info(i, p_No) ~= p_No    
                if collar_info(i, p_No) ~= 0
                    if abs(depth(i) - table(i_table - 1, 2)) <= (length+1*margin) && ...
                            abs(depth(i) - table(i_table - 1, 2)) >= (length-1*margin)
                        if collar_info(i, end)~=0 % not clamp
                            table(i_table, 2) = depth(i);
                            table(i_table, 3) = 2;
                            collar_info(i, :) = collar_info(i, :) - 1;
                            i_table = i_table + 1;
                        else
                            i = i + 1;
                            continue % is clamp
                        end       
                    else
                        if abs(depth(i) - table(i_table - 1, 2)) < (length-1*margin)
                            i = i + 1;
                            continue
                        elseif abs(depth(i) - table(i_table - 1, 2)) > (length+1*margin)
                            if abs(depth(i) - table(i_table - 1, 2)) <= (1.2 * length + margin)
                                i = i + 1;
                                continue
                            else
                                drct = 1;
                                [num_inrange, depth_1st, depth_2nd] = ...
                                    count_num(depth, table, i, p_No, collar_info, length, margin,...
                                    table(i_table - 1, 2), drct);
                                if num_inrange > 0
                                    % find index of depth_1st and depth_2nd
                                    idx_depth_1st = find(round(depth) == round(depth_1st));
                                    if depth_2nd ~= -99999
                                        idx_depth_2nd = find(round(depth) == round(depth_2nd));
                                    else
                                        idx_depth_2nd = NaN;
                                    end
                                    if num_inrange == 1
                                        if collar_info(idx_depth_1st, end)~=0 % not clamp
                                            table(i_table, 2) = depth_1st;
                                            table(i_table, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table = i_table + 1;
                                        elseif collar_info(idx_depth_1st, end) ==0 % clamp
                                            i = i + 1;
                                            continue
                                        end
                                    elseif num_inrange > 1
                                        if collar_info(idx_depth_1st, end)~=0 % not clamp
                                            table(i_table, 2) = depth_1st;
                                            table(i_table, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table = i_table + 1;
                                        elseif collar_info(idx_depth_1st, end) ==0 % clamp
                                            table(i_table, 2) = depth_2nd;
                                            table(i_table, 3) = 3;
                                            collar_info(idx_depth_2nd, :) = collar_info(idx_depth_2nd, :) - 1;
                                            i_table = i_table + 1;
                                        end
                                    end
                                elseif num_inrange <= 0
                                    table(i_table, 2) = table(i_table - 1, 2) + length;
                                    table(i_table, 3) = 3;
                                    i_table = i_table + 1;
                                end
                            end
                        end
                    end
                else 
                    i = i + 1;
                    continue
                end
            end
            if flag_repeat == 0
                i = i + 1;
            end
        end
    end
%%% count from anchor point to the top
    if p_No < 3 && p_No < size(collar_info, 2)
        i = anchor_idx;    
        table_len = 1;
        while i > 1 && table_len <= size(collar_info, 1)
            flag_repeat = 0;
            if collar_info(i, p_No) == p_No
                if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (length+1*margin) && ...
                        abs(depth(i) - table_r(i_table_r - 1, 2)) >= (length-1*margin)
                    if collar_info(i, end)~=0 % not clamp
                        table_r(i_table_r, 2) = depth(i);
                        table_r(i_table_r, 3) = 1;
                        collar_info(i, :) = collar_info(i, :) - 1;                        
                        % if next is high probablity to be green, the previous
                        % that match length rule is green as well
                        if table_r(i_table_r-1, 3) == 3 
                            table_r(i_table_r-1, 3) = 2;
                        elseif table_r(i_table_r-1, 3) == 2
                            table_r(i_table_r-1, 3) = 1;
                        end
                        i_table_r = i_table_r + 1;
                    else
                        i = i - 1;
                        continue % is clamp
                    end
                    
                else
                    if abs(depth(i) - table_r(i_table_r - 1, 2)) < (length-1*margin)
                        if abs(depth(i) - table_r(i_table_r - 1, 2)) <= 2* margin
                            if collar_info(i, end)~=0 % not clamp and replace last collar
                                i_restore = find(round(depth) == round(table_r(i_table_r - 1, 2)));
                                collar_info(i_restore, :) = collar_info(i_restore, :) + 1;
                                table_r(i_table_r - 1, 2) = depth(i);
                                table_r(i_table_r - 1, 3) = 2;
                                collar_info(i, :) = collar_info(i, :) - 1;
                            else
                                i = i - 1;
                                continue % is clamp
                            end
                        elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > 2* margin    
                            if table_r(i_table_r - 1, 3) ~= 1
                                if i_table_r > 2 && calc_mode(depth(i),table_r(i_table_r - 2, 2), length)...
                                        <= calc_mode(table_r(i_table_r - 1, 2),table_r(i_table_r - 2, 2), length)%
                                    if collar_info(i, end)~=0 % not clamp and replace last collar
                                        i_restore = find(round(depth) == round(table_r(i_table_r - 1, 2)));
                                        collar_info(i_restore, :) = collar_info(i_restore, :) + 1;
                                        table_r(i_table_r - 1, 2) = depth(i);
                                        table_r(i_table_r - 1, 3) = 2;
                                        collar_info(i, :) = collar_info(i, :) - 1;
                                    else
                                        i = i - 1;
                                        continue % is clamp
                                    end
                                elseif i_table_r <= 2 
                                    if collar_info(i, end)~=0 % not clamp and replace last collar
                                        table_r(i_table_r, 2) = depth(i);
                                        table_r(i_table_r, 3) = 2;
                                        collar_info(i, :) = collar_info(i, :) - 1;
                                        i_table_r = i_table_r + 1;
                                    else
                                        i = i - 1;
                                        continue % is clamp
                                    end
                                else
                                    i = i - 1;
                                    continue
                                end
                            elseif table_r(i_table_r - 1, 3) == 1
                                if collar_info(i, end)~=0 % not clamp and replace last collar
                                    table_r(i_table_r, 2) = depth(i);
                                    table_r(i_table_r, 3) = 2;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                    i_table_r = i_table_r + 1;
                                else
                                    i = i - 1;
                                    continue % is clamp
                                end                                
                            end
                        end
                    elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > (length+1*margin)
                        if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (1.2 * length + margin)
                            if collar_info(i, end)~=0 % not clamp
                                table_r(i_table_r, 2) = depth(i);
                                table_r(i_table_r, 3) = 2;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table_r = i_table_r + 1;
                            else
                                i = i - 1;
                                continue % is clamp
                            end 
                        elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > (1.2 * length + margin)
                            drct = -1;
                            [num_inrange, depth_1st, depth_2nd] = ...
                                count_num(depth, table_r, i, p_No, collar_info, length, margin,...
                                table_r(i_table_r - 1, 2), drct);
                            if num_inrange > 0
                                % find index of depth_1st and depth_2nd
                                idx_depth_1st = find(round(depth) == round(depth_1st));
                                if depth_2nd ~= -99999
                                    idx_depth_2nd = find(round(depth) == round(depth_2nd));
                                else
                                    idx_depth_2nd = NaN;
                                end
                                if num_inrange == 1
                                    if collar_info(idx_depth_1st, end)~=0 % not clamp
                                        if collar_info(i, end)~=0 % not clamp
                                            table_r(i_table_r, 2) = depth_1st;
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table_r = i_table_r + 1;
                                            flag_repeat = 1;
                                        else 
                                            table_r(i_table_r, 2) = depth_1st;
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table_r = i_table_r + 1;
                                        end
                                    elseif collar_info(idx_depth_1st, end) ==0 % clamp
                                        if collar_info(i, end)~=0 % not clamp
                                            table_r(i_table_r, 2) = depth(i);
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(i, :) = collar_info(i, :) - 1;
                                            i_table_r = i_table_r + 1;
                                        else
                                            i = i - 1;
                                            continue % is clamp
                                        end
                                    end             
                                elseif num_inrange > 1
                                    if collar_info(idx_depth_1st, end)~=0 % not clamp
                                        if collar_info(i, end)~=0 % not clamp
                                            table_r(i_table_r, 2) = depth_1st;
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table_r = i_table_r + 1;
                                            flag_repeat = 1;
                                        else 
                                            table_r(i_table_r, 2) = depth_1st;
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table_r = i_table_r + 1;
                                        end
                                    elseif collar_info(idx_depth_1st, end) ==0 % clamp
                                        if collar_info(i, end)~=0 % not clamp
                                            table_r(i_table_r, 2) = depth_2nd;
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(idx_depth_2nd, :) = collar_info(idx_depth_2nd, :) - 1;
                                            i_table_r = i_table_r + 1;
                                            flag_repeat = 1;
                                        else 
                                            table_r(i_table_r, 2) = depth_2nd;
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(idx_depth_2nd, :) = collar_info(idx_depth_2nd, :) - 1;
                                            i_table_r = i_table_r + 1;
                                        end
                                    end       
                                end    
                            elseif num_inrange <= 0
                                if collar_info(i, end) ~=0 % not clamp
                                    table_r(i_table_r, 2) = table_r(i_table_r - 1, 2) + length;
                                    table_r(i_table_r, 3) = 3;
                                    i_table_r = i_table_r + 1;
                                    flag_repeat = 1;
                                else
                                    table_r(i_table_r, 2) = table_r(i_table_r - 1, 2) + length;
                                    table_r(i_table_r, 3) = 3;
                                    i_table_r = i_table_r + 1;
                                end
                            end   
                        end
                    end
                end
            elseif collar_info(i, p_No) ~= p_No    
                if collar_info(i, p_No) ~= 0
                    if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (length+1*margin) && ...
                            abs(depth(i) - table_r(i_table_r - 1, 2)) >= (length-1*margin)
                        if collar_info(i, end)~=0 % not clamp
                            table_r(i_table_r, 2) = depth(i);
                            table_r(i_table_r, 3) = 2;
                            collar_info(i, :) = collar_info(i, :) - 1;
                            i_table_r = i_table_r + 1;
                        else
                            i = i - 1;
                            continue % is clamp
                        end       
                    else
                        if abs(depth(i) - table_r(i_table_r - 1, 2)) < (length-1*margin)
                            i = i - 1;
                            continue
                        elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > (length+1*margin)
                            if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (1.2 * length + margin)
                                i = i - 1;
                                continue
                            else
                                drct = -1;
                                [num_inrange, depth_1st, depth_2nd] = ...
                                    count_num(depth, table_r, i, p_No, collar_info, length, margin,...
                                    table_r(i_table_r - 1, 2), drct);
                                if num_inrange > 0
                                    % find index of depth_1st and depth_2nd
                                    idx_depth_1st = find(round(depth) == round(depth_1st));
                                    if depth_2nd ~= -99999
                                        idx_depth_2nd = find(round(depth) == round(depth_2nd));
                                    else
                                        idx_depth_2nd = NaN;
                                    end
                                    if num_inrange == 1
                                        if collar_info(idx_depth_1st, end)~=0 % not clamp
                                            table_r(i_table_r, 2) = depth_1st;
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table_r = i_table_r + 1;
                                        elseif collar_info(idx_depth_1st, end) ==0 % clamp
                                            i = i - 1;
                                            continue
                                        end
                                    elseif num_inrange > 1
                                        if collar_info(idx_depth_1st, end)~=0 % not clamp
                                            table_r(i_table_r, 2) = depth_1st;
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                            i_table_r = i_table_r + 1;
                                        elseif collar_info(idx_depth_1st, end) ==0 % clamp
                                            table_r(i_table_r, 2) = depth_2nd;
                                            table_r(i_table_r, 3) = 3;
                                            collar_info(idx_depth_2nd, :) = collar_info(idx_depth_2nd, :) - 1;
                                            i_table_r = i_table_r + 1;
                                        end
                                    end
                                elseif num_inrange <= 0
                                    table_r(i_table_r, 2) = table_r(i_table_r - 1, 2) + length;
                                    table_r(i_table_r, 3) = 3;
                                    i_table_r = i_table_r + 1;
                                end
                            end
                        end
                    end
                else 
                    i = i - 1;
                    continue
                end
            end
            if flag_repeat == 0
                i = i - 1;
            end
            table_len = size(find(table_r(:, 2)>0), 1);
        end
    end
    
    % p_No = 3, 4, 5 or outermost pipe, do not deal with root issue
    if p_No > 2 || p_No == size(collar_info, 2)
        % make decision from anchor to bottom
        i = anchor_idx;        
        while i <= size(collar_info, 1) && i >= 1
            flag_repeat = 0;
            if collar_info(i, p_No) == p_No
                if abs(depth(i) - table(i_table - 1, 2)) <= (length+1*margin) && ...
                        abs(depth(i) - table(i_table - 1, 2)) >= (length-1*margin)
                    table(i_table, 2) = depth(i);
                    table(i_table, 3) = 1;
                    collar_info(i, :) = collar_info(i, :) - 1;
                    % if next is high probablity to be green, the previous
                    % that match length rule is green as well
                    if table(i_table-1, 3) == 3
                        table(i_table-1, 3) = 2;
                    elseif table(i_table-1, 3) == 2
                        table(i_table-1, 3) = 1;
                    end
                    i_table = i_table + 1;
                    
                else
                    if abs(depth(i) - table(i_table - 1, 2)) < (length-1*margin)
                        if abs(depth(i) - table(i_table - 1, 2)) <= 2* margin
                            if i == anchor_idx
                                i = i + 1;
                                continue
                            else
                                i_restore = find(round(depth) == round(table(i_table - 1, 2)));
                                collar_info(i_restore, :) = collar_info(i_restore, :) + 1;
                                table(i_table - 1, 2) = depth(i);
                                table(i_table - 1, 3) = 2;
                                collar_info(i, :) = collar_info(i, :) - 1;
                            end
                        elseif abs(depth(i) - table(i_table - 1, 2)) > 2* margin
                            if table(i_table - 1, 3) ~= 1
                                if i_table > 2 && calc_mode(depth(i),table(i_table - 2, 2), length)...
                                        <= calc_mode(table(i_table - 1, 2),table(i_table - 2, 2), length)%
                                    i_restore = find(round(depth) == round(table(i_table - 1, 2)));
                                    collar_info(i_restore, :) = collar_info(i_restore, :) + 1;
                                    table(i_table - 1, 2) = depth(i);
                                    table(i_table - 1, 3) = 2;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                elseif i_table <= 2 
                                    table(i_table, 2) = depth(i);
                                    table(i_table, 3) = 2;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                    i_table = i_table + 1;
                                else
                                    i = i + 1;
                                    continue
                                end                        
                            elseif table(i_table - 1, 3) == 1
                                table(i_table, 2) = depth(i);
                                table(i_table, 3) = 2;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table = i_table + 1;
                            end
                        end
                    elseif abs(depth(i) - table(i_table - 1, 2)) > (length+1*margin)
                        if abs(depth(i) - table(i_table - 1, 2)) <= (1.2 * length + margin)
                            table(i_table, 2) = depth(i);
                            table(i_table, 3) = 2;
                            collar_info(i, :) = collar_info(i, :) - 1;
                            i_table = i_table + 1;
                        elseif abs(depth(i) - table(i_table - 1, 2)) > (1.2 * length + margin)
                            drct = 1;
                            [num_inrange, depth_1st, depth_2nd] = ...
                                count_num(depth, table, i, p_No, collar_info, length, margin,...
                                table(i_table - 1, 2), drct);
                            if num_inrange > 0
                                % find index of depth_1st and depth_2nd
                                idx_depth_1st = find(round(depth) == round(depth_1st));
                                if depth_2nd ~= -99999
                                    idx_depth_2nd = find(round(depth) == round(depth_2nd));
                                else
                                    idx_depth_2nd = NaN;
                                end
                                table(i_table, 2) = depth_1st;
                                table(i_table, 3) = 3;
                                collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                i_table = i_table + 1;
                                flag_repeat = 1;
                            elseif num_inrange <= 0
                                table(i_table, 2) = table(i_table - 1, 2) + length;
                                table(i_table, 3) = 3;
                                i_table = i_table + 1;
                                flag_repeat = 1;
                            end
                        end
                    end
                end
            elseif collar_info(i, p_No) ~= p_No
                if collar_info(i, p_No) ~= 0
                    if abs(depth(i) - table(i_table - 1, 2)) <= (length+1*margin) && ...
                            abs(depth(i) - table(i_table - 1, 2)) >= (length-1*margin)
                        table(i_table, 2) = depth(i);
                        table(i_table, 3) = 2;
                        collar_info(i, :) = collar_info(i, :) - 1;
                        i_table = i_table + 1;
                    else
                        if abs(depth(i) - table(i_table - 1, 2)) < (length-1*margin)
                            i = i + 1;
                            continue
                        elseif abs(depth(i) - table(i_table - 1, 2)) > (length+1*margin)
                            if abs(depth(i) - table(i_table - 1, 2)) <= (1.2 * length + margin)
                                i = i + 1;
                                continue
                            else
                                drct = 1;
                                [num_inrange, depth_1st, depth_2nd] = ...
                                    count_num(depth, table, i, p_No, collar_info, length, margin,...
                                    table(i_table - 1, 2), drct);
                                if num_inrange > 0
                                    % find index of depth_1st and depth_2nd
                                    idx_depth_1st = find(round(depth) == round(depth_1st));
                                    if depth_2nd ~= -99999
                                        idx_depth_2nd = find(round(depth) == round(depth_2nd));
                                    else
                                        idx_depth_2nd = NaN;
                                    end
                                    table(i_table, 2) = depth_1st;
                                    table(i_table, 3) = 3;
                                    collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                    i_table = i_table + 1;                                    
                                elseif num_inrange <= 0
                                    table(i_table, 2) = table(i_table - 1, 2) + length;
                                    table(i_table, 3) = 3;
                                    i_table = i_table + 1;
                                end
                            end
                        end
                    end
                else
                    i = i + 1;
                    continue
                end
            end
            if flag_repeat == 0
                i = i + 1;
            end
        end
    end
    %%% count from anchor point to the top
    if p_No > 2 || p_No == size(collar_info, 2)
        i = anchor_idx; 
        table_len = 1;
        while i > 1 && table_len <= size(collar_info, 1)
            flag_repeat = 0;
            if collar_info(i, p_No) == p_No
                if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (length+1*margin) && ...
                        abs(depth(i) - table_r(i_table_r - 1, 2)) >= (length-1*margin)
                    table_r(i_table_r, 2) = depth(i);
                    table_r(i_table_r, 3) = 1;
                    collar_info(i, :) = collar_info(i, :) - 1;
                    % if next is high probablity to be green, the previous
                    % that match length rule is green as well
                    if table_r(i_table_r-1, 3) == 3
                        table_r(i_table_r-1, 3) = 2;
                    elseif table_r(i_table_r-1, 3) == 2
                        table_r(i_table_r-1, 3) = 1;
                    end
                    i_table_r = i_table_r + 1;                    
                else
                    if abs(depth(i) - table_r(i_table_r - 1, 2)) < (length-1*margin)
                        if abs(depth(i) - table_r(i_table_r - 1, 2)) <= 2* margin
                            i_restore = find(round(depth) == round(table_r(i_table_r - 1, 2)));
                            collar_info(i_restore, :) = collar_info(i_restore, :) + 1;
                            table_r(i_table_r - 1, 2) = depth(i);
                            table_r(i_table_r - 1, 3) = 2;
                            collar_info(i, :) = collar_info(i, :) - 1;
                        elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > 2* margin
                            if table_r(i_table_r - 1, 3) ~= 1
                                if i_table_r > 2 && calc_mode(depth(i),table_r(i_table_r - 2, 2), length)...
                                        <= calc_mode(table_r(i_table_r - 1, 2),table_r(i_table_r - 2, 2), length)%
                                    i_restore = find(round(depth) == round(table_r(i_table_r - 1, 2)));
                                    collar_info(i_restore, :) = collar_info(i_restore, :) + 1;
                                    table_r(i_table_r - 1, 2) = depth(i);
                                    table_r(i_table_r - 1, 3) = 2;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                elseif i_table_r <= 2 
                                    table_r(i_table_r, 2) = depth(i);
                                    table_r(i_table_r, 3) = 2;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                    i_table_r = i_table_r + 1;
                                else
                                    i = i - 1;
                                    continue
                                end   
                            elseif table_r(i_table_r - 1, 3) == 1
                                table_r(i_table_r, 2) = depth(i);
                                table_r(i_table_r, 3) = 2;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table_r = i_table_r + 1;
                            end
                        end
                    elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > (length+1*margin)
                        if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (1.2 * length + margin)
                            table_r(i_table_r, 2) = depth(i);
                            table_r(i_table_r, 3) = 2;
                            collar_info(i, :) = collar_info(i, :) - 1;
                            i_table_r = i_table_r + 1;
                        elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > (1.2 * length + margin)
                            drct = -1;
                            [num_inrange, depth_1st, depth_2nd] = ...
                                count_num(depth, table_r, i, p_No, collar_info, length, margin,...
                                table_r(i_table_r - 1, 2), drct);
                            if num_inrange > 0
                                % find index of depth_1st and depth_2nd
                                idx_depth_1st = find(round(depth) == round(depth_1st));
                                if depth_2nd ~= -99999
                                    idx_depth_2nd = find(round(depth) == round(depth_2nd));
                                else
                                    idx_depth_2nd = NaN;
                                end
                                table_r(i_table_r, 2) = depth_1st;
                                table_r(i_table_r, 3) = 3;
                                collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                i_table_r = i_table_r + 1; 
                                i = i -1;
                            elseif num_inrange <= 0
                                table_r(i_table_r, 2) = table_r(i_table_r - 1, 2) - length;
                                table_r(i_table_r, 3) = 3;
                                i_table_r = i_table_r + 1;
                                i = i - 1;
                            end
                        end
                    end
                end
            elseif collar_info(i, p_No) ~= p_No
                if collar_info(i, p_No) ~= 0
                    if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (length+1*margin) && ...
                            abs(depth(i) - table_r(i_table_r - 1, 2)) >= (length-1*margin)
                        table_r(i_table_r, 2) = depth(i);
                        table_r(i_table_r, 3) = 2;
                        collar_info(i, :) = collar_info(i, :) - 1;
                        i_table_r = i_table_r + 1;
                    else
                        if abs(depth(i) - table_r(i_table_r - 1, 2)) < (length-1*margin)
                            i = i - 1;
                            continue
                        elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > (length+1*margin)
                            if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (1.2 * length + margin)
                                i = i - 1;
                                continue
                            else
                                drct = -1;
                                [num_inrange, depth_1st, depth_2nd] = ...
                                    count_num(depth, table_r, i, p_No, collar_info, length, margin,...
                                    table_r(i_table_r - 1, 2), drct);
                                if num_inrange > 0
                                    % find index of depth_1st and depth_2nd
                                    idx_depth_1st = find(round(depth) == round(depth_1st));
                                    if depth_2nd ~= -99999
                                        idx_depth_2nd = find(round(depth) == round(depth_2nd));
                                    else
                                        idx_depth_2nd = NaN;
                                    end
                                    table_r(i_table_r, 2) = depth_1st;
                                    table_r(i_table_r, 3) = 3;
                                    collar_info(idx_depth_1st, :) = collar_info(idx_depth_1st, :) - 1;
                                    i_table_r = i_table_r + 1;
                                elseif num_inrange <= 0
                                    table_r(i_table_r, 2) = table_r(i_table_r - 1, 2) - length;
                                    table_r(i_table_r, 3) = 3;
                                    i_table_r = i_table_r + 1;
                                end
                            end
                        end
                    end
                else
                    i = i - 1;
                    continue
                end
            end
            if flag_repeat == 0
                i = i - 1;
            end
            table_len = size(find(table_r(:, 2)>0), 1);
        end
    end

%%% compose two tables into one
        table_r_len = find(table_r(:, 2) == -999, 1, 'first') - 1;
        table_r = flip(table_r(1:table_r_len, :));
        table_len = find(table(:, 2) == -999, 1, 'first') - 1;
        table = table(1:table_len, :);
        table_final = [table_r(1:end, :); table];
        % sort according to column 2
        table_final = sortrows(table_final, 2);
        if p_No == 1
            table_P1 = table_final;
            % replaced the closest collar depth with anchor depth
            [~, anchr_idx] = min(abs(table_P1(:, 2) - anchor_depth(1)));
            achr_depth = table_P1(anchr_idx, 2);
            table_P1(table_P1(:, 2) == achr_depth, 3) = 1;
            table_P1(table_P1(:, 2) == achr_depth, 2) = anchor_depth(1);
        elseif p_No == 2
            table_P2 = table_final;
            [~, anchr_idx] = min(abs(table_P2(:, 2) - anchor_depth(2)));
            achr_depth = table_P2(anchr_idx, 2);
            table_P2(table_P2(:, 2) == achr_depth, 3) = 1;
            table_P2(table_P2(:, 2) == achr_depth, 2) = anchor_depth(2);
        elseif p_No == 3
            table_P3 = table_final;
            [~, anchr_idx] = min(abs(table_P3(:, 2) - anchor_depth(3)));
            achr_depth = table_P3(anchr_idx, 2);
            table_P3(table_P3(:, 2) == achr_depth, 3) = 1;
            table_P3(table_P3(:, 2) == achr_depth, 2) = anchor_depth(3);
        elseif p_No == 4
            table_P4 = table_final;
            [~, anchr_idx] = min(abs(table_P4(:, 2) - anchor_depth(4)));
            achr_depth = table_P4(anchr_idx, 2);
            table_P4(table_P4(:, 2) == achr_depth, 3) = 1;
            table_P4(table_P4(:, 2) == achr_depth, 2) = anchor_depth(4);
        elseif p_No == 5
            table_P5 = table_final;
            [~, anchr_idx] = min(abs(table_P5(:, 2) - anchor_depth(5)));
            achr_depth = table_P5(anchr_idx, 2);
            table_P5(table_P5(:, 2) == achr_depth, 3) = 1;
            table_P5(table_P5(:, 2) == achr_depth, 2) = anchor_depth(5);
        end        
        

end
                    
end    
                        
                        


