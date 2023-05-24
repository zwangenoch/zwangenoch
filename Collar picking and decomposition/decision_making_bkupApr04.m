clear;
close all; clc;

% input collar info from find_collar.m
collar_info = [0	2	3
1	2	3
1	2	3
1	2	3
0	0	0
0	0	3
0	2	3
0	0	0
1	2	3
0	0	3
0	2	3
0	0	0
1	2	3
0	0	3
0	0	0
0	2	3
0	0	0
1	2	3
0	0	3
0	2	3
0	0	0
1	2	3
0	0	3
0	2	3
0	0	0
0	0	0
1	2	3
0	0	3
0	2	3
0	0	0
0	0	0
0	0	0
];
depth = [180.03288
194.782821
195.699484
196.866146
198.616139
205.782777
217.699396
229.616015
234.865994
244.449289
255.615911
271.865846
273.532506
283.449133
287.199118
296.199082
309.282363
314.782341
321.198982
335.03226
346.948879
353.865518
360.698824
370.03212
379.448749
385.198726
393.198694
399.782001
409.781961
422.781909
432.115205
433.115201
];

%% initialization
% input specified pipe length-length between joints, and margin, in feet
len_spec = [40, 38, 39, 39, 39.6]; % from outer to inner
margin = 2;
anchor_depth = [234.9, 217.7, 205.8, 385.2, 198.6];

% browse according to pipe number
for p_No = 1:size(collar_info, 2)
    % assign length
    length = len_spec(p_No);
    

% plan to let user to enter anchor point depth for each pipe   
%     % find anchor point for 1-2 pipe and 3-5 pipe
%     anchor_idx = -999;
%     % create collar depth table
%     anchor_table = zeros(1200, 1);
%     i_anchor = 1;
    
%     for i = 1:size(collar_info, 1)
%         if p_No==1 && collar_info(i, p_No) == p_No
%             anchor_table(i_anchor) = depth(i, p_No);
%             i_anchor = i_anchor + 1;
%         elseif p_No<3 && p_No>1 && collar_info(i, p_No) == p_No && ...
%                 collar_info(i, p_No-1) == 0
%             anchor_table(i_anchor) = depth(i, p_No);
%             i_anchor = i_anchor + 1;
%         elseif p_No>2 && p_No>1 && collar_info(i, p_No) == p_No && ...
%                 collar_info(i, p_No-1) == 0
%             anchor_table(i_anchor) = depth(i, p_No);
%             i_anchor = i_anchor + 1;
%         end
%     end
%     for index = 1:size(collar_info, 1)
%         % anchor point for 1st and 2nd pipe
%         if p_No < 3 && index>1 && index<(size(collar_info, 1)-1) &&...
%                 abs(anchor_table(index-1)-anchor_table(index))<=(length+margin) && ...
%                 abs(anchor_table(index-1)-anchor_table(index))>=(length-margin) && ...
%                 abs(anchor_table(index+1)-anchor_table(index))<=(length+margin) && ...
%                 abs(anchor_table(index+1)-anchor_table(index))>=(length-margin) && ...
%                 ((abs(anchor_table(index+2)-anchor_table(index+1))<=((length+margin)) && ...
%                 abs(anchor_table(index+2)-anchor_table(index+1))>=((length-margin))) || ...
%                 (abs(anchor_table(index-2)-anchor_table(index-1))<=((length+margin)) && ...
%                 abs(anchor_table(index-2)-anchor_table(index-1))>=((length-margin))))
%             anchor_depth = anchor_table(index);
%             anchor_idx = find(depth == anchor_depth, 1);
%         elseif p_No > 2 && index>1 && index<(size(collar_info, 1)-1) &&...
%                 abs(anchor_table(index-1)-anchor_table(index))<=(length+margin) && ...
%                 abs(anchor_table(index-1)-anchor_table(index))>=(length-margin) && ...
%                 abs(anchor_table(index+1)-anchor_table(index))<=(length+margin) && ...
%                 abs(anchor_table(index+1)-anchor_table(index))>=(length-margin) && ...
%                 ((abs(anchor_table(index+2)-anchor_table(index+1))<=((length+margin)) && ...
%                 abs(anchor_table(index+2)-anchor_table(index+1))>=((length-margin))) || ...
%                 (abs(anchor_table(index-2)-anchor_table(index-1))<=((length+margin)) && ...
%                 abs(anchor_table(index-2)-anchor_table(index-1))>=((length-margin))))
%             anchor_depth = anchor_table(index);
%             anchor_idx = find(depth == anchor_depth, 1);
%         elseif anchor_table(index) == 0
%             break
%         end
%     end
%     for index = 1:size(collar_info, 1)
%         if anchor_idx == -999
%             if p_No < 3 && index>1 && index<(size(collar_info, 1)-1) &&...
%                     abs(anchor_table(index-1)-anchor_table(index))<=(length+margin) && ...
%                     abs(anchor_table(index-1)-anchor_table(index))>=(length-margin) && ...
%                     abs(anchor_table(index+1)-anchor_table(index))<=(length+margin) && ...
%                     abs(anchor_table(index+1)-anchor_table(index))>=(length-margin)
%                 anchor_depth = anchor_table(index);
%                 anchor_idx = find(depth == anchor_depth, 1);
%             elseif p_No > 2 && index>1 && index<(size(collar_info, 1)-1) &&...
%                     abs(anchor_table(index-1)-anchor_table(index))<=(length+margin) && ...
%                     abs(anchor_table(index-1)-anchor_table(index))>=(length-margin) && ...
%                     abs(anchor_table(index+1)-anchor_table(index))<=(length+margin) && ...
%                     abs(anchor_table(index+1)-anchor_table(index))>=(length-margin)
%                 anchor_depth = anchor_table(index);
%                 anchor_idx = find(depth == anchor_depth, 1);
%             elseif anchor_table(index) == 0
%                 break
%             end
%         end
%     end
%     if anchor_idx == -999
%         msgbox('cannot find anchor');
%     end
    
    anchor_idx = find(depth == anchor_depth(p_No), 1, 'first');
    % generate collar table filled with 0s
    table = zeros(1200, 3);
    table(:, 1) = (1:1200)';
    table(1, 2) = depth(anchor_idx);
    table(1, 3) = 1;
    i_table = 2;
    % generate reverse collar table filled with 0s from anchor to top
    table_r = zeros(1200, 3);
    table_r(:, 1) = (1:1200)';
    table_r(1, 2) = depth(anchor_idx);
    table_r(1, 3) = 1;
    i_table_r = 2;
        
    % Pipe 1-2: browse all possible collars along depth   
    if p_No < 3 
        % make decision from anchor to bottom
        for i = anchor_idx:size(collar_info, 1)
            if collar_info(i, p_No) == p_No
                if abs(depth(i) - table(i_table - 1, 2)) <= (length+1*margin) && ...
                        abs(depth(i) - table(i_table - 1, 2)) >= (length-1*margin)
                    if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                            collar_info(i, end)~=0 % not clamp
                        table(i_table, 2) = depth(i);
                        table(i_table, 3) = 1;
                        collar_info(i, :) = collar_info(i, :) - 1;                        
                        % if next is high probablity to be green, the previous
                        % that match length rule is green as well
                        if table(i_table-1, 3) > 0 && table(i_table-1, 3) < 0.5
                            table(i_table-1, 3) = 0.5;
                        elseif table(i_table-1, 3) > 0 && table(i_table-1, 3) < 1
                            table(i_table-1, 3) = 1;
                        end
                        i_table = i_table + 1;
                    elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                        table(i_table, 2) = depth(i);
                        table(i_table, 3) = 1;
                        collar_info(i, :) = collar_info(i, :) - 1;
                        if table(i_table-1, 3) > 0 && table(i_table-1, 3) < 0.5
                            table(i_table-1, 3) = 0.5;
                        elseif table(i_table-1, 3) > 0 && table(i_table-1, 3) < 1
                            table(i_table-1, 3) = 1;
                        end
                        i_table = i_table + 1;
                    else
                        continue % is clamp
                    end
                    
                else
                    if abs(depth(i) - table(i_table - 1, 2)) < (length-1*margin)
                        if abs(depth(i) - table(i_table - 1, 2)) < 2* margin
                            if table(i_table - 1, 3) ~= 1
                                if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                        collar_info(i, end)~=0 % not clamp and replace last collar
                                    table(i_table - 1, 2) = depth(i);
                                    table(i_table - 1, 3) = 0.5;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                    table(i_table - 1, 2) = depth(i);
                                    table(i_table - 1, 3) = 0.5;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                else
                                    continue % is clamp
                                end
                            elseif table(i_table - 1, 3) == 1
                                rmd_current = -999;
                                rmd_previous = -998;
                                % i_second_last = find(table(:,3)==1, 2, 'last');
                                % depth_second_last = table(i_second_last, 2);
                                rmd_current = abs(mod(abs(depth(i) - table(i_table-2, 2)), length));
                                rmd_previous = abs(mod(abs(table(i_table-1, 2)-table(i_table-2, 2)), length));
                                if rmd_current <= rmd_previous
                                    if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                            collar_info(i, end)~=0 % not clamp and replace last collar
                                        table(i_table - 1, 2) = depth(i);
                                        table(i_table - 1, 3) = 0.5;
                                        collar_info(i, :) = collar_info(i, :) - 1;
                                    elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                        table(i_table - 1, 2) = depth(i);
                                        table(i_table - 1, 3) = 0.5;
                                        collar_info(i, :) = collar_info(i, :) - 1;
                                    else
                                        continue % is clamp
                                    end
                                elseif rmd_current > rmd_previous
                                    if abs(table(i_table-1, 2)-table(i_table-2, 2)) < 20
                                        if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                                collar_info(i, end)~=0 % not clamp and replace last collar
                                            table(i_table - 1, 2) = depth(i);
                                            table(i_table - 1, 3) = 0.5;
                                            collar_info(i, :) = collar_info(i, :) - 1;
                                        elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                            table(i_table - 1, 2) = depth(i);
                                            table(i_table - 1, 3) = 0.5;
                                            collar_info(i, :) = collar_info(i, :) - 1;
                                        else
                                            continue % is clamp
                                        end
                                    else
                                        collar_info(i, :) = collar_info(i, :) - 1;
                                        continue
                                    end
                                end
                            end
                        else
                            if table(i_table - 1, 3) ~= 1
                                rmd_current = -999;
                                rmd_previous = -998;
                                % i_second_last = find(table(:,3)==1, 2, 'last');
                                % depth_second_last = table(i_second_last, 2);
                                rmd_current = abs(mod(abs(depth(i) - table(i_table-2, 2)), length));
                                rmd_previous = abs(mod(abs(table(i_table-1, 2)-table(i_table-2, 2)), length));
                                if rmd_current <= rmd_previous
                                    if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                            collar_info(i, end)~=0 % not clamp 
                                        table(i_table, 2) = depth(i);
                                        table(i_table, 3) = 0.5;
                                        collar_info(i, :) = collar_info(i, :) - 1;
                                        i_table = i_table + 1;
                                    elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                        table(i_table, 2) = depth(i);
                                        table(i_table, 3) = 0.5;
                                        collar_info(i, :) = collar_info(i, :) - 1;
                                        i_table = i_table + 1;
                                    else
                                        continue % is clamp
                                    end
                                elseif rmd_current > rmd_previous
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                    continue
                                end
                            elseif table(i_table - 1, 3) == 1
                                if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                        collar_info(i, end)~=0 % not clamp
                                    table(i_table, 2) = depth(i);
                                    table(i_table, 3) = 0.5;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                    i_table = i_table + 1;
                                elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                    table(i_table, 2) = depth(i);
                                    table(i_table, 3) = 0.5;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                    i_table = i_table + 1;
                                else
                                    continue % is clamp
                                end
                            end
                        end
                    elseif abs(depth(i) - table(i_table - 1, 2)) > (length+1*margin)
                        if abs(depth(i) - table(i_table - 1, 2)) <= (length+1.2*margin)
                            if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                    collar_info(i, end)~=0 % not clamp
                                table(i_table, 2) = depth(i);
                                table(i_table, 3) = 0.5;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table = i_table + 1;
                            elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                table(i_table, 2) = depth(i);
                                table(i_table, 3) = 0.5;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table = i_table + 1;
                            else
                                continue % is clamp
                            end                        
                        elseif abs(depth(i) - table(i_table - 1, 2)) > (length+1.2*margin)
                            num_inrange = count_num(depth, table, i, p_No, collar_info, length, margin, table(i_table - 1, 2))
                                
                                
                                
                                
                                
                                
                    end
                    
                        
                        
                        
                        
                        
                    elseif abs(depth(i) - table(i_table - 1, 2)) < (length-1*margin) && ...
                            abs(depth(i) - table(i_table - 1, 2)) < 2* margin
                        
                        if table
                            if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                collar_info(i, end)~=0 % not clamp
                            table(i_table, 2) = depth(i);
                            table(i_table, 3) = 0.5;
                            collar_info(i, :) = collar_info(i, :) - 1;
                            i_table = i_table + 1;
                        elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                            table(i_table, 2) = depth(i);
                            table(i_table, 3) = 0.5;
                            collar_info(i, :) = collar_info(i, :) - 1;
                            i_table = i_table + 1;
                        else
                            collar_info(i, :) = collar_info(i, :) * 0;
                            continue % is clamp
                        end
                    elseif abs(depth(i) - table(i_table - 1, 2)) > (length+1*margin)
                        if table(i_table-1, 3) ~= 1
                            if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                    collar_info(i, end)~=0 % not clamp
                                table(i_table, 2) = depth(i);
                                table(i_table, 3) = 0.5;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table = i_table + 1;
                            elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                table(i_table, 2) = depth(i);
                                table(i_table, 3) = 0.5;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table = i_table + 1;
                            else
                                collar_info(i, :) = collar_info(i, :) * 0;
                                continue % is clamp
                            end
                        elseif table(i_table-1, 3) == 1
                            % find second last probability==1
                            rmd_current = -999;
                            rmd_previous = -999;
%                             i_second_last = find(table(:,3)==1, 2, 'last');
%                             depth_second_last = table(i_second_last, 2);
                            rmd_current = abs(mod(abs(table(i_table-1, 2)-depth(i)), length));
                            rmd_previous = abs(mod(abs(table(i_table-1, 2)-table(i_table-1, 3)), length));
                            if rmd_current < rmd_previous
                                if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                        collar_info(i, end)~=0 % not clamp
                                    table(i_table, 2) = depth(i);
                                    table(i_table, 3) = 0.5;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                    table(i_table, 2) = depth(i);
                                    table(i_table, 3) = 0.5;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                else
                                    collar_info(i, :) = collar_info(i, :) * 0;
                                    continue % is clamp
                                end
                            else
                                collar_info(i, :) = collar_info(i, :) - 1;
                                continue
                            end
                        end
                    end
                end
                
            elseif collar_info(i, p_No) ~= p_No && collar_info(i, p_No)~=0
                if abs(depth(i) - table(i_table - 1, 2)) <= (length+1*margin) && ...
                        abs(depth(i) - table(i_table - 1, 2)) >= (length-1*margin)
                    if collar_info(i, end-1)~=0 && collar_info(i, end)~=0 % not clamp
                        table(i_table, 2) = depth(i);
                        table(i_table, 3) = 0.5;
                        collar_info(i, :) = collar_info(i, :) - 1;
                        i_table = i_table + 1;
                    else
                        collar_info(i, :) = collar_info(i, :) * 0;
                        continue % is clamp
                    end
                else
                    if abs(depth(i) - table(i_table - 1, 2)) < (length-1*margin)
                        collar_info(i, :) = collar_info(i, :) - 1;
                        continue % is clamp
                    elseif abs(depth(i) - table(i_table - 1, 2)) > (length+1*margin)
                        i_info_temp = 1;
                        rmd_min = 999;
                        for m = i-1:-1:1
                            if p_No+1<=size(collar_info, 2)
                                if abs(depth(m) - table(i_table - 1, 2)) <= (length+1*margin)&&...
                                        abs(depth(m) - table(i_table - 1, 2)) >= (length-1*margin) &&...
                                        (collar_info(m, p_No)>0 || collar_info(m, p_No+1)>0)
                                    rmd_temp = abs(mod(abs(table(i_table-1, 2)-depth(i)), length));
                                    if rmd_temp <= rmd_min
                                        rmd_min = rmd_temp;
                                        idx_pick1 = m;
                                    end
                                elseif abs(depth(m) - table(i_table - 1, 2)) < (length-1*margin)
                                    break
                                end
                            elseif p_No+1>size(collar_info, 2)
                                if abs(depth(m) - table(i_table - 1, 2)) <= (length+1*margin)&&...
                                        abs(depth(m) - table(i_table - 1, 2)) >= (length-1*margin) &&...
                                        collar_info(m, p_No)>0
                                    rmd_temp = abs(mod(abs(table(i_table-1, 2)-depth(i)), length));
                                    if rmd_temp <= rmd_min
                                        rmd_min = rmd_temp;
                                        idx_pick1 = m;
                                    end
                                elseif abs(depth(m) - table(i_table - 1, 2)) < (length-1*margin)
                                    break
                                end
                            end
                        end
                        if rmd_min == 999
                            table(i_table, 2) = depth(i);
                            table(i_table, 3) = 0.25;
                            collar_info(i, :) = collar_info(i, :) - 1;
                            i_table = i_table + 1;
                        else
                            table(i_table, 2) = depth(idx_pick1);
                            table(i_table, 3) = 0.25;
                            collar_info(idx_pick1, :) = collar_info(idx_pick1, :) - 1;
                            i_table = i_table + 1;
                        end
                    end
                end
            elseif collar_info(i, p_No) ~= p_No && collar_info(i, p_No)==0
                collar_info(i, :) = collar_info(i, :) - 1;
                continue
            end
        end
        % make decision from anchor to top
        for i = anchor_idx:-1:1
            if collar_info(i, p_No) == p_No
                if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (length+1*margin) && ...
                        abs(depth(i) - table_r(i_table_r - 1, 2)) >= (length-1*margin)
                    table_r(i_table_r, 2) = depth(i);
                    table_r(i_table_r, 3) = 1;
                    collar_info(i, :) = collar_info(i, :) - 1;
                    i_table_r = i_table_r + 1;
                    % if next is high probablity to be green, the previous
                    % that match length rule is green as well
                    if table_r(i_table_r-1, 3) > 0 && table_r(i_table_r-1, 3) < 0.5
                        table_r(i_table_r-1, 3) = 0.5;
                    elseif table_r(i_table_r-1, 3) > 0 && table_r(i_table_r-1, 3) < 1
                        table_r(i_table_r-1, 3) = 1;
                    end
                else
                    if abs(depth(i) - table_r(i_table_r - 1, 2)) < (length-1*margin)
                        if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                collar_info(i, end)~=0 % not clamp
                            table_r(i_table_r, 2) = depth(i);
                            table_r(i_table_r, 3) = 0.5;
                            collar_info(i, :) = collar_info(i, :) - 1;
                            i_table_r = i_table_r + 1;
                        elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                            table_r(i_table_r, 2) = depth(i);
                            table_r(i_table_r, 3) = 0.5;
                            collar_info(i, :) = collar_info(i, :) - 1;
                            i_table_r = i_table_r + 1;
                        else
                            collar_info(i, :) = collar_info(i, :) * 0;
                            continue % is clamp
                        end
                    elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > (length+1*margin)
                        if table_r(i_table_r-1, 3) ~= 1
                            if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                    collar_info(i, end)~=0 % not clamp
                                table_r(i_table_r, 2) = depth(i);
                                table_r(i_table_r, 3) = 0.5;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table_r = i_table_r + 1;
                            elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                table_r(i_table_r, 2) = depth(i);
                                table_r(i_table_r, 3) = 0.5;
                                collar_info(i, :) = collar_info(i, :) - 1;
                                i_table_r = i_table_r + 1;
                            else
                                collar_info(i, :) = collar_info(i, :) * 0;
                                continue % is clamp
                            end
                        elseif table_r(i_table_r-1, 3) == 1
                            % find second last probability==1
                            rmd_current = -999;
                            rmd_previous = -999;
                            %                             i_second_last = find(table_r(:,3)==1, 2, 'last');
                            %                             depth_second_last = table_r(i_second_last, 2);
                            rmd_current = abs(mod(abs(table_r(i_table_r-1, 2)-depth(i)), length));
                            rmd_previous = abs(mod(abs(table_r(i_table_r-1, 2)-table_r(i_table_r-1, 3)), length));
                            if rmd_current < rmd_previous
                                if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                                        collar_info(i, end)~=0 % not clamp
                                    table_r(i_table_r, 2) = depth(i);
                                    table_r(i_table_r, 3) = 0.5;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                                    table_r(i_table_r, 2) = depth(i);
                                    table_r(i_table_r, 3) = 0.5;
                                    collar_info(i, :) = collar_info(i, :) - 1;
                                else
                                    collar_info(i, :) = collar_info(i, :) * 0;
                                    continue % is clamp
                                end
                            else
                                collar_info(i, :) = collar_info(i, :) - 1;
                                continue
                            end
                        end
                    end
                end
                
            elseif collar_info(i, p_No) ~= p_No && collar_info(i, p_No)~=0
                if abs(depth(i) - table_r(i_table_r - 1, 2)) <= (length+1*margin) && ...
                        abs(depth(i) - table_r(i_table_r - 1, 2)) >= (length-1*margin)
                    if size(collar_info, 2)>1 && collar_info(i, end-1)~=0 &&...
                            collar_info(i, end)~=0 % not clamp
                        table_r(i_table_r, 2) = depth(i);
                        table_r(i_table_r, 3) = 0.5;
                        collar_info(i, :) = collar_info(i, :) - 1;
                        i_table_r = i_table_r + 1;
                    elseif size(collar_info, 2)==1 && collar_info(i, end)~=0 % not clamp
                        table_r(i_table_r, 2) = depth(i);
                        table_r(i_table_r, 3) = 0.5;
                        collar_info(i, :) = collar_info(i, :) - 1;
                        i_table_r = i_table_r + 1;
                    else
                        collar_info(i, :) = collar_info(i, :) * 0;
                        continue % is clamp
                    end
                else
                    if abs(depth(i) - table_r(i_table_r - 1, 2)) < (length-1*margin)
                        collar_info(i, :) = collar_info(i, :) - 1;
                        continue 
                    elseif abs(depth(i) - table_r(i_table_r - 1, 2)) > (length+1*margin)
                        i_info_temp = 1;
                        rmd_min = 999;
                        for m = (i+1):1:size(collar_info, 1)
                            if p_No+1<=size(collar_info, 2)
                                if abs(depth(m) - table_r(i_table_r - 1, 2)) <= (length+1*margin)&&...
                                        abs(depth(m) - table_r(i_table_r - 1, 2)) >= (length-1*margin) &&...
                                        (collar_info(m, p_No)>0 || collar_info(m, p_No+1)>0)
                                    rmd_temp = abs(mod(abs(table_r(i_table_r-1, 2)-depth(i)), length));
                                    if rmd_temp <= rmd_min
                                        rmd_min = rmd_temp;
                                        idx_pick1 = m;
                                    end
                                elseif abs(depth(m) - table_r(i_table_r - 1, 2)) < (length-1*margin)
                                    break
                                end
                            elseif p_No+1>size(collar_info, 2)
                                if abs(depth(m) - table_r(i_table_r - 1, 2)) <= (length+1*margin)&&...
                                        abs(depth(m) - table_r(i_table_r - 1, 2)) >= (length-1*margin) &&...
                                        collar_info(m, p_No)>0
                                    rmd_temp = abs(mod(abs(table_r(i_table_r-1, 2)-depth(i)), length));
                                    if rmd_temp <= rmd_min
                                        rmd_min = rmd_temp;
                                        idx_pick1 = m;
                                    end
                                elseif abs(depth(m) - table_r(i_table_r - 1, 2)) < (length-1*margin)
                                    break
                                end
                            end
                        end
                    end
                end
            elseif collar_info(i, p_No) ~= p_No && collar_info(i, p_No)==0
                collar_info(i, :) = collar_info(i, :) - 1;
                continue
            end
        end
        % discard tables zeros part
        h1 = find(table_r(:, 2)==0, 1, 'first');
        table_top = flip(table_r(1:h-1, 2:3), 1);
        h2 = find(table(:, 2)==0, 1, 'first');
        table_bot = flip(table(1:h-1, 2:3), 1);
        if p_No == 1
            table_1P = [table_top; table_bot];
        elseif p_No == 2
            table_2P = [table_top; table_bot];
        end
            
            
            
            
    %elseif p_No >= 3 &&  p_No <= 5
            
    end
    
end
aaa=1;

