function [num_inrange, depth_1st, depth_2nd] = ...
    count_num(depth, table, i, p_No, collar_info, length, margin,...
    depth_collar, drct)
% function to count how many number in range on valid channels
% rank top 1-2 numbers

% find valid channels for different cases - star_CH
if size(collar_info, 2) == 1
    star_CH = [1];
elseif size(collar_info, 2) == 2
    if p_No == 1
        star_CH = [1, 2];
    elseif p_No == 2
        star_CH = [1, 2];
    end
elseif size(collar_info, 2) == 3
    if p_No == 1
        star_CH = [1, 2];
    elseif p_No == 2
        star_CH = [1, 2, 3];
    elseif p_No == 3
        star_CH = [2, 3];
    end
elseif size(collar_info, 2) == 4
    if p_No == 1
        star_CH = [1, 2];
    elseif p_No == 2
        star_CH = [1, 2, 3];
    elseif p_No == 3
        star_CH = [2, 3, 4];
    elseif p_No == 4
        star_CH = [3, 4];
    end
elseif size(collar_info, 2) == 5
    if p_No == 1
        star_CH = [1, 2];
    elseif p_No == 2
        star_CH = [1, 2, 3];
    elseif p_No == 3
        star_CH = [2, 3, 4];
    elseif p_No == 4
        star_CH = [3, 4, 5];
    elseif p_No == 5
        star_CH = [3, 4, 5];
    end
end
if isempty(star_CH)
    a=1;
end
num_inrange = 1;
list1 = -99999 * ones(20, 1);
if drct == -1 % search upward (current depth > last collar depth)
    for k = (i : 1 : size(collar_info, 1))
        if abs(depth(k) - depth_collar) <= (length+margin) &&...
                abs(depth(k) - depth_collar) >= (length-margin) && ...
                sum(collar_info(k, star_CH)) > 0
            list1(num_inrange) = depth(k);
            num_inrange = num_inrange + 1;
        elseif abs(depth(k) - depth_collar) < (length - margin)
            num_inrange = num_inrange - 1;
            break
        end
    end
elseif drct == 1 % search downward (current depth < last collar depth)
    for k = (i : -1 : 1)
        if abs(depth(k) - depth_collar) <= (length+margin) &&...
                abs(depth(k) - depth_collar) >= (length-margin) && ...
                sum(collar_info(k, star_CH)) > 0
            list1(num_inrange) = depth(k);
            num_inrange = num_inrange + 1;
        elseif abs(depth(k) - depth_collar) < (length - margin)
            num_inrange = num_inrange - 1;
            break
        end
    end
end

for h = 1:size(list1, 1)
    if list1(h, 1) ~= -99999
        list1(h, 2) = mod(abs(list1(h) - depth_collar), length);
    else
        list1(h, 2) = 999;
    end
end

list2 = sortrows(list1, 2);

if size(find(list2(:, 2) ~= 999), 1) > 1
    depth_1st = list2(1, 1);
    depth_2nd = list2(2, 1);
elseif size(find(list2(:, 2) ~= 999), 1) == 1
    depth_1st = list2(1, 1);
    depth_2nd = -99999;
elseif isempty(find(list2(:, 2) ~= 999))
    num_inrange = 0;
    depth_1st = -99999;
    depth_2nd = -99999;
end

end