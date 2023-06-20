function dist = ...
    calc_mode(depth1, depth2, length)
% function to calculate the distance between two different depths with
% specific length of a joint

if abs(depth1 - depth2) < length
    dist = abs(abs(depth1 - depth2) - length);
elseif abs(depth1 - depth2) == length
    dist = 0;
elseif abs(depth1 - depth2) > length
    dist = abs(abs(depth1 - depth2) - floor(abs(depth1 - depth2)/length) * length);
end