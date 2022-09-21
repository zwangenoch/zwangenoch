function v_c = CombineSegments(t_seg, v_seg, rx_time)
%CombineSegments  Combine segmented mu contributions to a single decay.
% Input:
%  t_seg        time from segmented mu, 1 x n_seg cell
%  v_seg        emf from segmented mu, 1 x n_seg cell
%  rx_time      acq total length, unit: ms
%
% Output:
%  v_c          combined decay

n_seg = length(t_seg);
y_seg = cell(1, n_seg);
d1x_seg = cell(1, n_seg);
d1y_seg = cell(1, n_seg);

% calculate derivatives
for i = 1:n_seg
    t = t_seg{i};
    y_seg{i} = log(v_seg{i}(t));
    d1x_seg{i} = t(1:(end - 1));
    d1y_seg{i} = diff(y_seg{i}); 
end

% curve fit
d1x_c = [d1x_seg{:}];
d1y_c = [d1y_seg{:}];
p3d = polyfit(d1x_c, d1y_c, 2);
d1y_s = polyval(p3d, 1:(rx_time - 1));

% compensate for constant slope from analytical modeling
d1y_max = d1y_s(d1y_s == max(d1y_s));
d1y_min = d1y_s(d1y_s == min(d1y_s));

% if gap < 2.5e-3 and max < -0.01 do correction
if (d1y_max - d1y_min) < 2.5e-3 && d1y_max < -0.1e-3
    % linear transform max to near zero point: -0.0075
    ratio = (-0.0075 - d1y_min) / (d1y_max - d1y_min);
    d1y_s = d1y_min + (d1y_s - d1y_min) * ratio;
end

% output decay
y_c = zeros(1, rx_time);
y_c(d1x_seg{1}(1)) = y_seg{1}(1);
for i = d1x_seg{1}(1):-1:2
    y_c(i - 1) = y_c(i) - d1y_s(i - 1);
end

for i = d1x_seg{1}(1):(rx_time - 1)
    y_c(i + 1) = y_c(i) + d1y_s(i);
end

v_c = exp(y_c);

end