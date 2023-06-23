function [EMDdata_1D, EMDdata_2D] = EMD_2D(vdl, depth, CHs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
param.n_grid_1    = 500;  %Grid size in dimension 1
param.nimfs       = 10;    %Maximum number of IMFs that can be stored
param.type        = 5; %type of window size
param.tol         = 0.05; %sifting tolerance
param.plot        = 'off'; %plots on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=depth;
Signal = vdl - 1;
Results = EMD1DNV(Signal,t',param);

data = Results.IMF(:, :, 1) + Results.IMF(:, :, 2);

% cut blow zero
data(data<0) = 0;

% denoise data veritcally
for k = 1:size(data, 2)
    data(:, k) = signal_denoise(data(:, k));
end

% calculate instantaneous energy for selected IMFs
EMDdata_1D = zeros(size(data, 1), size(CHs, 2));
for ch = 1:size(CHs, 2)
    [~,~,~,~,EMDdata_1D(:, ch)] = hht(data(:, CHs(ch)));
end

% smooth horizontally for plotting
for i = 1:size(data, 1)
    data(i, :) = signal_denoise(data(i, :));
end

EMDdata_2D = data;

end
