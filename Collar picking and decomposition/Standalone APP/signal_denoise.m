function x_denoise=signal_denoise(x)
%% Moving average filter
% x is a matrix and filter x along rows
% windowSize =10; 
%  b = (1/windowSize)*ones(1,windowSize);
% x_denoise= (filter(b,1,x'))';
% data_old_denoise(:,1:windowSize)=data_old(:, 1:windowSize);

%% Savitzky-Golay finite impulse response (FIR) smoothing filter
r=2; % Polynomial order, must be smaller than framelen
f=29; % Frame length, specified as a positive odd integer.
x_denoise= (sgolayfilt(x', r, f))';

%% Outlier removal using Hampel identifier
% x_denoise=hampel(x',1);
% x_denoise=x_denoise'; 
end