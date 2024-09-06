clc; clear all;
% 20kHz
% ???????????????????
data = readtable('Data.xlsx');
depths = data{:,1}; % ???????? N
pressures = data{:,2}; % ??????????? N
R_test_values = linspace(0, 100, 1000); % ?????? R ????????????


% ???????????????
probabilities = zeros(length(R_test_values), length(depths)); % ???? R ?????

% ????????????
for i = 41:length(depths)-40
    % ????????????????
    d1 = depths(i-40) - depths(i); % d1 ????????????
    d2 = depths(i-20) - depths(i); % d2 ????????????
    d3 = 0; % ?????????
    d4 = depths(i+40) - depths(i); % d4 ????????????
    d5 = depths(i+20) - depths(i); % d5 ????????????

    % ??????
    P1 = 10^(pressures(i-40)/10);
    P2 = 10^(pressures(i-20)/10);
    P3 = 10^(pressures(i)/10); % ???????
    P4 = 10^(pressures(i+40)/10);
    P5 = 10^(pressures(i+20)/10);

    % ??????? R ?
    for j = 1:length(R_test_values)
        R_test = R_test_values(j); % ????? R ?
        
        % ????????????????
        cost1(j,i) = 1000 * sum([
            (P3/P1 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
            (exp(-0.01*sqrt(R_test^2 + d1^2)) / sqrt(R_test^2 + d1^2)))^2,%
            
%             (P3/P2 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d2^2)) / sqrt(R_test^2 + d2^2)))^2,
%             
%             (P3/P3 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d3^2)) / sqrt(R_test^2 + d3^2)))^2,
%             
%             (P3/P4 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d4^2)) / sqrt(R_test^2 + d4^2)))^2,
%             
%             (P3/P5 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d5^2)) / sqrt(R_test^2 + d5^2)))^2
        ]);
    
    cost2(j,i) = 1000 * sum([
%             (P3/P1 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d1^2)) / sqrt(R_test^2 + d1^2)))^2,%
            
            (P3/P2 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
            (exp(-0.01*sqrt(R_test^2 + d2^2)) / sqrt(R_test^2 + d2^2)))^2,
            
%             (P3/P3 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d3^2)) / sqrt(R_test^2 + d3^2)))^2,
%             
%             (P3/P4 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d4^2)) / sqrt(R_test^2 + d4^2)))^2,
%             
%             (P3/P5 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d5^2)) / sqrt(R_test^2 + d5^2)))^2
        ]);
    
    cost4(j,i) = 1000 * sum([
%             (P3/P1 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d1^2)) / sqrt(R_test^2 + d1^2)))^2,%
%             
%             (P3/P2 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d2^2)) / sqrt(R_test^2 + d2^2)))^2,
%             
%             (P3/P3 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d3^2)) / sqrt(R_test^2 + d3^2)))^2,
            
            (P3/P4 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
            (exp(-0.01*sqrt(R_test^2 + d4^2)) / sqrt(R_test^2 + d4^2)))^2,
            
%             (P3/P5 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d5^2)) / sqrt(R_test^2 + d5^2)))^2
        ]);
    
    cost5(j,i) = 1000 * sum([
%             (P3/P1 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d1^2)) / sqrt(R_test^2 + d1^2)))^2,%
%             
%             (P3/P2 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d2^2)) / sqrt(R_test^2 + d2^2)))^2,
%             
%             (P3/P3 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d3^2)) / sqrt(R_test^2 + d3^2)))^2,
%             
%             (P3/P4 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
%             (exp(-0.01*sqrt(R_test^2 + d4^2)) / sqrt(R_test^2 + d4^2)))^2,
            
            (P3/P5 - (exp(-0.01*sqrt(R_test^2)) / sqrt(R_test^2)) / ...
            (exp(-0.01*sqrt(R_test^2 + d5^2)) / sqrt(R_test^2 + d5^2)))^2
        ]);
        
        % ???????????????????????????
%         probabilities(j, i) = exp(-cost(j,i)); % ???????????cost ???????
        
    end
end
cst1 = log10(cost1);
cst2 = log10(cost2);
cst4 = log10(cost4);
cst5 = log10(cost5);
cst = cst1+cst2+cst4+cst5;
probability = 1./cst; 
% ?? R ???????
figure;
imagesc(depths, R_test_values, prob);
colorbar;
xlabel('Depth');
ylabel('R');
title('Probability Distribution of R at Different Depths');

