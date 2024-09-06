clc; clear all;
% 20kHz
% ?????P?d?
P_all = [];
d_all = [];
for m = 1:256
P_set = movmean(P_all(:,m), 120);
d_set = d_all;
for i = 1:90
    P_prime1 = P_set(i);
    P_prime2 = P_set(i+90);
    d_prime1 = d_set(i);
    d_prime2 = d_set(i+90);%(566-i)
    P1 = -90.707/P_prime1; % ???P?
    d1_1 = 0; 
    d1_2 = abs(d_prime1-9570.46);
    P2 = -90.707/P_prime2; % ???P?
    d2_1 = 0; 
    d2_2 = abs(d_prime2-9570.46);

%% ?????
% ??????
objective = @(params) sum([
    ((exp(-params(2) * sqrt(params(1)^2 + d1_1^2)) / sqrt(params(1)^2 + d1_1^2)) / ...
    (exp(-params(2) * sqrt(params(1)^2 + d2_1^2)) / sqrt(params(1)^2 + d2_1^2)) - P1)^2,
    
    ((exp(-params(2) * sqrt(params(1)^2 + d1_2^2)) / sqrt(params(1)^2 + d1_2^2)) / ...
    (exp(-params(2) * sqrt(params(1)^2 + d2_2^2)) / sqrt(params(1)^2 + d2_2^2)) - P2)^2
]);

% ?????
D0 = 1; % ?D?????
alpha0 = 0.1; % ?alpha?????

% ??fminsearch????????
params0 = [D0, alpha0];
options = optimset('Display', 'iter');
params_solution = fminsearch(objective, params0, options);

% ???
R_solution = params_solution(1);
alpha_solution = params_solution(2);

% ????
% disp(['D = ', num2str(D_solution)]);
% disp(['alpha = ', num2str(alpha_solution)]);
R_set(i) = R_solution;
alpha_set(i) = alpha_solution;
end
alphap(:,m)=transpose(alpha_set);
Rp (:,m)= transpose(R_set);
end

%% ?????
% % ?????
% fun = @(params) [
%     (exp(-params(2) * sqrt(params(1)^2 + d1_1^2)) / sqrt(params(1)^2 + d1_1^2)) / ...
%     (exp(-params(2) * sqrt(params(1)^2 + d2_1^2)) / sqrt(params(1)^2 + d2_1^2)) - P1;
%     (exp(-params(2) * sqrt(params(1)^2 + d1_2^2)) / sqrt(params(1)^2 + d1_2^2)) / ...
%     (exp(-params(2) * sqrt(params(1)^2 + d2_2^2)) / sqrt(params(1)^2 + d2_2^2)) - P2
% ];
% 
% % ????
% R0 = 1; % ?D?????
% alpha0 = 0.1; % ?alpha?????
% 
% % ??fsolve??
% params0 = [R0, alpha0];
% options = optimoptions('fsolve', 'Display', 'iter');
% params_solution = fsolve(fun, params0, options);
% 
% % ???
% R_solution = params_solution(1);
% alpha_solution = params_solution(2);
% 
% % ????
% disp(['R = ', num2str(D_solution)]);
% disp(['alpha = ', num2str(alpha_solution)]);