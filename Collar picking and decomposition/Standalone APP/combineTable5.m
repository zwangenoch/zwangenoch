function [table_out] = combineTable5(tp1, tp2, tp3, a_achr, range, len_spec, depthrange)

% this function is to combine and output table considering all three collar
% tables on pipe 5 based on up-to-3 given anchor depth

a1 = tp1(:, 2);
a2 = tp2(:, 2);
a3 = tp3(:, 2);

% replace a1 a2 a3 with a_achr if it is 0
if a2 == 0
    a2(1) = a_achr(1);
end
if a3 == 0
    a3(1) = a_achr(1);
end

a4 = a_achr;
% find repeating values
for i = 1:length(a1)
    if max(abs(a2 - a1(i)) <= range) || max(abs(a3 - a1(i)) <= range)
        a4 = [a4; a1(i)];
    end
end
for j = 1:length(a2)
    if max(abs(a3 - a2(j)) <= range) 
        a4 = [a4; a2(j)];
    end
end
a4 = unique(a4);

% insert anchor points to a1 2 3
a1=[a1; a_achr];
a2=[a2; a_achr];
a3=[a3; a_achr];

% generate a new array by filling Table 4
a4 = sort(a4);
a41 = a4;
k = 1;
diff = 0;
while k <= (length(a4) - 1)
    diff = abs(a4(k) - a4(k+1));
    n_intervals = round(diff / len_spec);
    if (diff / n_intervals) <= (1.4*len_spec) || (diff / n_intervals) >= (0.6*len_spec)        
        if n_intervals > 1
            n_points = n_intervals - 1;
            points = transpose(linspace(a4(k), a4(k+1), n_points+2));
            a41 = [a41; points(2:end-1)];
        end
        k = k + 1;
    else
        a41(a41==a4(k+1)) = [];
        diff = abs(a4(k) - a4(k+2));
        n_intervals = round(diff / len_spec);
        if n_intervals > 1
            n_points = n_intervals - 1;
            points = transpose(linspace(a4(k), a4(k+1), n_points+2));
            a41 = [a41; points(2:end-1)];
        end
        k = k + 2;
    end
end

% create table 5 with repeat depth
a5 = a_achr;
for m = 1:length(a41)
    if any(abs(a1 - a41(m)) <= range) || any(abs(a2 - a41(m)) <= range) || any(abs(a3 - a41(m)) <= range)
        a5 = [a5; a41(m)];
    end
end
a5 = unique(a5);

j = 1;
diff = 0;
a51 = a5;
while j < (length(a5) - 1)
    diff = abs(a5(j) - a5(j+1));
    n_intervals = round(diff / len_spec);
    if (diff / n_intervals) <= (1.4*len_spec) || (diff / n_intervals) >= (0.6*len_spec)
        n_intervals = round(diff / len_spec);
        if n_intervals > 1
            n_points = n_intervals - 1;
            points = transpose(linspace(a5(j), a5(j+1), n_points+2));
            a51 = [a51; points(2:end-1)];
        end
        j = j + 1;
    else
        a51(a51==a5(j+1)) = [];
        diff = abs(a5(j) - a5(j+2));
        n_intervals = round(diff / len_spec);
        if n_intervals > 1
            n_points = n_intervals - 1;
            points = transpose(linspace(a5(j), a5(j+1), n_points+2));
            a51 = [a51; points(2:end-1)];
        end
        j = j + 2;
    end
end

a_final = a51;
if a51(1) > depthrange(1)
    diff = abs(a51(1) - depthrange(1));
    n_intervals = round(diff / len_spec);
        if n_intervals > 1
            n_points = n_intervals - 1;
            points = transpose(linspace(a51(1), depthrange(1), n_points+2));
            a_final = [a51; points(2:end-1)];
        end
end
if a51(end) < depthrange(2)
    diff = abs(max(a51) - depthrange(2));
    n_intervals = round(diff / len_spec);
        if n_intervals > 1
            n_points = n_intervals - 1;
            points = transpose(linspace(max(a51), depthrange(2), n_points+2));
            a_final = [a_final; points(2:end-1)];
        end
end

a_final = unique(a_final);

table_out = zeros(size(a_final));
table_out(:,2) = sort(a_final);
table_out(:, 1) = 1:size(a_final, 1);
table_out(:, 3) = 2;

for n = 1:size(table_out, 1)
    if round(table_out(n, 2)) == round(a_achr(1))
        table_out(n, 3) = 1;
    end
end
if length(a_achr) > 1
    for n = 1:size(table_out, 1)
        if round(table_out(n, 2)) == round(a_achr(2))
            table_out(n, 3) = 1;
        end
    end
end
if length(a_achr) > 1
    for n = 1:size(table_out, 1)
        if round(table_out(n, 2)) == round(a_achr(3))
            table_out(n, 3) = 1;
        end
    end
end
if max(any(table_out == 0))
    [~, idx] = min(table_out(:, 2));
    table_out(idx, :) = [];

end
end
% 
% clc; clear all;
% a1=[57.365
% 90.035
% 128.95
% 167.365
% 206.2
% 244.78
% 283.78
% 321.7
% 361.035
% 400.03
% 438.865
% 478.03
% 516.195
% 555.865
% 595.865
% 635.865
% 675.865
% 715.86
% 755.865
% 795.865
% 829.945
% 871.11
% 910.11
% 950.695
% ];
% a2=[57.365
% 97.365
% 137.365
% 177.37
% 217.365
% 257.365
% 297.37
% 337.365
% 385.45
% 423.115
% 463.115
% 503.115
% 543.115
% 592.615
% 631.865
% 671.615
% 711.11
% 750.78
% 788.78
% 828.365
% 865.695
% 904.61
% 943.36
% ];
% a3=[34.53
% 76.7
% 116.2
% 155.865
% 195.865
% 234.615
% 273.615
% 314.28
% 354.03
% 393.285
% 433.45
% 472.45
% 512.035
% 551.865
% 592.53
% 632.285
% 671.365
% 711.28
% 750.945
% 788.11
% 830.535
% 868.11
% 909.95
% 950.36
% ];
% depthrange=[1, 1000];
% len_spec = 39;
% a_achr = [234.615; 385.45; 829.945];
% range = 3.0;



