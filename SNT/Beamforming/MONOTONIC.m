x = linspace(0, 10000, 100); % ??x???
y = (0.1 - (exp(-0.1 * sqrt(x.^2)) ./ sqrt(x.^2)) ./ ...
            (exp(-0.1 * sqrt(x.^2 + 5^2)) ./ sqrt(x.^2 + 5^2))).^2;

% ??????????
dy = diff(y) ./ diff(x);

% ??????
if all(dy > 0)
    disp('1');
elseif all(dy < 0)
    disp('2');
else
    disp('3');
end