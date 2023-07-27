function [x_ms, data_ms] = CHtoTime(data_CH, modeinfo)

data_ms = zeros(size(data_CH, 1), modeinfo.chendtime(end));
FM_ms = zeros(size(data_CH, 1), modeinfo.chendtime(end));

for i = 1:3
    if modeinfo.sensor(i) == 1
        x_1 = zeros(1, (modeinfo.adec_idx_end(i) - modeinfo.adec_idx_start(i)+1));
        n = 1;
        for j = modeinfo.adec_idx_start(i):modeinfo.adec_idx_end(i)
            x_1(n) = ceil(mean([modeinfo.chstarttime(j), modeinfo.chendtime(j)]));
            n = n + 1;
        end
        y_1 = data_CH(modeinfo.adec_idx_start(i):modeinfo.adec_idx_end(i));
        p = polyfit(x_1, y_1, 10);
        x = (modeinfo.chstarttime(modeinfo.adec_idx_start(i))+1):(modeinfo.chendtime(modeinfo.adec_idx_end(i)));
        y = polyval(p, x);
        
        FM_ms((modeinfo.chstarttime(modeinfo.adec_idx_start(i))+1):(modeinfo.chendtime(modeinfo.adec_idx_end(i))))=y;
    else
        continue;
    end
end

data_ms = FM_ms;

x_ms = (1:modeinfo.chendtime(end));
end
