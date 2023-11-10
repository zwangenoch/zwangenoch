function [x_ms, data_ms] = CHtoTime(data_CH, modeinfo, sensor_flag)

% sensor index will be according to sensorflag
% switch sensor_flag 
%     case 4
%         istart = 1;
%         iend = sum(modeinfo.sensor);
%     case 1
%         istart = 1;
%         iend = 1;
%     case 2
%         istart = 2;
%         iend = 2;
%     case 3
%         istart = 3;
%         iend = 3;
% end
istart = 1;
iend = sum(modeinfo.sensor);

data_ms = zeros(size(data_CH, 1), modeinfo.chendtime(modeinfo.adec_idx_end(iend)));
FM_ms = zeros(size(data_CH, 1), modeinfo.chendtime(modeinfo.adec_idx_end(iend)));

for i = istart:iend
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
        
        % correct y if y has a increase trend at the tail part
        diff_y = diff(y);
        idx_increase = find(diff_y>0, 1, 'first');
        if ~isempty(idx_increase)
            for m = idx_increase:size(y,2)
                y(m) = y(m-1) - (y(m-10)-y(m-1))/9;
            end
        end
        
        FM_ms((modeinfo.chstarttime(modeinfo.adec_idx_start(i))+1):(modeinfo.chendtime(modeinfo.adec_idx_end(i))))=y;
    else
        continue;
    end
end

data_ms = FM_ms;

x_ms = (1:modeinfo.chendtime(modeinfo.adec_idx_end(iend)));
end
