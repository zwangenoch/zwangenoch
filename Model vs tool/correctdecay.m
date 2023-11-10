function [decay_out] = correctdecay(decay, CHflag, logflag, config)

% this function is correct decay from logging or modeling by replace the
% later chs with straight line in log domain
% the steps are conducted in time domain

% mapping the data into time domain
if CHflag == 1
    [~, decay_t] = CHtoTime(decay, config, 4);
else
    decay_t = decay;
end

% find first value below 9
if logflag == 0
    decay_t = log10(decay_t);
end

% correct all available sensors
for i = 1:3
    if config.sensor(i) == 1
        idx = 0;
        idx = find(decay_t < log10(9),1, 'first');
        
        if  isempty(idx) || idx == 0
            idx = floor(4/5*config.chendtime(config.adec_idx_end(i)));
        end
        
        % find current sensor and compare with first ch of current sensor
        %sensor_No = find(config.chstarttime(config.adec_idx_end) >= idx, 1, 'last');
        if idx < (config.chendtime(config.adec_idx_start(i))+3)
            msgbox("Nominal value in bad quality, please re-select.");
            return
        end
        
        decay_corrected(1:idx-1) = decay_t(1:idx-1);
        slope = (decay_t(idx-3) - decay_t(idx-1))/2;
        for i = idx:size(decay_t, 2)
            decay_corrected(i) = decay_corrected(i-1) - slope;
        end
    end
end

% recover data if not logged
if logflag == 0
    decay_corrected = 10.^decay_corrected;
end

if CHflag == 1
    decay_out = zeros(1, max(config.adec_idx_end));
    for i = 1:3
        if i == 1 && config.sensor(i)
            for j = config.adec_idx_start(i):config.adec_idx_end(i)
                decay_out(j) = mean(decay_corrected(...
                    (config.chstarttime(j) + 1):config.chendtime(j)), 2);
            end
        elseif i == 2 && config.sensor(i)
            for j = config.adec_idx_start(i):config.adec_idx_end(i)
                decay_out(j) = mean(decay_corrected(...
                    (config.chstarttime(j) + 1):config.chendtime(j)), 2);
            end
        elseif i == 3 && config.sensor(i)
            for j = config.adec_idx_start(i):config.adec_idx_end(i)
                decay_out(j) = mean(decay_corrected(...
                    (config.chstarttime(j) + 1):config.chendtime(j)), 2);
            end
        end
    end
else
    decay_out = decay_corrected;
end

end