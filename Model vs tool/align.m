function [aligned_data, D, newCH, alignCH] = align(raw_data, nom_data, CH, testConfig)

% this function is to align the raw ADEC data with nominal data at selected CH
% D is for calculating k

% find the first CH to be aligned
diff = ceil(max(testConfig.adec_idx_end)/20);

if CH == testConfig.adec_idx_start(1) || CH == testConfig.adec_idx_start(2) ||...
        CH == testConfig.adec_idx_start(3)
    idx = CH;
    CH = CH + 2;
elseif CH < testConfig.adec_idx_start(1)
    msgbox("Picked CH must be greater than 1");
elseif CH > testConfig.adec_idx_start(1) && CH <= testConfig.adec_idx_end(1)
    if CH - testConfig.adec_idx_start(1) >= diff
        idx = CH - diff + 1;
    else
        idx = testConfig.adec_idx_start(1) + 1;
    end    
elseif CH > testConfig.adec_idx_start(2) && CH <= testConfig.adec_idx_end(2)
    if CH - testConfig.adec_idx_start(2) >= diff
        idx = CH - diff + 1;
    else
        idx = testConfig.adec_idx_start(2) + 1;
    end  
elseif CH > testConfig.adec_idx_start(3)
    if CH - testConfig.adec_idx_start(3) >= diff
        idx = CH - diff + 1;
    else
        idx = testConfig.adec_idx_start(3) + 1;
    end  
else
    msgbox("Picked CH cannot be greater than last CH");
end

% alignment in log scale
aligned_data(1:idx-1) = nom_data(1:idx-1);
for i = idx:max(testConfig.adec_idx_end)
    aligned_data(i) = raw_data(i) / (raw_data(idx) / nom_data(idx));    
end

% calculate D
D = abs(aligned_data(CH) - nom_data(CH));

alignCH = idx;
newCH = CH;
end