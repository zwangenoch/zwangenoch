function [emf_combined] = combine(emf, testConfig)

% this function is to combine segmented emf to a smooth emf

% S&M sensors
if testConfig.sensor(2) == 1
    emf_log = log(emf);
    emf_log(testConfig.adec_idx_start(2)) = emf_log(testConfig.adec_idx_start(2) - 1) - ...
        (emf_log(testConfig.adec_idx_start(2) - 2) - emf_log(testConfig.adec_idx_start(2) - 1));
    for i = (testConfig.adec_idx_start(2) + 1) : testConfig.adec_idx_end(2)
        emf_log(i) = log(emf(i)) - (log(emf(testConfig.adec_idx_start(2))) -...
            emf_log(testConfig.adec_idx_start(2)));
    end
end

% M&L sensors
if testConfig.sensor(3) == 1
    emf_log(testConfig.adec_idx_start(3)) = emf_log(testConfig.adec_idx_start(3) - 1) - ...
        (emf_log(testConfig.adec_idx_start(3) - 2) - emf_log(testConfig.adec_idx_start(3) - 1));
    for i = (testConfig.adec_idx_start(3) + 1) : testConfig.adec_idx_end(3)
        emf_log(i) = log(emf(i)) - (log(emf(testConfig.adec_idx_start(3))) -...
            emf_log(testConfig.adec_idx_start(3)));
    end
end

emf_combined = exp(emf_log);



end