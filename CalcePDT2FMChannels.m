function [ v_chs ] = CalcePDT2FMChannels(n_p, od_p, tk_p, sig_p, mur_p, mode)

n_ch = max(mode.adec_idx_end);
v_chs = zeros(1, n_ch);

coils = 'SML';
for i = 1:3
    if i == 1
        msgbox("Calculating S Sensor");
    elseif i == 2
        msgbox("Calculating M Sensor");
    elseif i == 3
        msgbox("Calculating L Sensor");
    end
    if mode.sensor(i)
        v_ms = CalcePDTBSegMuAvgEMF(n_p, od_p, tk_p, coils(i), ...
            sig_p, mur_p, mode.txcurrent(i), mode.txduration(i), ...
            mode.chdelaytime, mode.chendtime(mode.adec_idx_end(i)));
        for j = mode.adec_idx_start(i):mode.adec_idx_end(i)
            v_chs(j) = mean(v_ms(...
                (mode.chstarttime(j) + 1):mode.chendtime(j)), 2);
        end
    end
end

end