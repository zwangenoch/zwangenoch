function [acc_calibration_para_new]=shrinkrange(acc_calibration_para, acc_calibration_temperature)

table = [1,0.995380000000000,0.990710000000000,0.986040000000000,0.981370000000000,0.976700000000000,0.972030000000000,0.967360000000000,0.962690000000000,0.958020000000000,0.953350000000000,0.948680000000000,0.944010000000000,0.939340000000000,0.934670000000000,0.930000000000000
    ];
acc_calibration_para_new=acc_calibration_para;
for i = 1:size(acc_calibration_temperature, 2)
    acc_calibration_para_new(i, 2) = acc_calibration_para(i, 1)+table(i) * (acc_calibration_para(i, 2) - acc_calibration_para(i, 3))/2;
    acc_calibration_para_new(i, 3) = acc_calibration_para(i, 1)-table(i) * (acc_calibration_para(i, 2) - acc_calibration_para(i, 3))/2;
    acc_calibration_para_new(i, 5) = acc_calibration_para(i, 4)+table(i) * (acc_calibration_para(i, 5) - acc_calibration_para(i, 6))/2;
    acc_calibration_para_new(i, 6) = acc_calibration_para(i, 4)-table(i) * (acc_calibration_para(i, 5) - acc_calibration_para(i, 6))/2;
    acc_calibration_para_new(i, 8) = acc_calibration_para(i, 7)+table(i) * (acc_calibration_para(i, 8) - acc_calibration_para(i, 9))/2;
    acc_calibration_para_new(i, 9) = acc_calibration_para(i, 7)-table(i) * (acc_calibration_para(i, 8) - acc_calibration_para(i, 9))/2;
end

end