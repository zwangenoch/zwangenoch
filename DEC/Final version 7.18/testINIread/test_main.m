file_path = 'C:\Users\ze.wang\OneDrive - GOWell OilField Technology\Ze Wang\Github Code\zwangenoch\DEC\Final version 7.18\testINIread\GWDECTool_Cal 17.ini';
instrument_name = '23993'; 
[acc_calibration_para, acc_calibration_temperature, gyro_calibration_para, gyro_calibration_temperature, SC] = read_ini_file(file_path, instrument_name);
[acc_calibration_para1,acc_calibration_temperature1,gyro_calibration_temperature1,gyro_calibration_para1,SC1]=loadclibdat(str2double(instrument_name));

d1 = acc_calibration_para - acc_calibration_para1;
d2 = acc_calibration_temperature - acc_calibration_temperature1;
d3 = gyro_calibration_para - gyro_calibration_para1;
d4 = gyro_calibration_temperature - gyro_calibration_temperature1;
d5 = SC - SC1;
Ad = sum(d1, 'all') + sum(d2, 'all') + sum(d3, 'all') + sum(d4, 'all') + sum(d5, 'all');
a=1;