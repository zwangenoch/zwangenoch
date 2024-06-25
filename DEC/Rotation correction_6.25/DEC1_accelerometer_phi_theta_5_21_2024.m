function [PHiDE,ThetaDE,theta1,phi1]= DEC1_accelerometer_phi_theta_5_21_2024(resampled_accel_readings, resampled_accel_temp, ...
                               xi01, xi901, xi_901, yi01, yi901, yi_901, zi01, zi901, zi_901, ...
                               CalibTAccel, data_flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The script computes the displacement using accelerometer raw
% readings
%
% The generated outputs  are:
% - the inclination angle (ThetaDE)
% - the rotation angle (PHiDE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input variables 
% - resampled accelerometer raw readings (resampled_accel_readings)
% - resampled accelerometer temperature (resampled_accel_temp)
% - calibration parameters for accelerometer, refer to main script (GyroCalib, 
%   CalibTGyro, xi01, xi901, xi_901, yi01, yi901, yi_901, zi01, zi901, zi_901, x01, x901, 
%   x_901, y01, y901, y_901, z01, z901, z_901, CalibTAccel)
% - logging direction (data_flag)

% get first reading from the second accelerometer chip
% 1st column is ACC1-z, 2nd col is ACC1-x, 3rd is ACC1-y
FFx= resampled_accel_readings(1,3);
FFy= resampled_accel_readings(1,2);
FFz= resampled_accel_readings(1,1);

% loop through calibration temperature of accelerometer
% to get the current and previous values of index 'jj'
% when accelerometer temperature reading is not higher than
% calibration temperature of accelerometer
hhhd= 1;
jj= 1;

while hhhd && jj < length(CalibTAccel)
    jj= jj + 1;

    if resampled_accel_temp > CalibTAccel(jj)
        hhhd= 1;
    else
        % condition to break the while loop
        hhhd= 0;
    end
end
% get calibration parameters for the second accelerometer chip
% along three axes based on indexes 'jj' and 'jj-1'
% index 'jj' --> denotes accelerometer temperature < calibration
% temperature
% index 'jj-1' --> denotes accelerometer temperature > calibration
% temperature
xxi0= xi01(jj);
xxi02= xi01(jj - 1);

xxi90= xi901(jj);
xxi902= xi901(jj - 1);

xxi_90= xi_901(jj);
xxi_902= xi_901(jj - 1);

zzi0= zi01(jj);
zzi02= zi01(jj - 1);

zzi90= zi901(jj);
zzi902= zi901(jj - 1);

zzi_90= zi_901(jj);
zzi_902= zi_901(jj - 1);

yyi0= yi01(jj);
yyi02= yi01(jj - 1);

yyi90= yi901(jj);
yyi902= yi901(jj - 1);

yyi_90= yi_901(jj);
yyi_902= yi_901(jj - 1);

T2= CalibTAccel(jj);
T1= CalibTAccel(jj - 1);

% combine calibration temperature assuming a
% piecewise interpolation for the second accelerometer chip
xi0= (xxi0 - xxi02) * (resampled_accel_temp - T1)/(T2 - T1) + xxi02;
yi0= (yyi0 - yyi02) * (resampled_accel_temp - T1)/(T2 - T1) + yyi02;
zi0= (zzi0 - zzi02) * (resampled_accel_temp - T1)/(T2 - T1) + zzi02;

% combine calibration temperature assuming a
% piecewise interpolation for the second accelerometer chip
xi90= (xxi90 - xxi902) * (resampled_accel_temp - T1)/(T2 - T1) + xxi902;
yi90= (yyi90 - yyi902) * (resampled_accel_temp - T1)/(T2 - T1) + yyi902;
zi90= (zzi90 - zzi902) * (resampled_accel_temp - T1)/(T2 - T1) + zzi902;

% combine calibration temperature assuming a
% piecewise interpolation for the second accelerometer chip
xi_90= (xxi_90 - xxi_902) * (resampled_accel_temp - T1)/(T2 - T1) + xxi_902;
yi_90= (yyi_90 - yyi_902) * (resampled_accel_temp - T1)/(T2 - T1) + yyi_902;
zi_90= (zzi_90 - zzi_902) * (resampled_accel_temp - T1)/(T2 - T1) + zzi_902;

%%%Same PArt as last time
%%Nothing changed in this part

% get first reading from the second accelerometer chip
xin= FFx;
yin= FFy;
zin= FFz;

% the x-axis of the second acceletemer chip is linearly calibrated
if FFx > xi0
    xaccel= (xin - xi0)/(xi90 - xi0); 
    
% or by the line created by [-g,0]
else
    xaccel= -(xin -  xi0)/(xi_90 - xi0);
end
% if xaccel is greater than 1, round it to 1
% because sometime the calibration abs max or min is less than
% measurement abs
if xaccel > 1
    xaccel = 1;
elseif xaccel < -1
    xaccel = -1;
end

% the y-axis of the second acceletemer chip is linearly calibrated
if FFy > yi0
    yaccel= (yin - yi0)/(yi90 - yi0);
    
% or by the line created by [-g,0]
else
    yaccel= -(yin - yi0)/(yi_90 - yi0);
end
% if yaccel is greater than 1, round it to 1
if yaccel > 1
    yaccel = 1;
elseif yaccel < -1
    yaccel = -1;
end

% the z-axis of the second acceletemer chip is linearly calibrated
if FFz > zi0
    zaccel= (zin - zi0)/(zi90 - zi0);
    
% or by the line created by [-g,0]
else
    zaccel= -(zin - zi0)/(zi_90 - zi0);
end
% if zaccel is greater than 1, round it to 1
if zaccel > 1
    zaccel = 1;
elseif zaccel < -1
    zaccel = -1;
end

% the second accelerometer chip is tilted 45 degrees
% this causes all three axes have a vertical component.
% in the extreme cases: (1) logging tool is horizontal, 
% the x-axis has the strongest vertical component, (2)
% logging tool is vertical the y-axis has the strongest
% vertical component
% only the x- and y-axes corrected
xaccel1= xaccel;
yaccel1= yaccel;

% compute inclination angle (ThetaDE)
theta= atan(((yaccel1^2 + zaccel^2))^0.5/xaccel1)*180/pi;

%theta = acos(xaccel1/1)*180/pi;
% 

% use lookup table to do arctan
theta1 = search_in_table(((yaccel1^2 + zaccel^2))^0.5, xaccel1);

if xaccel1 < 0
    theta= theta + 180;
end
% if xaccel1 < 0
%     theta1= theta1 + 180;
% end
% compute rotation angle (PhiDE)
if zaccel >= 0
    phi= atan(yaccel1/zaccel)*180/pi;    
    % use lookup table to do arctan
    phi1 = search_in_table(yaccel1, zaccel);

elseif yaccel1 >= 0
    phi= atan(yaccel1/zaccel)*180/pi + 180;
    phi1 = search_in_table(yaccel1, zaccel);

else
    phi= atan(yaccel1/zaccel)*180/pi - 180;
    phi1 = search_in_table(yaccel1, zaccel);
end

PHiDE= phi;

if data_flag == 1
    PHiDE = -PHiDE;
    phi1 = -phi1;
end

ThetaDE= theta;

end
