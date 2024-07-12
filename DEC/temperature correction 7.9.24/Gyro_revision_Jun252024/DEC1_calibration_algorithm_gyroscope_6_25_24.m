function [GyroProCalib,CalibTGyroPro]= DEC1_calibration_algorithm_gyroscope_11_3_2023(resampled_gyro_readings, ...
                                resampled_accel_readings, resampled_gyro_temp, resampled_accel_temp, ...
                                resampled_seq, SC, Gyrocalib, CalibTGyro, xi01, xi901, xi_901, ...
                                yi01, yi901, yi_901, zi01, zi901, zi_901, x01, x901, x_901, y01, y901, ...
                                y_901, z01, z901, z_901, CalibTAccel, Stationdata, data_flag, ...
                                resampled_realtime_angle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The script computes the calibration for the drift in the
% gyroscope raw readings
%
% The generated outputs  are:
% - the dynamic calibration of gyroscope (GyroProCalib)
% - the dynamic calibration temperature (CalibTGyroPro)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input variables 
% - resampled gyroscope raw readings (resampled_gyro_readings)
% - resampled accelerometer raw readings (resampled_accel_readings)
% - resampled gyroscope temperature (resampled_gyro_temp)
% - resampled accelerometer temperature (resampled_accel_temp)
% - resample sequence number of raw readings(resample_seq)
% - scaling coefficient for gyroscope (SC)
% - calibration parameters for gyroscope and accelerometer, refer to main script (Gyrocalib, 
%   CalibTGyro, xi01, xi901, xi_901, yi01, yi901, yi_901, zi01, zi901, zi_901, x01, x901, 
%   x_901, y01, y901, y_901, z01, z901, z_901, CalibTAccel)
% - addtional calibration for gyroscope (Stationdata)
% - logging direction (data_flag)

%%%%%%%%%%%%%%Start of the Dynamic calibration: Basic C Variance%%%%%%%%%
% C is also used in warrior for realtime processing

% set data_flag to 0, because currently we have not found the correct
% way for downpass to get prominence
data_flag = 0;

% initialization block
Ko= 0;
CC= [];

% loop through gyroscope temperature readings
for i= 81:length(resampled_gyro_temp)
    % extract a sub-matrix of 80 rows x all columns 
    % from gyroscope readings to convert it into a vector
    X= reshape(resampled_gyro_readings(i-80:i,:), [], 1);
    
    % check the variance of the extracted readings
    if (var(X) < 15)
        Ko= Ko + 1;

        % compute mean from extracted gyroscope data
        CC(Ko,1)= mean(X);

        % compute mean from 80 data points from
        % gyroscope temperature readings
        CC(Ko,2)= mean(resampled_gyro_temp(i-80:i));
    end
end
% f=figure('Name', 'accel temp','Visible', 'off');
% % plot(CC(:,1))
% legend('Average Gyro readings')
% saveas(f,'DEC1_prominence_average_gyro_readings.png')
% f=figure('Name', 'accel temp','Visible', 'off');
% plot(CC(:,2))
% legend('Average Gyro temp')
% saveas(f,'DEC1_prominence_average_gyro_temperature.png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% New prominance Algorithm%%%%%%%%%%%%%%%%%%%%%

% apply a series of filters to remove normal
% and (potentially) impulse noises from the first accelerometer
% chip 
AccData(:,1)= movmedian(resampled_accel_readings(:,1),50);
AccData(:,2)= movmedian(resampled_accel_readings(:,2),50);

AccData(:,1)= movmean(resampled_accel_readings(:,1),100);
AccData(:,2)= movmean(resampled_accel_readings(:,2),100);

ACCx= smoothdata(AccData(:,1),'gaussian',50);
ACCy= smoothdata(AccData(:,2),'gaussian',50);
% f=figure('Name', 'accel temp','Visible', 'off');
% plot(ACCx, 'k')
% hold on
% plot(ACCy, 'r')
% legend('ACCx smooth no noise','ACCy smooth no noise')
% saveas(f,'DEC1_prominence_accx-accy.png')

% find prominence peaks/throughs and their locations
[pkxs,locxs]= findpeaks(ACCx,'MinPeakProminence',5);
[pkxsn,locxsn]= findpeaks(-ACCx,'MinPeakProminence',5);   
[pkys,locys]= findpeaks(ACCy,'MinPeakProminence',5);
[pkysn,locysn]= findpeaks(-ACCy,'MinPeakProminence',5);

% count how many peaks/throughs were found
L1= length(locxs);
L2= length(locxsn);
L3= length(locys);
L4= length(locysn);
L= max([L1,L2,L3,L4]);

% generate a matrix with 'inf' values. Then replace the 'inf' value
% by the index of the 4 prominence peaks
MM= inf(L,4);
MM(1:length(locxs),1)= locxs;
MM(1:length(locxsn),3)= locxsn;
MM(1:length(locys),2)= locys;
MM(1:length(locysn),4)= locysn;

% generate plots with found prominence peaks
f=figure('Visible','off');
plot(ACCx)
hold on 
plot(locxs,pkxs,'r+')
plot(locxsn,-pkxsn,'k+')
hold off
legend('ACCx')
saveas(f,'DEC1_prominence_peaks_accx.png')
f=figure('Visible','off');
plot(ACCy)
hold on 
plot(locys,pkys,'r+')
plot(locysn,-pkysn,'k+')
hold off
legend('ACCy')
saveas(f,'DEC1_prominence_peaks_accy.png')

% convert matrix into a vector such that the vector
% consists of a sequence of 4 elements extracted 
% row by row from the source matrix
X= reshape(MM',[],1)';

% loop through prominence peaks to get their location and their respective index location
LL= length(X);
X1= [];
X2= [];

for i= 1:LL
    % find the smallest prominence peak/through location and its index location
    [DD,ind]= min(X);

    % stop loop if all prominence peak/through locations are replaced by 'inf'
    if (DD == inf)
        break
    end 

    % get the smallest prominence peak/through value
    X1(i)= DD;

    % compute the modulus of the associated index location
    % to assign values between [0, 1, 2, 3]
    X2(i)= mod(ind,4);

    % replace smallest prominence peak location by 'inf'
    X(ind)= inf;
end 

% search for 4 indexes of prominence peak locations representing
% maximum sensitivity to acceleration
k= 0;
i= 1;
closest_prominence_peaks= [];

while (i <= length(X2) - 4) 
    % make use of 4 consecutive index locations 'i'
    % to get the computed modulus of the prominence's index location
    FF= [X2(i) X2(i + 1) X2(i + 2) X2(i + 3)];
    
    % and get the separation between these locations
    DF= X2(i) - X2(i + 1) + X2(i + 2) - X2(i + 3);
    
    % look for a group of 4 consecutive index locations
    % that have (1) unique modulus values and (2) they are 2 units apart
    if (length(unique(FF)) == 4) && abs(DF) == 2
        % get the index locations of the closest 4 peaks/throughs
        k= k + 1;
        closest_prominence_peaks(k,:)= X1(i:i + 4);

        % evaluate if the extracted prominence's index locations
        % correspond to clockwise or anticlockwise rotation
        % by computing the slope
        if mod(X2(i) - X2(i+1), 4) == 3
            % clockwise rotation was found
            closest_prominence_peaks(k,5)= 1;

        else 
            % anticlockwise rotation was found
            closest_prominence_peaks(k,5)= -1;
        end 

        i= i + 4;

    else 
        i= i + 1;
    end 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Old Accelermeter C1
% C1 is also used in warrior for real time processing

% compute the average row wise from gyroscope readings
LLdd= mean(resampled_gyro_readings,2);

% initialization block
CN= [];
CAcc= [];
DDVa= 0;
kk= 0;
ii= 0;

gyro_angle= resampled_realtime_angle;

FullOrientation= zeros(length(LLdd), 1);
Seqelap= zeros(length(LLdd), 1);
PHiDE= zeros(length(LLdd), 1);
ThetaDE= zeros(length(LLdd), 1);

% compute angles PhiDE (spinning) and ThetaDE (inclination) from accelerometer data
[PHiDE(1),ThetaDE(1)]= DEC1_accelerometer_phi_theta_5_21_2024(resampled_accel_readings(1,:), resampled_accel_temp(1), ...
                               xi01, xi901, xi_901, yi01, yi901, yi_901, zi01, zi901, zi_901, ...
                               CalibTAccel, data_flag);
% [PHiDE(1),ThetaDE(1)]= AccDE01_7_21_2021(resampled_accel_readings(1,:), resampled_accel_temp(1),xi01,xi901,xi_901,yi01,yi901,yi_901,zi01,zi901,zi_901,x01,x901,x_901,y01,y901,y_901,z01,z901,z_901,CalibTAccel,data_flag);

% loop through gyroscope readings
for i= 2:length(LLdd)
    % the modulus reduces the separation between current and 
    % previous sequence to (largely) 1 or a higher value. How higher
    % the difference, it depends on how apart are the consecutive 
    % sequence numbers. The final outcome is a new sequence

    % direction of the logging: top to bottom
    if data_flag == 1
        Seqelap(i) = Seqelap(i-1) + mod(resampled_seq(i-1) - resampled_seq(i), 255);
    else
        Seqelap(i) = Seqelap(i-1) + mod(resampled_seq(i) - resampled_seq(i-1), 255);
    end

    % compute angles PhiDE (rotation) and ThetaDE (inclination) from accelerometer data
    [PHiDE(i),ThetaDE(i)]= DEC1_accelerometer_phi_theta_5_21_2024(resampled_accel_readings(i,:), resampled_accel_temp(i), ...
                               xi01, xi901, xi_901, yi01, yi901, yi_901, zi01, zi901, zi_901, ...
                               CalibTAccel, data_flag);
%     [PHiDE(i),ThetaDE(i)]= AccDE01_7_21_2021(resampled_accel_readings(i,:), resampled_accel_temp(i),xi01,xi901,xi_901,yi01,yi901,yi_901,zi01,zi901,zi_901,x01,x901,x_901,y01,y901,y_901,z01,z901,z_901,CalibTAccel,data_flag);
    % check if the wellbore is vertical or inclined/horizontal
    if (abs(ThetaDE(i))>5)
        % when wellbore is inclined the full rotation is measured
        % with accelerometer data

        % check if full rotation is over +/- 180 and apply
        % correction if necessary
        if (PHiDE(i) - PHiDE(i-1)) > 180
            FullOrientation(i)= FullOrientation(i - 1) + (PHiDE(i) - PHiDE(i - 1)) - 360;
        
        elseif (PHiDE(i) - PHiDE(i - 1)) < -180
            FullOrientation(i)= FullOrientation(i - 1) + (PHiDE(i) - PHiDE(i - 1)) + 360;
        
        else 
            FullOrientation(i)= FullOrientation(i - 1) + (PHiDE(i) - PHiDE(i - 1));
        end 

        % direction of the logging: top to bottom
        % compute variable used to avoid an impulse noise because logging tool
        % is inclined
        if data_flag == 1
            DDVa= DDVa + mod(resampled_seq(i-1) - resampled_seq(i), 255);

        else
            DDVa= DDVa + mod(resampled_seq(i) - resampled_seq(i-1),255);
        end
        
        kk= kk + 1;

        % in case the computed variable is higher than ~100 seconds
        if DDVa > 400
            ii= ii + 1;

            % compute correction from accelerometer to gyroscope
            % NOTE: the FullOrientation is divided by a factor that
            % represents how often data is taken from gyroscope (5 msec)
            CAcc(ii,1)= (((FullOrientation(i) - FullOrientation(i - kk))) - ...
                        ((gyro_angle(i) - gyro_angle(i - kk))))/DDVa/0.005/SC;

            CAcc(ii,2)= mean(resampled_gyro_temp(i - kk:i));
            
            % reset variables
            DDVa= 0;
            kk= 0;
        end 
    else 
        % when wellbore is vertical the full rotation is measured
        % with gyroscope data
        kk= 0;
        DDVa= 0;
        FullOrientation(i)= FullOrientation(i - 1) + (gyro_angle(i) - gyro_angle(i - 1));
    end
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Continuation of the prominance Algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop through retained prominence peak locations to compute the
% offset correction
nrows= size(closest_prominence_peaks, 1);

CN= zeros(nrows, 13);

for i=1:nrows 
    % assuming a unit circle: LL1 ... LL4 correspond
    % how far are the peaks/throughs from each other
    % following the idea how sine/cosine function fluctuates
    % North (LL1), East (LL2), South (LL3), West (LL4)
    LL1= gyro_angle(closest_prominence_peaks(i,1));
    LL2= gyro_angle(closest_prominence_peaks(i,2));
    LL3= gyro_angle(closest_prominence_peaks(i,3));
    LL4= gyro_angle(closest_prominence_peaks(i,4));

    % direction of the logging: top to bottom
    if data_flag == 0
        % based on estimated rotation orientation
        % if closest_prominence_peaks(i,5) is 1 --> clockwise direction, angle is
        % substracted
        % if closest_prominence_peaks(i,5) is -1 --> anticlockwise direction, angle is
        % added
        % CNN --> normalized by sequence, constant, configuration value
        CN(i,1)= ((LL4 - LL2) + 180 * -closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,4)) - Seqelap(closest_prominence_peaks(i,2)))/0.005/SC;
        CN(i,2)= ((LL3 - LL1) + 180 * -closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,3)) - Seqelap(closest_prominence_peaks(i,1)))/0.005/SC;
        CN(i,3)= ((LL4 - LL1) + 270 * -closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,4)) - Seqelap(closest_prominence_peaks(i,1)))/0.005/SC;
        CN(i,4)= ((LL2 - LL1) + 90 * -closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,2)) - Seqelap(closest_prominence_peaks(i,1)))/0.005/SC;
        CN(i,5)= ((LL3 - LL2) + 90 * -closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,3)) - Seqelap(closest_prominence_peaks(i,2)))/0.005/SC;
        CN(i,6)= ((LL4 - LL3) + 90 * -closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,4)) - Seqelap(closest_prominence_peaks(i,3)))/0.005/SC;
    else
        % based on estimated rotation orientation
        % if closest_prominence_peaks(i,5) is 1 --> clockwise direction, then angle is
        % substracted
        % if closest_prominence_peaks(i,5) is -1 --> anticlockwise direction, angle is
        % added
        % CNN --> normalized by sequence, constant, scaling factor
        CN(i,1)= ((LL4 - LL2) + 180 * closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,4)) - Seqelap(closest_prominence_peaks(i,2)))/0.005/SC;
        CN(i,2)= ((LL3 - LL1) + 180 * closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,3)) - Seqelap(closest_prominence_peaks(i,1)))/0.005/SC;
        CN(i,3)= ((LL4 - LL1) + 270 * closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,4)) - Seqelap(closest_prominence_peaks(i,1)))/0.005/SC;
        CN(i,4)= ((LL2 - LL1) + 90 * closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,2)) - Seqelap(closest_prominence_peaks(i,1)))/0.005/SC;
        CN(i,5)= ((LL3 - LL2) + 90 * closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,3)) - Seqelap(closest_prominence_peaks(i,2)))/0.005/SC;
        CN(i,6)= ((LL4 - LL3) + 90 * closest_prominence_peaks(i,5))/ ...
                 (Seqelap(closest_prominence_peaks(i,4)) - Seqelap(closest_prominence_peaks(i,3)))/0.005/SC;
    end
    if i < nrows
        % assuming a unit circle: LL5 ... LL8 correspond
        % how far are the peaks/throughs from each other
        % following the idea how sine/cosine function fluctuates
        % North (LL5), East (LL6), South (LL7), West (LL8)
        % with respect to previous positions, the new ones are on top
        % of each other
        LL5= gyro_angle(closest_prominence_peaks(i + 1,1));
        LL6= gyro_angle(closest_prominence_peaks(i + 1,2));
        LL7= gyro_angle(closest_prominence_peaks(i + 1,3));
        LL8= gyro_angle(closest_prominence_peaks(i + 1,4));

        % direction of the logging: top to bottom
        if data_flag == 0
            % based on estimated rotation orientation
            % if closest_prominence_peaks(i,5) is 1 --> clockwise direction, angle is
            % substracted
            % if closest_prominence_peaks(i,5) is -1 --> anticlockwise direction, angle is
            % added
            % CNN --> normalized by sequence, constant, scaling factor
            CN(i,7)= ((LL5 - LL1) + 360 * -closest_prominence_peaks(i,5))/ ...
                     (Seqelap(closest_prominence_peaks(i + 1,1)) - Seqelap(closest_prominence_peaks(i,1)))/0.005/SC;
            CN(i,8)= ((LL6 - LL2) + 360 * -closest_prominence_peaks(i,5))/ ...
                     (Seqelap(closest_prominence_peaks(i + 1,2)) - Seqelap(closest_prominence_peaks(i,2)))/0.005/SC;
            CN(i,9)= ((LL7 - LL3) + 360 * -closest_prominence_peaks(i,5))/ ...
                     (Seqelap(closest_prominence_peaks(i + 1,3)) - Seqelap(closest_prominence_peaks(i,3)))/0.005/SC;
            CN(i,10)= ((LL8 - LL4) + 360 * -closest_prominence_peaks(i,5))/ ...
                      (Seqelap(closest_prominence_peaks(i + 1,4)) - Seqelap(closest_prominence_peaks(i,4)))/0.005/SC;

        else
            % based on estimated rotation orientation
            % if closest_prominence_peaks(i,5) is 1 --> clockwise direction, angle is
            % substracted
            % if closest_prominence_peaks(i,5) is -1 --> anticlockwise direction, angle is
            % added
            % CNN --> normalized by sequence, constant, scaling factor
            CN(i,7)= ((LL5 - LL1) + 360 * closest_prominence_peaks(i,5))/ ...
                     (Seqelap(closest_prominence_peaks(i + 1,1)) - Seqelap(closest_prominence_peaks(i,1)))/0.005/SC;
            CN(i,8)= ((LL6 - LL2) + 360 * closest_prominence_peaks(i,5))/ ...
                     (Seqelap(closest_prominence_peaks(i + 1,2)) - Seqelap(closest_prominence_peaks(i,2)))/0.005/SC;
            CN(i,9)= ((LL7 - LL3) + 360 * closest_prominence_peaks(i,5))/ ...
                     (Seqelap(closest_prominence_peaks(i + 1,3)) - Seqelap(closest_prominence_peaks(i,3)))/0.005/SC;
            CN(i,10)= ((LL8 - LL4) + 360 * closest_prominence_peaks(i,5))/ ...
                      (Seqelap(closest_prominence_peaks(i + 1,4)) - Seqelap(closest_prominence_peaks(i,4)))/0.005/SC;
        end

        % get the 1st and 4th index locations of prominence peaks
        CN(i,11)= closest_prominence_peaks(i,1);
        CN(i,12)= closest_prominence_peaks(i,4);

        % get the mean of gyroscope temperature readings based on the 1st and 
        % 4th index locations of prominence peaks
        CN(i,13)= mean(resampled_gyro_temp(closest_prominence_peaks(i,1):closest_prominence_peaks(i,4)));
    else 
         % make sure that the CN matrix is full
         CN(i,7:10)= 0;
         CN(i,11)= closest_prominence_peaks(i,1);
         CN(i,12)= closest_prominence_peaks(i,4);
         CN(i,13)= mean(resampled_gyro_temp(closest_prominence_peaks(i,1):closest_prominence_peaks(i,4)));
    end     
end 

% loop through computed offset correction to find how close are
% to each other in order to have an idea of the drift and then 
% to decide to apply a correction or not
[LF,LG]= size(CN);
CNF= [];
LG= LG - 3;
dd= 0;

for i= 1:LF
    k= 0;
    for j= 1:4
        % from columns 7th to 10th compute the mean on current row index 'i' 
        % and compare it against a value on current row index 'i' from
        % column 7th to 10th
        if abs(CN(i,LG-4+j) - mean(CN(i,LG-3:LG))) < 2
            k= k + 1;
        else 
            k= 0;
        end 
    end

    kk= 0;
    ddk= 0;

    for j= 1:6
        % from columns 7th to 10th compute the mean on current row index 'i'
        % and compare it against a value on current row index 'i' from
        % column 1st to 6th
        if abs(CN(i,j) - mean(CN(i,LG-3:LG))) < 2
            kk= kk + 1; 
        end 

        % from columns 1st to 6th compute the median on current row index 'i'
        % and compare it against a value on current row index 'i' from
        % column 1st to 6th
        if abs(CN(i,j) - median(CN(i,1:6))) < 1.5
            ddk= ddk + 1;
        end 
    end
    
    if k == 4 && kk >= 1 && i < LF
        % for current row index 'i' 
        % --> the corrections for CN[columns 7th to 10th] worked on 
        % all columns
        % --> the corrections for CN[columns 1sth to 6th] worked at least
        % on one column
        % CNF: (a) get average of realtime gyro angle, (b-c) get 1st and 4th
        % locations of peaks/throughs, (d) get average of gyroscope
        % temperature readings based on the 1st and 4th locations of
        % peaks/throughs
        dd= dd + 1;
        CNF(dd,1)= mean(CN(i,LG-3:LG));
        %CNF(dd,2)= CN(i,11);
        %CNF(dd,3)= CN(i,12);
        CNF(dd,4)= CN(i,13);

    elseif ddk >= 6
        % for current row index 'i' 
        % --> the corrections for CN[columns 1st to 6th] worked on 
        % all columns
        % CNF: (a) get median of realtime gyro angle, (b-c) get 1st and 4th
        % locations of peaks/throughs, (d) get average of gyroscope
        % temperature readings based on the 1st and 4th locations of
        % peaks/throughs
        dd= dd + 1;
        CNF(dd,1)= median(CN(i,1:6));
        %CNF(dd,2)= CN(i,11);
        %CNF(dd,3)= CN(i,12);
        CNF(dd,4)= CN(i,13);
    end 
end 

% %%% Revised by Ze Wang Jun19.2024 due to DL can have large portion of outlier causing problem on calib table
% % use lab calibration parameter as a standard to remove outlier in all components
% % get max and min limit from Gyrocalib
% if max(Gyrocalib) > 0
%     maxlimit = 3*max(Gyrocalib);
% else
%     maxlimit = 80; % empirical number
% end
% if min(Gyrocalib) < 0
%     minlimit = 3*min(Gyrocalib);
% else
%     minlimit = -80;
% end
% 
% % remove outlier in all components
% if ~isempty(CC)
%     idx = (CC(:, 1) > maxlimit) + (CC(:, 1) < minlimit);
%     if ~isempty(idx)
%         idx = idx>0;
%         CC(idx, :) = [];
%     end
% end
% if ~isempty(CAcc)
%     idx = (CAcc(:, 1) > maxlimit) + (CAcc(:, 1) < minlimit);
%     if ~isempty(idx)
%         idx = idx>0;
%         CAcc(idx, :) = [];
%     end
% end
% if ~isempty(CNF)
%     idx = (CNF(:, 1) > maxlimit) + (CNF(:, 1) < minlimit);
%     if ~isempty(idx)
%         idx = idx>0;
%         CNF(idx, :) = [];
%     end
% end
% if ~isempty(Stationdata)
%     idx = (Stationdata(:, 1) > maxlimit) + (Stationdata(:, 1) < minlimit);
%     if ~isempty(idx)
%         idx = idx>0;
%         Stationdata(idx, :) = [];
%     end
% end

% combine computed corrections to remove the drift in the gyroscope data
if isempty(CNF)
    DL= [CC;CAcc]; %Stationdata

else
    DL= [CC;-CNF(:,1) CNF(:,4);CAcc]; %Stationdata
end

% sort corrections according to temperature readings
FF= (sortrows(DL,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Preparing the new Calib%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get unique values according to temperature readings
DD= unique(FF(:,2));

% loop through unique gyroscope temperature readings
DS= zeros(size(DD,1),2);

for i= 1:length(DD)
    % get the average of computed corrections based on unique
    % temperature readings
    DS(i,1)= mean(FF(DD == DD(i),1));

    % get the unique temperature reading
    DS(i,2)= DD(i);
end

% initialization block
FinDS= size(DS, 1);

% correct the calibrated data of gyroscope with the help of gyroscope 
% readings and temperature
if FinDS == 0   
    % DS variable is empty
    CalibTGyroPro= CalibTGyro;
    GyroProCalib= Gyrocalib;
    
elseif FinDS == 1
    % DS variable only contains a single data point

    % initialization block
    hhhd= 1;
    jj= 1;

    % loop through calibration gyroscope temperature
    while hhhd && jj < length(CalibTGyro)
        jj= jj + 1;

        % check temperature reading in DS
        if DS(1,2) > CalibTGyro(jj)
            hhhd= 1;

        else 
            % condition to break the while loop
            hhhd= 0;
        end 
    end

    % compute calibrated data and temperature for gyroscope
    Gyrostat2= Gyrocalib(jj);
    Gyrostat1= Gyrocalib(jj - 1);
    T2= CalibTGyro(jj);
    TGyro= CalibTGyro(jj - 1);

    Gyrostat= (Gyrostat2 - Gyrostat1) * (DS(1,2) - TGyro)/(T2 - TGyro) + Gyrostat1;
    
    % get the corrected gyroscope data
    GyroProCalib= Gyrocalib - Gyrostat+DS(1,1);
    CalibTGyroPro= CalibTGyro;
    
else
    % DS variable contains several single data points

    % loop through calibration gyroscope temperature
    hhhd= 1;
    jj= 1;

    while hhhd && jj < length(CalibTGyro)
        jj= jj + 1;

        % check temperature reading in FF
        if FF(1,2) > CalibTGyro(jj)
            hhhd= 1;

        else 
            % condition to break the while loop
            hhhd= 0;
        end 
    end

    % compute calibrated data and temperature for gyroscope
    % first section
    Gyrostat2= Gyrocalib(jj);
    Gyrostat1= Gyrocalib(jj - 1);
    T2= CalibTGyro(jj);
    TGyro= CalibTGyro(jj - 1);
    T1forcalib1= TGyro;
    
    Gyrostat= (Gyrostat2 - Gyrostat1) * (DS(1,2) - TGyro)/(T2 - TGyro) + Gyrostat1;
    
    calib2min= -Gyrostat + DS(1,1);    

    % loop through calibration gyroscope temperature
    hhhd= 1;
    jj= 1;

    while hhhd && jj<length(CalibTGyro)
        jj= jj + 1;

        % check temperature reading in FF
        if DS(end,2) > CalibTGyro(jj)
            hhhd= 1;

        else 
            % condition to break the while loop
            hhhd= 0;
        end 
    end

    % compute calibrated data and temperature for gyroscope
    % second section
    Gyrostat2= Gyrocalib(jj);
    Gyrostat1= Gyrocalib(jj - 1);
    T2= CalibTGyro(jj);
    TGyro= CalibTGyro(jj - 1);
    T1forcalib2= T2;

    Gyrostat=(Gyrostat2 - Gyrostat1) * (DS(end,2) - TGyro)/(T2 - TGyro) + Gyrostat1;
    
    calib2max=-Gyrostat+DS(end,1); 
    
    GyroCalib(CalibTGyro<=T1forcalib1)= Gyrocalib(CalibTGyro <= T1forcalib1) + calib2min;
    GyroCalib(CalibTGyro>T1forcalib2)= Gyrocalib(CalibTGyro > T1forcalib2) + calib2max;
%     GyroCalib= Gyrocalib + calib2max;

    % get the corrected gyroscope data
    GyroProCalib= [GyroCalib(CalibTGyro < DS(1,2)) DS(:,1)' GyroCalib(CalibTGyro > DS(end,2))];
    CalibTGyroPro= [CalibTGyro(CalibTGyro < DS(1,2)) DS(:,2)' CalibTGyro(CalibTGyro > DS(end,2))];
end
%%% revised by Ze Wang Jun19.2024
% use moving median to smooth GyroProCalib because calibration table cannot
% sudden change which against physics
GyroProCalib = movmedian(GyroProCalib, 5);

% f= figure('Visible', 'off');
% plot(GyroProCalib,'k')
% legend('GyroProCalib')
% saveas(f,'DEC1__prominence_calibgyropro.png')
% f= figure('Visible', 'off');
% plot(CalibTGyroPro,'k')
% legend('CalibTGyroPro')
% saveas(f,'DEC1_prominence_calibtempgyropro.png')
end
