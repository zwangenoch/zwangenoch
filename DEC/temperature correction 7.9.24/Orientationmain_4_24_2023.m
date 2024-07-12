function [ThetaDE,FullOrientation,... %BCasing1, BTubing1,
    BCasinOriented, BTubinOriented,Depth_new] = Orientationmain_4_24_2023( depth, gyro_readings, Dataace22,...
    gyro_temperature, ccb_temp, seq_num, casing_original,tubing_original,realtime_angle_original2,Tool_num)
% depth(realtime_angle_original2==0)=[];
% gyro_readings(realtime_angle_original2==0,:)=[];
%  Dataace22(realtime_angle_original2==0,:)=[];
%     gyro_temperature(realtime_angle_original2==0)=[];
%  ccb_temp(realtime_angle_original2==0)=[];
% seq_num(realtime_angle_original2==0)=[];
% casing_original(realtime_angle_original2==0,:)=[];
% tubing_original(realtime_angle_original2==0,:)=[];
% realtime_angle_original2(realtime_angle_original2==0,:)=[];
% realtime_angle_original2=realtime_angle_original2-realtime_angle_original2(1);
%outputs 
%	ThetaDE_out output inclination
% 	FullOrientationDec Full rotation 
% inputs 
%CS Calibration data
%TS Temperature calibration data
%gyro_readings Gyro raw data 
%Dataace22 accelerometer raw data
% gyro_temperature Gyrotemperature sensor raw data
%ccb_temp CCb board temoperature raw data
%seq_num sequence number raw data
%
% [file,path] = uigetfile(".csv");
% DD=readmatrix (strcat(path,file));
% DD(1,:)=[];
% depth=DD(:,1);
% gyro_readings=DD(:,50:99);
% Dataace22=DD(:,100:105);
% gyro_temperature=DD(:,107);
% ccb_temp=DD(:,106);
% seq_num=DD(:,108);
% tubing_original=DD(:,2:25);
% casing_original=DD(:,26:49);
% realtime_angle_original2=DD(:,109);
TTDG2 = gyro_temperature;
GyroData2 = gyro_readings;
seq = seq_num;
% Tool_num=21945;
K=1;
[acc_calibration_para,acc_calibration_temperature,gyro_calibration_temperature,gyro_calibration_para,SC]=loadclibdat(Tool_num);
xi01        = acc_calibration_para(:, 1);
xi901       = acc_calibration_para(:, 2);
xi_901      = acc_calibration_para(:, 3);
yi01        = acc_calibration_para(:, 4);
yi901       = acc_calibration_para(:, 5);
yi_901      = acc_calibration_para(:, 6);
zi01        = acc_calibration_para(:, 7);
zi901       = acc_calibration_para(:, 8);
zi_901      = acc_calibration_para(:, 9);
x01         = acc_calibration_para(:, 10);
x901        = acc_calibration_para(:, 11);
x_901       = acc_calibration_para(:, 12);
y01         = acc_calibration_para(:, 13);
y901        = acc_calibration_para(:, 14);
y_901       = acc_calibration_para(:, 15);
z01         = acc_calibration_para(:, 16);
z901        = acc_calibration_para(:, 17);
z_901       = acc_calibration_para(:, 18);
CalibTAccel = acc_calibration_temperature;
CalibTGyro  = gyro_calibration_temperature;
Gyrocalib   = gyro_calibration_para;

AccData2 = Dataace22;
TTCCB2 = ccb_temp;
BTubingOrg2 = tubing_original;
BCasingOrg2 = casing_original;
SeqNum2 = seq_num;
Gyro(1,:) =GyroData2(1,:);   % Gyro data, total 25, sampling at 8ms
Dataace2(1,:) = AccData2(1,:);     % two accelerometer with 3 axis, sampling at 100ms
TGyro(1,:) = TTDG2(1,:);      % Temperature Data at Digital Gyro.
TAccel(1,:) = TTCCB2(1,:);      % Temperature Data at Controller Board.
seq(1,:) = SeqNum2(1,:); % The seqnum should be read from the dataset, here I generated the synthetic data to test the algorithm.              
k=1;
highside=0;
adc=0;
kk=0;
C1=0;
Cini=0;
M=0;
DF=0;
dt=0.25;
sigmax=50;
sigma=0.05;
F=1;
Q=dt^2*sigma^2 ;
H=1;
P=1;
gyro_angle(1)=0;
seqp=SeqNum2(1);
FullOrientation(1)=0;
offset(1)=0;
k=0;
drd=0;
for i=2:length(SeqNum2)
    if (SeqNum2(i)==SeqNum2(i-1))||(SeqNum2(i)<0)
    else 
        drd=drd+1;
    end 
end

Dataace2 = zeros(drd,6);     % two accelerometer with 3 axis, sampling at 100ms
TGyro= zeros(drd,1);      % Temperature Data at Digital Gyro.
TAccel= zeros(drd,1);      % Temperature Data at Controller Board.
seq =zeros(drd,1); % The seqnum should be read from the dataset, here I generated the synthetic data to test the algorithm.                     
BCasingZigzag=zeros(drd,24);
BTubingZigzag=zeros(drd,24);
BCas_original = zeros(drd,24);
BTub_original= zeros(drd,24);
Depth_new = zeros(drd,1);
realtime_angle_original=zeros(drd,1);
PHiDE=zeros(1,drd);
ThetaDE=zeros(1,drd);
FullOrientation=zeros(1,drd);
gyro_angle=zeros(1,drd);
offset=zeros(1,drd);
for i=2:length(SeqNum2)
    if (SeqNum2(i)==SeqNum2(i-1))||(SeqNum2(i)<0)
    else 
        k=k+1;
           % Gyro data, total 25, sampling at 8ms
        for jj=1:50
            if GyroData2(i,jj)>2^15 
                Gyro(k,jj) =GyroData2(i,jj)-2^16;
            else 
                Gyro(k,jj) =GyroData2(i,jj);
            end 
        end 
        Dataace2(k,:) = AccData2(i,:);     % two accelerometer with 3 axis, sampling at 100ms
        TGyro(k,:) = TTDG2(i,:);      % Temperature Data at Digital Gyro.
        TAccel(k,:)= TTCCB2(i,:);      % Temperature Data at Controller Board.
        seq(k,:) = SeqNum2(i,:); % The seqnum should be read from the dataset, here I generated the synthetic data to test the algorithm.                     
        BCas_original(k,:) = BCasingOrg2(i,:);
        BTub_original(k,:) = BTubingOrg2(i,:);
        Depth_new(k,:) = depth(i,:);
        realtime_angle_original(k,:)=realtime_angle_original2(i,:);
    end
end

%

% [file,path] = uigetfile(".csv");
file = 0;  %% skip load stationary process, Shan May12

% if file==0
%     Stationdata=[];
% else
%     DD1=readmatrix (strcat(path,file));
% gyro_readings=DD1(:,50:99);
for i=1:length(gyro_readings)
   for j=1:50
    if gyro_readings(i,j)>2^15
        gyro_readings(i,j)=gyro_readings(i,j)-2^16;
    end 
   end 
end 
Stationdata=mean(gyro_readings,'all');
Stationdata(2)=mean(gyro_temperature,'all');
% end

    
cas_mean_original = mean(casing_original(:));
tub_mean_original = mean(tubing_original(:));

%Stationdata=[];

for jj=1:6
    Dataace2(:,jj)= movmedian(Dataace2(:,jj),100);
end
for i=1:length(TGyro)
    %% accelerometer data perperation 
    for j=1:50
        if abs((median(Gyro(i,:))-Gyro(i,j)))>1e3
            Gyro(i,j)=median(Gyro(i,:));
        end
    end
end

%% revision by Ze Wang, Jun 26, 2024, flags for adjust highside angle difference when switching to ACC with previous using of ACC and gyro
% ori_diff is a array that has same dimension of zeros(drd,1)
% when the tool is using ACC->gyro->ACC for orientation, the ori_diff array
% will correct highside angle difference at switching
ori_diff = zeros(drd,1);
% 0: never using 
% 1: currently using
% 2: previous used, currently not using
ACC_ori_flag = 0; 
Gyro_ori_flag = 0;
status(1)=0;

%% revision by Ze Wang, function DEC1_accelerometer_phi_theta_5_21_2024 is used to calculate phi and theta
data_flag = 0;
%[GyroProCalib,CalibTGyroPro]=ProminanceAlgorithmupdates_UPDATE(Gyro,Dataace2,TGyro,TAccel,seq,SC,Gyrocalib,CalibTGyro,xi01,xi901,xi_901,yi01,yi901,yi_901,zi01,zi901,zi_901,x01,x901,x_901,y01,y901,y_901,z01,z901,z_901,CalibTAccel,Stationdata,data_flag,realtime_angle_original);
%%% comment out gyro calibration adjusting for stability, June 18. 24 by Ze Wang
[GyroProCalib,CalibTGyroPro]= DEC1_calibration_algorithm_gyroscope_6_25_24(Gyro,Dataace2,TGyro,TAccel,seq,SC,Gyrocalib,CalibTGyro,xi01,xi901,xi_901,yi01,yi901,yi_901,zi01,zi901,zi_901,x01,x901,x_901,y01,y901,y_901,z01,z901,z_901,CalibTAccel,Stationdata,data_flag,realtime_angle_original);
% GyroProCalib = Gyrocalib;
% CalibTGyroPro = CalibTGyro;  % keep gyro calibration temp and para as same as lab set
data_flag = 1;
for i=1:length(TGyro)
   if i==1
%         [PHiDE(i),ThetaDE(i)] = AccDE01_7_21_2021(Dataace2(i,:),TAccel(i),xi01,xi901,xi_901,yi01,yi901,yi_901,zi01,zi901,zi_901,x01,x901,x_901,y01,y901,y_901,z01,z901,z_901,CalibTAccel,data_flag);
% revised by Ze Jun 4. 2024, revised algorithm for angle calculation
        [PHiDE(i),ThetaDE(i)] = DEC1_accelerometer_phi_theta_5_21_2024(Dataace2(i,:),TAccel(i),x01,x901,x_901,y01,y901,y_901,z01,z901,z_901,CalibTAccel,data_flag);
   else 
%         [PHiDE(i),ThetaDE(i)] = AccDE01_7_21_2021(Dataace2(i,:),TAccel(i),xi01,xi901,xi_901,yi01,yi901,yi_901,zi01,zi901,zi_901,x01,x901,x_901,y01,y901,y_901,z01,z901,z_901,CalibTAccel,data_flag);
% revised by Ze Jun 4. 2024, revised algorithm for angle calculation
        [PHiDE(i),ThetaDE(i)] = DEC1_accelerometer_phi_theta_5_21_2024(Dataace2(i,:),TAccel(i),x01,x901,x_901,y01,y901,y_901,z01,z901,z_901,CalibTAccel,data_flag);
        [gyro_angle(i),offset(i),Gyrostat(i), modvalue(i)]=gyroalgorithm_6_25_24(Gyro(i,:),TGyro(i),seq(i),offset(i-1),seq(i-1),CalibTGyroPro,GyroProCalib,SC,data_flag,realtime_angle_original(i));
        % buffer 5-7 is added to avoid orientation jumping between gyro and ACC frequently
        if (Gyro_ori_flag ~= 1 && abs(ThetaDE(i))>5) || (Gyro_ori_flag == 1 && abs(ThetaDE(i))>7 && ACC_ori_flag==2)...
                || (Gyro_ori_flag == 1 && abs(ThetaDE(i))>5 && ACC_ori_flag==0)
            status(i)=1;
            if (PHiDE(i)-PHiDE(i-1))>180
                FullOrientation(i)=FullOrientation(i-1)+(PHiDE(i)-PHiDE(i-1))-360;
            elseif (PHiDE(i)-PHiDE(i-1))<-180
                FullOrientation(i)=FullOrientation(i-1)+(PHiDE(i)-PHiDE(i-1))+360;
            else 
                FullOrientation(i)=FullOrientation(i-1)+(PHiDE(i)-PHiDE(i-1));
            end
            if ACC_ori_flag==0
                ACC_ori_flag=1; % turn the currently using flag on
                highside=PHiDE(i)-mod(FullOrientation(i),360)+195; %Counter for high side 
                ori_diff(1:i) = 0; % all ori_diff no need to change before ACC first use
                if data_flag == 1
                    highside=-(-PHiDE(i)+mod(FullOrientation(i),360)-165); 
                end
            elseif ACC_ori_flag==1
                ori_diff(i) = ori_diff(i-1);
            elseif ACC_ori_flag==2 && Gyro_ori_flag==1 % meaning ACC was used before and reused after using gyro
                ACC_ori_flag=1; % turn the currently using flag on
                ori_diff(i)=(PHiDE(i)-mod(FullOrientation(i),360)+195)-highside; %compensate difference of highside change
                if data_flag == 1
                    ori_diff(i)=-(-PHiDE(i)+mod(FullOrientation(i),360)-165)-highside; 
                end
                if ori_diff(i) < -180
                    ori_diff(i) = ori_diff(i) + 360; % avoid ori_diff<-180 because FullOrientation is never less than -180
                end
                % adjust previous gyro angle data
                ori_diff_grad = (ori_diff(i)-ori_diff(i_gyro_start))/(i - i_gyro_start);
                for n = i_gyro_start:(i-1)
                    ori_diff(n) = ori_diff(n) + (n-i_gyro_start+1) * ori_diff_grad;
                end   
            end                   
            % reset gyro ori status
            if Gyro_ori_flag==1
                Gyro_ori_flag = 2;
            end
        else 
            status(i)=2;
            Gyro_ori_flag = 1;
            ori_diff(i) = ori_diff(i-1);
            if ACC_ori_flag==1
                ACC_ori_flag = 2;
                i_gyro_start = i;
            end
            FullOrientation(i)=(FullOrientation(i-1)+(gyro_angle(i)-gyro_angle(i-1)));
        end
    % Gyroscope data perparation
   end 
end
 
FullOrientation=FullOrientation + highside + transpose(ori_diff); % use ori_diff to compensate highside rotation

%% revision of high resolution to 1 degree and rotation shift in this section
% BCasinOriented=zeros(size(BCas_original));
% BTubinOriented=zeros(size(BCas_original));

[nrows, mcols] = size(BCas_original);
BCas_originalex= zeros(nrows, mcols + 10);
BTub_originalex= zeros(nrows, mcols + 10);

% add five columns to both sides of the orignal data
for icol= 1:5
    BCas_originalex(:,icol)= BCas_original(:,(mcols - icol) + 1);
    BTub_originalex(:,icol)= BTub_original(:,(mcols - icol) + 1);
end
BCas_originalex(:,icol + 1:mcols + 5)= BCas_original;
BTub_originalex(:,icol + 1:mcols + 5)= BTub_original;
for icol= 1:5
    BCas_originalex(:,(mcols + 5) + icol)= BCas_original(:,icol);
    BTub_originalex(:,(mcols + 5) + icol)= BTub_original(:,icol);
end

sample_points = 1:34;
query_points = 1:1/15:34;

cubic_interpolator_casing = zeros(nrows, (mcols + 10 - 1) * 15 + 1);
cubic_interpolator_tubing = zeros(nrows, (mcols + 10 - 1) * 15 + 1);

for i = 1:nrows
    cubic_interpolator_casing(i,:) = interp1(sample_points, BCas_originalex(i,:), query_points, 'cubic');
    cubic_interpolator_tubing(i,:) = interp1(sample_points, BTub_originalex(i,:), query_points, 'cubic');
end

% drop extra columns
BCas_originalexsuper_spline = cubic_interpolator_casing(:, 76:end-61); 
BTub_originalexsuper_spline = cubic_interpolator_tubing(:, 76:end-61); 
% prepare for rotation correction
BCasinOriented = zeros(size(BCas_originalexsuper_spline));
BTubinOriented = zeros(size(BCas_originalexsuper_spline));

% rotation correction
for i = 1:length(FullOrientation)
    ShiftNum = round(FullOrientation(i) / (15 / 15));
    if data_flag == 1
        BCasinOriented(i,:) = circshift(BCas_originalexsuper_spline(i,:), -ShiftNum);  
        BTubinOriented(i,:) = circshift(BTub_originalexsuper_spline(i,:), -ShiftNum);
    else
        BCasinOriented(i,:) = circshift(BCas_originalexsuper_spline(i,:), ShiftNum);  
        BTubinOriented(i,:) = circshift(BTub_originalexsuper_spline(i,:), ShiftNum);
    end
end

%         for i=1:length(FullOrientation)
%             ShiftNum = round((FullOrientation(i)/15)); %added
%             if data_flag == 1
%                 BCasinOriented(i,:) = circshift(BCas_original(i,:), -ShiftNum);  %added
%                 BTubinOriented(i,:) = circshift(BTub_original(i,:), -ShiftNum); %added
%             else
%                 BCasinOriented(i,:) = circshift(BCas_original(i,:), ShiftNum);  %added
%                 BTubinOriented(i,:) = circshift(BTub_original(i,:), ShiftNum); %added
%             end
%         end
 
end