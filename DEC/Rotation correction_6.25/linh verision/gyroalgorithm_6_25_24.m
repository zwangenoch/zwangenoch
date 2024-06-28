function [gyro_angle,gyro_anglep1]=gyroalgorithm_7_21_2021(DGyro,TTr,seq,gyro_anglep,seqp,calibT,calib,SC,data_flag,realtime_angle_original)
    %Raw data
% output: gyro_angle: anlge measured after 
% input :DGyro Gyro raw data 
%	TTr Gyro temperautre 
%	seq sequence number of the current frame 
%	gyro_anglep anlge estimated at the previous frame 
%	seqp sequence number of the previous frame 
%	calibT Temperature of calibration 
%	calib calibration 
%	SC Scaling factor 
%	data_flag direction of the logging 
     %coefficient to adjsut for each tool depending on the calibration 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt=0.005;
    %Important coefficient for Kalman
    %change the matrix of the data to a single cloumn 
    % processing using Kalman filter 
    hhhd=1;
    jj=1;
    i=1;
    %Intial temperature calibration 
    while (hhhd)&&(jj<length(calibT))
        jj=jj+1;
        if TTr(i)>calibT(jj)
            hhhd=1;
        else 
            hhhd=0;
        end 
    end
    Gyrostat2=calib(jj);
    Gyrostat1=calib(jj-1);
    T2=calibT(jj);
    T1=calibT(jj-1);
    Gyrostat=(Gyrostat2-Gyrostat1)*(TTr(i)-T1)/(T2-T1)+Gyrostat1;% warning
    %X(1+50*(i-1):50*i)=X(1+50*(i-1):50*i)-Gyrostat;
    
    if data_flag == 1
%         x=mod((seq-seqp),255);% need to remove
        gyro_anglep1=gyro_anglep-Gyrostat*dt*SC*mod((seq-seqp),255);%((seqp-seq)+255); %mod(-(seqp-seq),255);
        gyro_angle=-realtime_angle_original+gyro_anglep1;
    else
%         x=mod((seq-seqp),255);% need to remove
        gyro_anglep1=gyro_anglep-Gyrostat*dt*SC*mod((seq-seqp),255);%(seq-seqp); %mod((seq-seqp),255);
        gyro_angle=realtime_angle_original-gyro_anglep1;
    end
%     x=x; % need to remove
end 
