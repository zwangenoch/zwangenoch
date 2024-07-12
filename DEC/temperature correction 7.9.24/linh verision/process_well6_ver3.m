
%% Processing the well number 6. Version 2 is with Wang's angle correction method version 06/25/2024 (latest one)


% Load the file path
% file_path = 'D:\Linh\Commercial_jobs_DEC\Malaysia\6th_job\data\WELL NO.6 (EBA-07) DEC_JADESTONE (3.5in Tubing Interval)\WELL 6 (EBA-07) ORGINAL CSV FILES\WELL 6 (EBA-07) UpLog-2.csv';
file_path = 'D:\Linh\Commercial_jobs_DEC\Malaysia\6th_job\data\WELL NO.6 (EBA-07) DEC_JADESTONE (3.5in Tubing Interval)\WELL 6 (EBA-07) ORGINAL CSV FILES\WELL 6 (EBA-07) DownLog-2.csv';


[~, fileName, fileExt] = fileparts(file_path);
disp(fileName);

if (contains(fileName, 'UP') || contains(fileName, 'Up'))
    log_flag = 0;  % For uppass
elseif (contains(fileName, 'DN') || contains(fileName, 'DOWN') || contains(fileName, 'Down'))
    log_flag = 1;  % For downpass
end

Tool_num = 23916;
rawTable = readtable(file_path);
disp(height(rawTable))


%% Remove the row having more than 3 zero-columns 
% Convert table to an array for easier manipulation


% Count the number of zeros in each row
rowsWithZeros = any(rawTable{:, 1:48} == 0, 2);


% Identify rows with any column value less than -900
rowsWithNegativeValues = any(rawTable{:,:} < -900, 2);

% Combine both conditions
rowsToRemove = rowsWithZeros | rowsWithNegativeValues;

% Remove these rows from the table
rawTable(rowsToRemove, :) = [];
disp(height(rawTable))

temp1 = rawTable.('DE_DGTemp_degF');
%'DE_DGTemp_degF'  and 'DE_CCBTemp_degF'
% temp1 = rawTable.(selectColumns{1});
temp = temp1;
mean_value = mean(temp1);
std_value = std(temp1);

z_scores = (temp1 - mean_value) / std_value;
threshold = 3;
spike_indices = find(abs(z_scores) > threshold);

% Correct the temperature as mean of two neighbors
if ~isempty(spike_indices)
    for i = 1:length(spike_indices)
        if spike_indices(i) > 1 && spike_indices(i) < length(temp1)
            temp1(spike_indices(i)) = mean([temp1(spike_indices(i)-1), temp1(spike_indices(i)+1)]);
        end
    end

end
rawTable.('DE_DGTemp_degF') = temp1; % Reassign the corrected values


temp2 = rawTable.('DE_CCBTemp_degF');
%'DE_DGTemp_degF'  and 'DE_CCBTemp_degF'
% temp1 = rawTable.(selectColumns{1});
% temp = temp2;
mean_value = mean(temp2);
std_value = std(temp2);

z_scores = (temp2 - mean_value) / std_value;
threshold = 3;
spike_indices = find(abs(z_scores) > threshold);

% Correct the temperature as mean of two neighbors
if ~isempty(spike_indices)
    for i = 1:length(spike_indices)
        if spike_indices(i) > 1 && spike_indices(i) < length(temp2)
            temp2(spike_indices(i)) = mean([temp2(spike_indices(i)-1), temp2(spike_indices(i)+1)]);
        end
    end

end
rawTable.('DE_CCBTemp_degF') = temp2; % Reassign the corrected values

% % Step to sort the data in the descending order of Depth
columntoSort = "DEPTH_m";
sortedTable = sortrows(rawTable, columntoSort, "descend");

%% Continue with all the steps from latest version
raw_data = table2array(sortedTable);

%%%%% read from Warrior
BCas_raw = raw_data(:, 26:49)/(8*4096)*2.5; %% channel 6 is facing down in initial
BTub_raw = raw_data(:, 2:25)/(8*4096)*2.5;
Gyro_raw = raw_data(:,50:99);
ACC_raw = raw_data(:,100:105);
Dtemp = raw_data(:,106);
Btemp = raw_data(:,107);
Seq_num = raw_data(:,108);
Gyro_ofst = raw_data(:,109);
Depth = raw_data(:,1);

% Check if there are any NaN values
BTub_raw(isnan(BTub_raw(:,1)),:)=[];
BCas_raw(isnan(BCas_raw(:,1)),:)=[];
Depth(isnan(Depth(:,1)),:)=[];

% Use the tool responses as normalization factors
% tool_bias = median(BCas_raw);
% tool_bias_tub = median(BTub_raw);

load DEC21916_35_7_C_T_voltage.mat   % Should have tool_bias and tool_bias_tub by loading this file

BC_bdCr = BCas_raw - tool_bias + mean(tool_bias);
BT_bdCr = BTub_raw - tool_bias_tub + mean(tool_bias_tub);

BCas_raw = BC_bdCr;
BTub_raw = BT_bdCr ;
%% New folder destination
rotCrc_folder = 'D:\Linh\Commercial_jobs_DEC\Malaysia\5th_job\Gyro_revision_Jun252024';
current_folder = pwd;

cd(rotCrc_folder);

[ThetaDE,FullOrientation,...
    BC_rot, BT_rot,Depth_new] = Orientationmain_6_25_24( Depth, Gyro_raw, ACC_raw,...
    Dtemp, Btemp, Seq_num, BCas_raw,BTub_raw,Gyro_ofst,Tool_num, log_flag);
Warrio_rot_ang = mod(FullOrientation',360);

cd(current_folder);

Depth = Depth_new;
dep_len = length(Depth);
%% Pay attention here
idx_take = 1:15:360;
for dep_id = 1:dep_len

test = BC_rot(dep_id,idx_take);
bc_test = test/mean(test)-1;

test = BT_rot(dep_id,idx_take);
bc_test_tub = test/mean(test)-1;

test_fft = fft(bc_test');
C2_amp_rot(dep_id) = abs(test_fft(2,:));

test_fft_tub = fft(bc_test_tub');
T2_amp_rot(dep_id) = abs(test_fft_tub(2,:));
end

%% Part 2: Get from the BC_rot2
for dep_id = 1:dep_len   % Depth are the same for two methods: They have the same resolution 
    temp_part2 = BC_rot(dep_id,:);
    temp_sm2 = smooth_GOWELL(temp_part2,60);
    temp_sm2 = temp_sm2(121:240);
  [~,high_amp_pk(dep_id)] = max(temp_sm2);      % For each depth, having one highest peak for casing
end
high_amp_pk_sm = smooth_GOWELL(high_amp_pk,20);

%%%%%  test another method to find highest energy zone %%%%%
for dep_id = 1:dep_len
    temp_part = BT_rot(dep_id,:);    
    temp_sm = smooth_GOWELL(temp_part,60);
    temp_sm = temp_sm(121:240);
  [~,high_amp_tub_pk(dep_id)] = max(temp_sm);
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

high_amp_tub_pk_sm = smooth_GOWELL(high_amp_tub_pk,20);

%% tool position correction correction

tp_fct = 2;

C2_cos = cos(high_amp_pk_sm*3/180*pi);
C2_sin = sin(high_amp_pk_sm*3/180*pi);

T2_cos = cos((high_amp_tub_pk_sm*3)/180*pi);  %% when C2 T2 angle has an offset, + Angle, for Basker is 20
T2_sin = sin((high_amp_tub_pk_sm*3)/180*pi);  %% when C2 T2 angle has an offset, + Angle, for Basker is 20

ecc_tp_cos = C2_amp_rot.*C2_cos' - tp_fct*T2_amp_rot.*T2_cos';  %% ecc tool correction
ecc_tp_sin = C2_amp_rot.*C2_sin' - tp_fct*T2_amp_rot.*T2_sin';

ecc_tp_ang = angle(ecc_tp_cos + 1i* ecc_tp_sin)/pi*180;
ecc_tp_amp = abs(ecc_tp_cos + 1i* ecc_tp_sin);

% ecc_lab_rng = [0,0.3,0.6,1];a
ecc_lab_rng = [0.2,0.6,0.9,1];
% c2_amp_tp = [0.24,0.64,0.92,2.07];  %% from 23992 lab reading  4+9
c2_amp_tp = [0.8,0.9,1.1,1.2]*4.-3.1;  %% from 23992 lab reading  3+7

ecc_cnv = polyfit(c2_amp_tp,ecc_lab_rng,2);

ecc_factor = polyval(ecc_cnv,ecc_tp_amp);


ecc_factor_full = smooth_GOWELL(ecc_factor,50);
ecc_factor_full(ecc_factor_full < 0) = 0.01;      % This is after Karim's comment on some negative values of ecc_factor


tool_pst_rng = [0,2,5];
% c2_amp_tp = [0.24,0.64,0.92,2.07];  %% from 23992 lab reading 4+9
% t2_amp_tp = [0.17,0.6,1.3];  %% from 23991 lab reading 4+9
 t2_amp_tp = [0.02,0.9,1.9];  %% from 23991 lab reading 3+7
tp_cnv = polyfit(t2_amp_tp,tool_pst_rng,1);

% ecc_cnv(2) = ecc_cnv(2) + 0.25;

tool_pst = polyval(tp_cnv,T2_amp_rot);



tool_offset = smooth_GOWELL(tool_pst,60);

%% ECC VDL
idx_take_vdl = 1:3:360;
for dep_id = 1:length(Depth)
temp_c = BC_rot(dep_id,idx_take_vdl);
fft_temp = fft(temp_c);
amp_fct = ecc_tp_amp(dep_id)/abs(fft_temp(2));
fft_raw(2:2:length(fft_temp)) = amp_fct.*fft_temp(2:2:end);

ecc_vdl(dep_id,:) = ifft(fft_raw);

[~,max_id_ecc] = max(ecc_vdl(dep_id,:));

angle_offset = max_id_ecc - round(ecc_tp_ang(dep_id)/3);

% ecc_vdl(dep_id,:) = circshift(ecc_vdl(dep_id,:),angle_offset);

end
ecc_vdl = sign(real(ecc_vdl)).*abs(ecc_vdl);


%% updated ecc angle (Taking ecc_vdl as input)
[amp,ang_id]=max(ecc_vdl');
 
ecc_tp_ang = ang_id *3;
smooth_angle = nanfastsmoothAngle_shan(ecc_tp_ang/180*pi,20,3,0.1);
ecc_ang_full = smooth_angle'/pi*180;
ecc_ang_full(ecc_ang_full<0) = ecc_ang_full(ecc_ang_full<0)+360;


figure
plot(Depth,ecc_ang_full,'ro')
view(90,90)
grid on
box on
xlabel( 'Depth [m]')
ylabel ('ECC angle')
xlim([min(Depth),max(Depth)]);
ylim([0,360])


figure
plot(Depth,smooth_GOWELL(ecc_factor,50))            % 30 controls how much smooth of ECC factor. Should have some connection with the previous 50.
view(90,90)
grid on
box on
ylim([0,1])
xlabel( 'Depth [m]')
ylabel ('ECC factor')
xlim([min(Depth),max(Depth)]);


figure
imagesc(0:3:360,Depth,[ecc_vdl,ecc_vdl(:,1)])     %x, y and C: such that the Depth should be y axis
ylabel( 'Depth [m]')
xlabel( 'ECC VDL ')
clim([-0.03,.03]);
 colormap(jet);


figure
plot(Depth,smooth_GOWELL(tool_pst,60))
view(90,90)
grid on
box on
xlabel( 'Depth [m]')
ylabel( 'tool offset [mm]')
xlim([min(Depth),max(Depth)]);

% figure
% plot(Depth, high_amp_tub_pk_sm*3)
% view(90,90)
% grid on
% box on
% xlabel( 'Depth [m]')
% ylabel( 'tool offset angle [deg]')
% xlim([min(Depth),max(Depth)]);
% ylim([0,360])



figure
% plot(Depth,smooth(tool_pst,60))
plot(Depth,Warrio_rot_ang)
view(90,90)
grid on
box on
xlabel( 'Depth [m]')
ylabel( 'Tool rotation angle [deg]')
xlim([min(Depth),max(Depth)]);