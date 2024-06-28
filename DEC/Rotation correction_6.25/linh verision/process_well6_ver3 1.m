
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
Warrior_rot_ang = mod(FullOrientation',360);

cd(current_folder);

