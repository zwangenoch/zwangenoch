function [acc_calibration_para, acc_calibration_temperature, gyro_calibration_para, gyro_calibration_temperature, SC] = read_ini_file(file_path, instrument_name)
    ini_data = ini2struct(file_path);

    section_name = ['SN_', instrument_name];

    if ~isfield(ini_data, section_name)
        error('Section %s not found in ini file.', section_name);
    end

    section_data = ini_data.(section_name);

    acc_calibration_para = [];
    fieldnames_section = fieldnames(section_data);
    for i = 1:length(fieldnames_section)
        if startsWith(fieldnames_section{i}, 'Accel_cal_')
            line_content = section_data.(fieldnames_section{i});
            if endsWith(line_content, ',')
                line_content = line_content(1:end-1);
            end
            acc_calibration_para = [acc_calibration_para; str2double(strsplit(line_content, ','))];
        end
    end
    
    acc_calibration_temperature = str2double(strsplit(strtrim(section_data.Accel_calibT_1), ','));
    
    gyro_calibration_para = str2double(strsplit(strtrim(section_data.DigitCalib_Gyro_1), ','));
    
    gyro_calibration_temperature = str2double(strsplit(strtrim(section_data.Calib_GyroT_1), ','));
    
    SC = str2double(strtrim(section_data.SC_Val));
end