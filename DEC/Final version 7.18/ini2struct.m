function ini_data = ini2struct(file_path)
    ini_data = struct();
    fid = fopen(file_path);
    current_section = '';
    while ~feof(fid)
        tline = strtrim(fgetl(fid));
        if isempty(tline) || startsWith(tline, ';')
            continue;
        elseif startsWith(tline, '[') && endsWith(tline, ']')
            current_section = tline(2:end-1);
            ini_data.(current_section) = struct();
        else
            kv = strsplit(tline, '=');
            if length(kv) == 2
                ini_data.(current_section).(strtrim(kv{1})) = strtrim(kv{2});
            end
        end
    end
    fclose(fid);
end