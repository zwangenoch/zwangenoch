clear;
close all; clc;

% this is a signal decomposition process using ensembeled EMD. 
% you may set a break point at vdl=vdl-1; and past 5 chs around the 
% target ch into vdl and do following calculation

% input selected vdl (2 channels on each side of target CH, total 5)
vdl=1;
vdl = vdl-1;

IMF_compress = zeros(size(vdl));
imfinse = zeros(size(vdl));

for CH = 1:size(vdl, 2)     
    
    Kurtosis = zeros(1, size(vdl, 2));
    % emd for each CH
    [imf, residual] = emd(vdl(:, CH));
    
    % calculate Kurtosis for each IMF of this CH
    for n = 1:size(imf, 2)
        avg_IMF(1, n) = mean(imf(:, n), 1);
        std_IMF(1, n) = std(imf(:, n), 1);
        for m = 1:size(imf, 1)
            diff_IMF(m, n) = (imf(m, n) - avg_IMF(1, n))^4;
        end
        avg_diff(1, n) = mean(diff_IMF(:, n), 1);
        Kurtosis(n) = avg_diff(n)/(std_IMF(n)^4);
    end
    
    % select proper IMF to sum according to Kuritosis
    index1 = 1; 
    index2 = 1; 
    index3 = 1;
    Kurt_temp = Kurtosis(3:end);
    [~, index1] = max(Kurt_temp);
    Kurt_temp(index1) = -Inf;
    [~, index2] = max(Kurt_temp);
    Kurt_temp(index2) = -Inf;
    [~, index3] = max(Kurt_temp);
    Kurt_temp(index3) = -Inf;
    index1 = index1 + 2;
    index2 = index2 + 2;
    index3 = index3 + 2;
    index1 = 3; 
    index2 = 4;

    if (Kurtosis(index1) + Kurtosis(index2)) >= 2*sum(Kurtosis(3:end))/2
        IMF_compress(:, CH) = imf(:, index1) + imf(:, index2);
    else
        IMF_compress(:, CH) = imf(:, index1) + imf(:, index2) + imf(:, index3);    
        IMF_compress(:, CH) = imf(:, index1) + imf(:, index2);
    end
%     
%     if CH>5
%         for i = 1:size(IMF_com_avg, 1)
%             IMF_com_avg(i, CH) = mean(IMF_compress(i, CH-5:CH));
%         end
%     else
%         IMF_com_avg(:, CH) = IMF_compress(:, CH);
%     end
% %     for i=1:size(IMF_com_avg, 1)
% %         if IMF_com_avg(i) <= 0
% %             IMF_com_avg(i) = 0;
% %         end
% %     end
% %     
%     [~,~,~,~,imfinse(:, CH)] = hht(IMF_com_avg(:, CH));
    
end
IMF_com_avg = mean(IMF_compress, 2);

% set signal<=0 to 0
for i=1:size(IMF_com_avg, 1)
   if IMF_com_avg(i) <= 0
       IMF_com_avg(i) = 0;
   end   
end

% find instantaneous energy for IMF_compress
[~,~,~,~,imfinse] = hht(IMF_com_avg);
a=1;
