function [emf] = cosumek(emf_thin, emf_nom, emf_adj, k, pickedCH, alignCH, testConfig, sensorNom)
% this function is to consume k to modify emf_thin to proper numbers
% emf_thin is FM thin without abcd
% emf_nom is FM nom without abcd
% emf_adj is FM nom with abcd

% align at alignCH in ms
emf_align = emf_nom;
alignT = testConfig.chstarttime(alignCH);
for i = alignT:testConfig.chendtime(testConfig.adec_idx_end(sensorNom))
    emf_align(i) = emf_thin(i) / (emf_thin(alignT) / emf_nom(alignT));   
end

% apply k at pickedCH in ms
pickedT = testConfig.chstarttime(pickedCH);
%D = (emf_align(pickedT) - emf_nom(pickedT))/k;

emf = emf_adj;
for j = alignT:testConfig.chendtime(testConfig.adec_idx_end(sensorNom))
    emf(j) = emf_adj(j) - (abs(emf_align(j) - emf_nom(j))/k);   
    a(j) = ((emf_align(j) - emf_nom(j)));
    if emf(j) <= 0
        emf(j) = 0.01;
    end       
end

end