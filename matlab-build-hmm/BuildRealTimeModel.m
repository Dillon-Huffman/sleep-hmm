%This script prompts the user to select a csv file saved by LabVIEW during
%a baseline, which is used to build a model for real-time classification.

%Select .CSV file to open
[file_o,path_o]=uigetfile('*.csv','Select CSV containing baseline features:'); 
%Read .CSV file
data=csvread([path_o,file_o]);
%Select 8 hours of data, ignoring the first 20 sec
data = data(20:7*3600+20,:);

%Model features
[S,HMM,Feat]=model_LabVIEW_baseline(data);
figure; plot([Feat S'*3-15]);

% Write model parameters to CSV file
Model = [cat(1,HMM.mu{:}); cat(1,HMM.Sigma{:}); HMM.P0'; HMM.TR; HMM.norm];
[file_s,path_s] = uiputfile('*.csv','Save Model to File:','Model_.csv');
dlmwrite([path_s,file_s],Model, 'delimiter', ',', 'precision', 16);

