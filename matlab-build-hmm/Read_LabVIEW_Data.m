% This code reads raw signals saved by the ReadRawSignals VI. The EMG and
% EEG data are saved as a 2-column bin file, which is read here and
% plotted.  This code also prompts for the accompanying EDF file, which you
% can plot to validate signals.

%%Binary File 
close all; clc; clear all;

%% Load Binary file
[filename,pathname]=uigetfile('*.bin','Select BIN file containing IR data:');
filename_data=[pathname,filename];
fid=fopen(filename_data,'r','ieee-be');


%Read timestamp data
data = fread(fid,[2 inf],'double')';
EEG = data(:,1); EMG = data(:,2);
% EEG(find(EEG<1E-20),:)=[]; EMG(find(EMG<1E-20),:)=[];

%% Load EDF
% Load signals
[filename,pathname]=uigetfile('*.edf','Select EDF file:');
fn=[pathname,filename];
HDR = sopen_noann(fn,'r',[2:3]);    % Read signal header 
nmin = HDR.NRec*HDR.Dur/60;         % Recording duration in minutes
S = sread_nothresh(HDR);   % Signals: 2 EEG, EMG, Piezo
lim=[0 3E4];
figure(2); subplot(411);
plot(EEG);
title('EEG (LabVIEW)');
xlim(lim);
subplot(412);
plot(EMG);
title('EMG (LabVIEW)');
xlim(lim);
subplot(413);
plot(S(:,1));
title('EEG (EDF)');
xlim(lim);
subplot(414);
plot(S(:,2));
title('EMG (EDF)');
xlim(lim);

fclose(fid);