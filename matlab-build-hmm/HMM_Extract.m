function HMM_extract()
% This script opens an EDF file, extracts EEG/EMG/piezo features, fits a
% HMM to the data, writes the model parameters to a CSV, and saves
% estimated sleep scores, proportions, and bout information to mat files

%% Global settings, filters
fs = 400;   % data sampling rate
nsec = 1;   % epoch size
N = nsec*fs;    % samples/epoch
[NBf,Hf_EEG,Hf_EMG,Hf_Piezo,Hf_EMG2] = makefilters(fs);

%% Load EEG, EMG, Piezo signals and sleep scores
% Load signals
fn = ['NSL16_0519_2'];
HDR = sopen_noann([fn '.edf'],'r',[1:4]);    % Read signal header 
nmin = HDR.NRec*HDR.Dur/60;         % Recording duration in minutes

k = 1;           % Segment counter
seglength = 60;  % 60 min segments

%% Compute signal features for each segment
while k <= 2%nmin/seglength    
    disp(k); 

    S = sread_nothresh(HDR,seglength*60,(k-1)*seglength*60);   % Signals: 2 EEG, EMG, Piezo
    %Fill gaps
    sf = find(isnan(S));
    if ~isempty(sf), S = inpaint_nans(S,4); end
    %Remove DC trend
    S = detrend(S);
    %Remove 60Hz noise
    S = filter(NBf,S);
    ts = ([(k-1)*size(S,1)+1:k*size(S,1)]'-1)/fs; % Time Vector
    Ts{k} = ts;

    %% Computations for entire segment (one feature vector extracted per epoch)
    
    % EEG band power estimates by filtering
    for i=1:length(Hf_EEG)
        sigf = filter(Hf_EEG{i},S(:,1:2));
        sigf = filter(ones(N,1),N,sigf.^2);
        BP{i} = sigf(N:N:end,:);
    end    
    R{k} = cat(3,BP{:});
    EEGf = filter(Hf_EEG{end},S(:,1:2));   % 80Hz lowpass-filtered EEG (for  TE, LL calculations)
    
    % Teager energy and line length of EEG 
    te = filter(ones(N,1),N,EEGf(2:end-1,:).^2-EEGf(1:end-2,:).*EEGf(3:end,:)); 
    te = [te(1,:);te;te(end,:)]; TE_EEG{k} = te(N:N:end,:);
    ll = filter(ones(N,1),N,abs(diff(EEGf)));
    ll = [ll;ll(end,:)];LL_EEG{k} = ll(N:N:end,:);   % Line length
        
    % Estimated EMG power:
    % Original
    EMG{k} = filter(ones(N,1),N,S(:,3).^2); EMG{k} = EMG{k}(N:N:end);
    % Filtered
    EMGf{k} = filter(Hf_EMG,S(:,3)); EMGf{k} = filter(ones(N,1),N,EMGf{k}.^2); EMGf{k} = EMGf{k}(N:N:end);
    EMGf2{k} = filter(Hf_EMG2,S(:,3)); EMGf2{k} = filter(ones(N,1),N,EMGf2{k}.^2); EMGf2{k} = EMGf2{k}(N:N:end);

    % Piezo signal
    Pzof = filter(Hf_Piezo{end},S(:,4)); % 0.5-20 Hz for TE computation
    Pzol = filter(Hf_Piezo{1},S(:,4));   % 0.5-5 Hz for breath rate, regularity computation
   
    % Piezo power:
     % Original (raw signal)
    Pz{k} = filter(ones(N,1),N,S(:,4).^2);
    Pz{k} = Pz{k}(N:N:end);
     % Filtered (0.5-20 Hz for TE, AE computation)
    Pzf{k}= filter(ones(N,1),N,Pzof.^2);
    Pzf{k} = Pzf{k}(N:N:end);
 
    % Teager energy of piezo signal ((0.5-20Hz)
     te = filter(ones(N,1),N,Pzof(2:end-1).^2-Pzof(1:end-2).*Pzof(3:end)); 
     te = [te(1);te;te(end)]; TE_Piezo{k}(:,2) = te(N:N:end);
     ll = filter(ones(N,1),N,abs(diff(Pzof)));
     ll = [ll;ll(end)];LL_Piezo{k} = ll(N:N:end,:);   % Line length
        
     
     % Features extracted from the Hilbert transform
     [BR{k},REG{k},ENV{k}] = hilbert_features(-Pzol,fs,N,ts);
         
     % Features extracted from the Fourier transform
     [BR{k}(:,3),REG{k}(:,3)] = fourier_features(Pzol,fs,N);    
     
    k = k+1;
    
end
EEG=cat(1,R{:}); EMG=cat(1,EMGf{:}); Pz = cat(1,Pz{:}); EMG2=cat(1,EMGf2{:});
for i=1:7; eeg(:,i)=EEG(:,2,i); end; EEG=eeg; 
EMG_norm=EMG./EMG2; % (80-100Hz)/(10-20 Hz)


[Scores, sleep_metrics, HMM, delta, Feat] = GOHMM(EEG,EMG_norm); %Change EMG to Pz for Piezo feat


fn=['VAR16_0913'];
% [Scores2, sleep_metrics2, HMM2, delta2, Feat2] = GOHMM(EEG,EMG_norm); %Change EMG to Pz for Piezo feat
% figure; plot([Feat(:,3) Feat2(:,3)]);
save(['HMM_' fn '.mat'],'HMM');
Model = [cat(1,HMM.mu{:}); cat(1,HMM.Sigma{:}); HMM.P0'; HMM.TR; HMM.norm];
dlmwrite(['Model_' fn '.csv'],Model, 'delimiter', ',', 'precision', 16);
save([fn '_sleepMetrics.mat'],'Scores','sleep_metrics');
% save('Features_MATLAB.mat','Feat','delta','Scores','HMM');
% save([fn '_feat.mat'],'HDR','R','TE_EEG','LL_EEG','EMG','EMGf','TE_Piezo','LL_Piezo','Pz','Pzf','BR','REG','ENV');
return

function [NBf,Hf_EEG,Hf_EMG,Hf_Piezo,Hf_EMG2] = makefilters(fs)
%Filter parameters

% Narrowband 60Hz bandstop filter as away of removing line noise
[nn,wn] = buttord([55 65]/(fs/2), [59 61]/(fs/2),3,5);
[f2,f1] = butter(nn,wn,'stop');
NBf = dfilt.df1(f2,f1);   
isstable(NBf) 
NBf.PersistentMemory = 1; 
% Delta (low, high), theta, alpha, sigma, beta, gamma (and total) limits
% for band power computation from EEG/EMG
flim_EEG_p = [0.5 2 6 9 13 30 0.5; 2 4 9 13 30 80 80]';
flim_EEG_s = [0.1 1 4 8 11 25 0.1; 3 5 10 15 32 85 90]';

for i=1:size(flim_EEG_p,1)
     [nn,wn] = buttord(flim_EEG_p(i,:)/(fs/2), flim_EEG_s(i,:)/(fs/2),3,5);
     ord(i)=nn;
     [b{i},a{i}] = butter(nn,wn);
     Hf_EEG{i} = dfilt.df1(b{i},a{i});  
     isstable(Hf_EEG{i})
     Hf_EEG{i}.PersistentMemory = 1;
end

% EMG filters
flim_EMG_p = [80 100];
flim_EMG_s = [75 110];
flim_EMG_p2 = [10 20];
flim_EMG_s2 = [8 22];
[nn,wn] = buttord(flim_EMG_p/(fs/2), flim_EMG_s/(fs/2),3,5);
[nn2,wn2] = buttord(flim_EMG_p2/(fs/2), flim_EMG_s2/(fs/2),3,5);
[d,c] = butter(nn,wn);
[d2,c2] = butter(nn2,wn2);
Hf_EMG = dfilt.df1(d,c); 
Hf_EMG2 = dfilt.df1(d2,c2); 
isstable(Hf_EMG)
isstable(Hf_EMG2)
Hf_EMG.PersistentMemory = 1;
Hf_EMG2.PersistentMemory = 1;
% Piezo filters
flim_Piezo_p = [0.5 0.5; 5 20]';
flim_Piezo_s = [0.1 0.1; 6 22]';
for i=1:size(flim_Piezo_p,1)
    [nn,wn] = buttord(flim_Piezo_p(i,:)/(fs/2), flim_Piezo_s(i,:)/(fs/2),3,5);
     [f{i},e{i}] = butter(nn,wn);
     Hf_Piezo{i} = dfilt.df1(f{i},e{i});  
     isstable(Hf_Piezo{i})
     Hf_Piezo{i}.PersistentMemory = 1;
end
return

% Features extracted from the Hilbert transform
function [br,reg,env] = hilbert_features(x,fs,N,ts)
hx = hilbert(x);%Piezo

% Coeff of variation in signal envelope (to estimate CV of breath amplitude)
env = abs(hx); envm = filter(ones(N,1),N,env); 
envstd = filter(ones(N,1),N,sqrt((env-envm).^2));
env = envstd(N:N:end)./envm(N:N:end);    

% Breath rate estimated from number of rotations
ph = unwrap(angle(hx)+pi);  % unwrapped Hilbert phase
breath_ind = [0;diff(round(ph/(2*pi)))>0];
ind1 = find(breath_ind); tb = ts(ind1); % breath times
gw = gausswin(N);
br1 = filter(gw,sum(gw)/N,breath_ind)/(N/fs);
%figure; plot([Pzol*10 breath_ind(:) angle(hx)+pi br1]);
br1 = br1(N:N:end);

% Breath regularity, estimated as mean phase coherence with respect to
% a reference signal (same size epoch, 1s back)
PH = [[ph(1:fs); ph] [ph; ph(end-fs+1:end)]];
reg1 = abs(filter(ones(N,1),N,exp(1j*(PH(:,1)-PH(:,2))))); reg1(1:fs,:) = [];
reg1 = reg1(N:N:end);
%br1 = [ph(1);ph(N:N:end)]/(2*pi); BR1{k} = diff(br1)/(N/fs);
% Another breath rate and Rayleigh index of regularity
nb = cumsum(breath_ind); nb = diff([0;nb(N:N:end)]);  % no. of breaths in each epoch
dt = [0;diff(tb)];  % interbreath intervals
ibi = breath_ind; ibi(ind1) = dt;    % interbreath interval at each breath
ibi = filter(ones(N,1),1,ibi); ibi = ibi(N:N:end); % sum of interbreath intervals in each epoch
br2 = nb./ibi;   % breath rate as inverse of mean interbreath interval
brext = br2(:,ones(N,1))'; brext = brext(:);
phasor = exp(1j*2*pi*brext.*ts); phasor = phasor.*breath_ind;
rr = filter(ones(N,1),1,phasor); rr = rr(N:N:end)./nb;
reg2 = rr.*conj(rr);

br = [br1 br2];
reg = [reg1 reg2];
return

% Features extracted from the Fourier transform
function [br,reg] = fourier_features(x,fs,N)
[s,f,t] = spectrogram(detrend(x),N,0,[],fs);
%[ss,ff,tt] = spectrogram(detrend(rem(ph,2*pi)),N,0,[],fs); % or from Hilbert phase signal
sum_spec = sum(abs(s)); sum_spec = sum_spec(ones(size(s,1),1),:);
[reg,fmax] = max(abs(s)./sum_spec); reg = reg(:); 
br = f(fmax);
%REG2{k} = exp(mean(log(abs(diff(abs(ss))))))...
%    ./mean(abs(diff(abs(ss)))); REG2{k} = REG2{k}(:);  % spectral flatness (Weiner entropy)
return
