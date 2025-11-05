% This function fits a GOHMM (using GMM clusters as initial
% guess)to EEG/EMG features and estimates sleep scores and metrics. Farid
% Yaghouby 11/10/2015
function [Scores,sleep_metrics,HMM,delta,Feat] = GOHMM(EEG,EMG)

% nmin = length(EMGf)*nsec/60;
% t = [1:nmin/(nsec/60)]'*nsec;

R = EEG; EMGf = EMG;
f = 20; % moving average window length
nsec=1;
Feat = [sum(R(:,1:2),2)./R(:,3)  sum(R(:,4:6),2)./sum(R(:,1:3),2) EMGf];%feature set: delta/theta; hi/lo and EMG power
Feat = filter(ones(f,1),f,Feat);% 20s moving average smoothing
Feat = 10*log10(Feat);% features in log scale

%Unsupervised HMM
[HMM,S] = fitHMM(Feat,3);
norm=[0 0 HMM.mu{2}(3)];
Feat(:,3)=Feat(:,3)-HMM.mu{2}(3);
HMM.norm=norm;


% delta=0;
% co = cat(1,HMM.mu{:});
% S = matchstate(co,S);
%%  Get final scores and sleep metrics
Scores = S; 

% Estimate  sleep metrics for SegWay scores
[Prop,Bouts,Dur] = compute_sleep_metrics(S,nsec);
sleep_metrics = [Prop; Bouts; Dur];

end

%% Estimate sleep metrics from sleep scores
function [Prop,Bouts,Dur] = compute_sleep_metrics(scores,nsec)
for i=1:3
    h = (scores==i);
    Prop(i) = mean(h);
    Bouts(i) = sum(diff(h)==1);
    Dur(i) = sum(h)*nsec/Bouts(i);
end
end

function [HMM,score] = fitHMM(X,n)
 options = statset('Display','final','maxiter',5000,'TolFun',1e-15);% GMM options
 err = 0;counter=0;
        % Fit GMM to training data
        while err==0 && counter<10
        try
% Sleep-wake model            
%sw_GMM = fitGMM(SlFeat(:,end),2);
objsw = gmdistribution.fit(X(:,end),2,'replicates',15,'options',options);
sw = cluster(objsw,X(:,end));
sw = sw+3;
% %Assign states feature trends
C = objsw.mu;
Stmp = sw*NaN;
[dummy,ii] = sort(C); % sort by increasing EMG 
Stmp(sw==ii(1)+3) = 2;Stmp(sw==ii(end)+3) = 1;

S = zeros(size(Stmp,1),1);S(Stmp==1)=1;
%%Sleep models
objnr = gmdistribution.fit(X(Stmp==2,1),2,'replicates',15,'options',options);
nr = cluster(objnr,X(Stmp==2,1));
nr = nr+3;
 %Assign states feature trends
C = objnr.mu;
Stm = nr*NaN;
[dummy,ii] = sort(C); % sort by increasing delta/tetha 
Stm(nr==ii(1)+3) = 3;Stm(nr==ii(end)+3) = 2;
S(S==0)=Stm;
score = S(:);

% Fit HMM using preliminary scores to determine an initial guess
for i=1:n
    p0(i) = mean(score==i);         % Initial prior probabilities
    mu{i}= mean(X(score==i,:));     % Mean feature vector in each state
    sigma{i} = cov(X(score==i,:));  % Covariance of features in each state
    % Transition probabilities between states
    for q=1:n
        TR(i,q) = sum([score;Inf]==i & [Inf;score]==q)/sum([score;Inf]==i);
    end
end
%TR_GMM = TR;
 %Apply restrictions to transition matrix
TR(TR==0)=0.00001;
RM = ones(3,3);
RM(1,3)=0;
RM(3,2)=0;
TR= TR.*RM;
for i=1:3 TR(i,:)=TR(i,:)./sum(TR(i,:)); end
H0.mu = mu; H0.Sigma = sigma; H0.P0 = p0'; H0.TR = TR;
[HMM,score] = hmmbmwlch(X,H0);


            err = 1;
        catch
            err = 0;
        end
        counter = counter+1;
        end
end
       
% % function S = matchstate(C,S0)
% % %matching model state labels based on teh feature distribution. Farid
% % %Yaghouby 10/17/2013
% % %   Detailed explanation goes here
% % S = S0*NaN;
% % %ns = 3;
% % [dummy,ii2] = sort(C(:,end)); % sort by increasing EMG power ratio
% % S(S0==ii2(3)) = 1;
% % C(ii2(3),1)=Inf;
% % [dummy,ii1] = sort(C(:,1)); % sort by increasing delta/theta power ratio
% % S(S0==ii1(1)) = 3;
% % S(S0==ii1(2)) = 2;
% %    
% % end