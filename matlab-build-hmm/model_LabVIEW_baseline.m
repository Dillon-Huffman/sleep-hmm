% This function fits an HMM to features previously recorded by LabVIEW and
% returns the model and unsupervised classifications.  The input parameter 
% is data saved to a csv by LabVIEW, in which the first three columns are 
% the delta theta ratio, high-low ratio, and EMG feature respectively.
% (Dillon Huffman, 11/10/2016)

function [S,HMM,Feat] = model_LabVIEW_baseline(Data);
%Model original features into clusters
Feat = Data(:,1:3);
[HMM,S] = fitHMM(Feat,3);
%Extract mean EMG power during NREM to serve as Normalization factor
norm=[0 0 HMM.mu{2}(3)]; 
%Normalize EMG feature
Feat(:,3)=Feat(:,3)-norm(3);
%Rebuild model on new feature set
[HMM,S] = fitHMM(Feat,3);
%Assign normalization factor to HMM to use in real time modeling
HMM.norm=norm;
%Pull out Q values during NREM
Q_NREM = Feat(find(S==2),1);
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
        TR(i,q) = sum([Inf;score]==i & [score;Inf]==q)/sum([Inf;score]==i);
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