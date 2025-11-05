% HMMVITERBI.M (c) 2009 by Sridhar Sunderam
% Viterbi algorithm for hidden Markov models
% Estimate most probable sequence of hidden states for a given sequence of
% observations and hidden Markov model H.
% Follows exposition by Andrew Fraser in "Hidden Markov models and
% dynamical systems", SIAM 2009.

%Lines 60-66 disabled to remove backtracking (DMH 11/10/2016)

function [S,delta] = hmmviterbi(y,H)
% load vittest.mat
% y=Feat_VI; H=HMM;
switch nargin
    case 0
        % HMM
        mu{1} = [1 -1]; Sigma{1} = [.9 .4; .4 .3];
        mu{2} = [0 0]; Sigma{2} = [1 0; 0 1];
        P0 = [0.25 0.75]';           % State priors
        TR = [0.3 0.7; 0.6 0.4];    % Transition probabilities
        
        % Simulated observations
        x1 = mvnrnd(mu{1}, Sigma{1}, 100);
        x2 = mvnrnd(mu{2}, Sigma{2}, 100);
        y = [x1;x2];
    case 1
        % HMM
        mu{1} = [1 -1]; Sigma{1} = [.9 .4; .4 .3];
        mu{2} = [0 0]; Sigma{2} = [1 0; 0 1];
        P0 = [0.25 0.75]';           % State priors
        TR = [0.3 0.7; 0.6 0.4];    % Transition probabilities
    otherwise
        % Extract variables from struct H
        fn = fieldnames(H);
        for i=1:length(fn)
            eval([fn{i} ' = getfield(H,fn{i});']);
        end
end

% Forward algorithm to compute the likelihood of an observation sequence
% given the model.
nS = size(P0,1);    % no. of states
T = size(y,1);     % no. of time steps

% Precompute values of P(y|s) given the form of the pdf
Py = zeros(nS,T);  % preallocate
for i=1:nS
    Py(i,:) = mvnpdf(y,mu{i},Sigma{i});
end

delta = zeros(nS,T);
psi = zeros(nS,T);
delta(:,1) = log(Py(:,1).*P0);
for t=2:T
    for s=1:nS
        [D,psi(s,t)] = max(delta(:,t-1)+ log(TR(:,s)));        
        delta(s,t) = D + log(Py(s,t));
    end
end

[dummy,S] = max(delta);

% This section was removed to more closely match the current 
% LabVIEW implementation, which does not backtrack (DMH 11/10/2016)
% {
% [dummy,S(T)] = max(delta(:,T));
% for t=T-1:-1:1
%     S(t) = psi(S(t+1),t+1);
% end
% S = S(:);
% }
%{
% The code below gives a staggered optimal state sequence: must correct.
u(:,1) = log(Py(:,1).*P0);
v(:,1) = u(:,1);
for t=1:T-1
    for s=1:nS
        w(:,s,t) = v(:,t) + log(TR(:,s)) + log(Py(:,t+1));
        [dummy,B(s,t+1)] = max(w(:,s,t));
        v(s,t+1) = w(B(s,t+1),s,t);
    end
end

[dummy,S(T)] = max(v(:,T));
for t=T-1:-1:1
    S(t) = B(S(t+1),t+1);
end
%}
return