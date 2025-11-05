% HMMBMWLCH.M (c) 2009 by Sridhar Sunderam
% Implementation of the Baum-Welch EM algorithm for estimation of HMM
% parameters from an observation set.
% Follows exposition by Andrew Fraser in "Hidden Markov models and
% dynamical systems", SIAM 2009.
%function [Hest,Sest,Htrue,Strue,y] = hmmbmwlch(y,H0)
function [Hest,Sest] = hmmbmwlch(y,H0)

switch nargin
    case 0       
        % Simulated observations
        [y,S,H0] = hmmsim;
    case 1
        % Get an initial guess for the model
        [dummy1,dummy2,H0] = hmmsim;
    otherwise
        if ~isstruct(H0)
            % Only the number of states is given; generate a toy model
            [dummy1,dummy2,Htrue] = hmmsim(H0);
            if isempty(y)
                y = dummy1; Strue = dummy2;
            end
            % Initial guess for model
            [dummy1,dummy2,H0] = hmmsim(H0);
        end
end

[Ly,alpha,gamma,beta] = hmmfb(y,H0);
H1 = reestimate(y,alpha,beta,gamma,H0);

n = 1;
while max(abs(H1.P0-H0.P0)) > 1e-8
    H0 = H1;
    [Ly,alpha,gamma,beta] = hmmfb(y,H0);
    H1 = reestimate(y,alpha,beta,gamma,H0);
    disp(H1.P0(:)');
    n = n+1;
end

Hest = H1;
Sest = hmmviterbi2(y,Hest);

return

function Hnew = reestimate(y,alpha,beta,gamma,Hold)
mu = Hold.mu;
Sigma = Hold.Sigma;
P0 = Hold.P0;
TR = Hold.TR;

nS = size(P0,1);    % no. of states
[T,nvar] = size(y);     % no. of time steps, variables

% Precompute values of P(y|s) given the form of the pdf
Py = zeros(nS,T);  % preallocate
for i=1:nS
    Py(i,:) = mvnpdf(y,mu{i},Sigma{i});
end

w = alpha.*beta;
%P0 = w(:,1);    % state priors
P0 = sum(w,2)/sum(w(:));    % state priors

W = zeros(nS,nS,T-1);
for t=1:T-1
    for s = 1:nS
        W(:,s,t) = alpha(:,t).*TR(:,s)*Py(s,t+1)*beta(s,t+1)/gamma(t+1);
%        w(:,s,t) = v(:,t) + log(TR(:,s)) + log(Py(:,t+1));
    end
end

for s=1:nS
%    mu{s} = w(s,:)*y;
    mu{s} = w(s,:)*y/sum(w(s,:));
    ydev = y-mu{s}(ones(T,1),:);
%    Sigma{s} = w(s*ones(nvar,1),:).*ydev'*ydev;
    Sigma{s} = w(s*ones(nvar,1),:).*ydev'*ydev/sum(w(s,:));
    Ws = sum(W(s,:,:),3);
    TR(s,:) = Ws/sum(Ws);
end

Hnew.mu = mu;
Hnew.Sigma = Sigma;
Hnew.P0 = P0;
Hnew.TR = TR;
return
