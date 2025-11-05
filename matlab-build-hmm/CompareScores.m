Scores_VI=csvread('Features_d2.csv');
Scores_VI=Scores_VI(:,11);
load Scores.mat; Scores=Scores';

% Scores=Scores(:,ones(4,1))';Scores=Scores(:);
% Scores=Scores(1:length(Scores_VI));

w1=Scores==1; w2=Scores_VI==1;
cm=confusionmat(w1,w2);
spec(1,1) = cm(1,1)/(cm(1,1)+ cm(1,2));
sens(1,1) = cm(2,2)/(cm(2,1)+ cm(2,2));

n1=Scores==2; n2=Scores_VI==2;
cm=confusionmat(n1,n2);
spec(1,2) = cm(1,1)/(cm(1,1)+ cm(1,2));
sens(1,2) = cm(2,2)/(cm(2,1)+ cm(2,2));

r1=Scores==3; r2=Scores_VI==3;
cm=confusionmat(r1,r2);
spec(1,3) = cm(1,1)/(cm(1,1)+ cm(1,2));
sens(1,3) = cm(2,2)/(cm(2,1)+ cm(2,2));

figure; stairs([Scores Scores_VI+4]);
ylim([-1 9]);

cm=confusionmat(Scores,Scores_VI)