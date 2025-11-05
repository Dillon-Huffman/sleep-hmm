function [Class,Feat,class,feat] = ExtractClassifications(data_file,start_time,hours_to_read)
% This function is used to load data from LabVIEW saved files during 
% REM Sleep Restriction trials that were conducted during the STTR Phase II
% from 2016-2018. (March 2020, DH)

% Classifications are reduced from 1-second resolution to 4-second
% resolution for comparison to manual scores. This is done by using the
% timestamp to determine which classification correspond to 4-second
% epochs, and taking the mode of the classification for that time

%Inputs:
%   file: path\file to load (i.e. '\\directory\subdirectory\file.csv'
%   DayID: numeric input to identify baseline (1) or stimulation (2)
%   start_time: 1x3 vector containing time to start reading ([HH MM SS])
%   hours_to_read: number of hours to read, starting from start_time
%Outputs:
%   Class: subset of classifications for desired period
%   Feat: subset of classification features for desired period

end_time = datevec(datenum(start_time)+hours_to_read/24);

data = dlmread(data_file,',');

tvec = data(:,[9 7 8 10 11 12 13]);
secs = sum(tvec(:,end-1:end),2);
dt = diff(secs);

[~,ind_start]=intersect(tvec(:,1:6),start_time,'rows'); 
ind_start=ind_start+1; %Classifications are stored for the previous second
[~,ind_end]=intersect(tvec(ind_start:end,1:6),end_time,'rows');
ind_end = ind_end+1;
ind_end = ind_end + ind_start-1;

tvec2 = tvec(ind_start:ind_end,:);
dt = sum(tvec2(:,end-1:end),2);
dtr = [1;round(diff(rem(dt,4)))];

x=find(diff(dtr)>=2)+1;
x = [1; x];
dx = diff(x);

if sum(rem(dx(find(dx~=4)),4)==0)~= 0; %If dist betw epochs are div by 4..
    ix = find(dx~=4);
    ix2 = find((rem(dx(ix),4)~=0)==0);
    ix = ix(ix2);
    dx2 = dx(dx~=4)/4;
    dx2 = dx2(ix2);
    for ii = 1:length(ix)
        x_add = x(ix(ii))+flipud(cumsum(ones(dx2(ii)-1,1))*4);
        x = [x(1:ix(ii)); x_add; x(ix(ii)+1:end)];
        if ii+1 <= length(ix)
            ix(ii+1) = ix(ii+1)+length(x_add);
        end
    end
end

xe = x(find(diff(x)~=4)); %1-min widows to exclude
 
class = data(ind_start:ind_end,4);
feat = data(ind_start:ind_end,1:3);
for epoch = 1:length(x)-1
     Class(epoch,1) = mode(class(x(epoch):x(epoch+1)-1));
     Feat(epoch,:) = mean(feat(x(epoch):x(epoch+1)-1,:),1);
end
% 
% % nEx = length(xe); %number of excluded epochs
% class = Class;
% feat = Feat;
return
