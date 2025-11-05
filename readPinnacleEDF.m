function readPinnacleEDF
% opens EDF files epxorted using Sirenia Sleep v1.5.0 (missing annotations)
%
% Dillon Huffman 2025
s
    [f,p] = uigetfile('*.edf');
    if isequal(p,0) && isequal(f,0)
        error('File selection aborted by user');
    end
    file = fullfile(p,f);

    % Read EDF header
    [hdr,hdrLength] = readEDFHeader(file);
    
    nChan = hdr.channels-1; %number of signals in file (excluding annotations)
    chanRead = [1 2 3 4 5];  %channels to read

    fid = fopen(file,'r','ieee-le');
    fseek(fid,hdrLength,"bof");
    
    tmp = zeros(hdr.duration*hdr.samplerate(1),numel(chanRead)); %assumes all channels have same sample rate
    wb = waitbar(0,'Reading EDF File.. ')
    for i = 1:hdr.records
        waitbar(i/hdr.records,wb);
        k=1;
        for j = 1:nChan 
            if ismember(j,chanRead)
                % data in file is interleved by channel for each data
                % record, written as int16 (2 bytes)
                % nBytes per data record is hdr.duration * header.sampleRate * 2 bytes (int16)
                tmp(:,k) = fread(fid,hdr.duration*hdr.samplerate(j),'int16'); 
                k=k+1;
            else
                fseek(fid,hdr.duration*hdr.samplerate(j)*2,'cof'); %skip data for channels not read                        fseek(fid,hdr.samplerate(end)*2,'cof');
            end
        end
        fseek(fid,hdr.samplerate(end)*2,'cof'); %skip annotations
        sig{i,1} = tmp*((hdr.physmax(j)-hdr.physmin(j))/(hdr.digmax(j) - hdr.digmin(j)));
    end
    fclose(fid);
end
function [header,endPos] = readEDFHeader(filename)
    % code extracted from ReadEDF function uploaded to MathWorks exchange by
    % Shapkin Andrey. Modified to only read header fields.
    
    
    % filename - File name
    % data - Contains a signals in structure of cells
    % header  - Contains header
     
    fid = fopen(filename,'r','ieee-le');
    
    %%% HEADER LOAD
    % PART1: (GENERAL) 
    hdr = char(fread(fid,256,'uchar')'); 
    header.ver=str2num(hdr(1:8));               % 8 ascii : version of this data format (0)
    header.patientID  = char(hdr(9:88));        % 80 ascii : local patient identification
    header.recordID  = char(hdr(89:168));       % 80 ascii : local recording identification
    header.startdate=char(hdr(169:176));        % 8 ascii : startdate of recording (dd.mm.yy)
    header.starttime  = char(hdr(177:184));     % 8 ascii : starttime of recording (hh.mm.ss)
    header.length = str2num (hdr(185:192));     % 8 ascii : number of bytes in header record
    reserved = hdr(193:236); % [EDF+C       ]   % 44 ascii : reserved
    header.records = str2num (hdr(237:244));    % 8 ascii : number of data records (-1 if unknown)
    header.duration = str2num (hdr(245:252));   % 8 ascii : duration of a data record, in seconds
    header.channels = str2num (hdr(253:256));   % 4 ascii : number of signals (ns) in data record
    
    header.labels=cellstr(char(fread(fid,[16,header.channels],'char')')); % ns * 16 ascii : ns * label (e.g. EEG FpzCz or Body temp)
    header.transducer =cellstr(char(fread(fid,[80,header.channels],'char')')); % ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode)
    header.units = cellstr(char(fread(fid,[8,header.channels],'char')')); % ns * 8 ascii : ns * physical dimension (e.g. uV or degreeC)
    header.physmin = str2num(char(fread(fid,[8,header.channels],'char')')); % ns * 8 ascii : ns * physical minimum (e.g. -500 or 34)
    header.physmax = str2num(char(fread(fid,[8,header.channels],'char')')); % ns * 8 ascii : ns * physical maximum (e.g. 500 or 40)
    header.digmin = str2num(char(fread(fid,[8,header.channels],'char')')); % ns * 8 ascii : ns * digital minimum (e.g. -2048)
    header.digmax = str2num(char(fread(fid,[8,header.channels],'char')')); % ns * 8 ascii : ns * digital maximum (e.g. 2047)
    header.prefilt =cellstr(char(fread(fid,[80,header.channels],'char')')); % ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
    header.samplerate = str2num(char(fread(fid,[8,header.channels],'char')')); % ns * 8 ascii : ns * nr of samples in each data record
    
    % Hot fix for incorrect sample rate reading (assuming Pinnacle had error in EDF
    % writing). Signal sample rates appear to be 10x higher than actual.
    % Assuming this is somehow related to the EDF record duration being 10s.
    header.samplerate(1:end-1) = header.samplerate(1:end-1)/header.duration;
    
    reserved = char(fread(fid,[32,header.channels],'char')'); % ns * 32 ascii : ns * reserved
    endPos = ftell(fid);
    fclose(fid);
end
function [scaleFactor,offset] = getScaling(physicalMax,physicalMin,...
    digitalMax,digitalMin,scalFactorFlag)
    %getScaling calculate required scaling and dc offset for the signals stored
    %in EDF/EDF+ files.
    %
    %   This function is for internal use only. It may change or be removed.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    if scalFactorFlag
        scaleFactor = (physicalMax - physicalMin) ./ (digitalMax - digitalMin);
        offset = physicalMin - scaleFactor .* digitalMin;
    else
        reqLength = length(physicalMax);
        scaleFactor = ones(reqLength,1);
        offset = zeros(reqLength,1);
    end
end