% This function comes from MRC Cognition & Brain Sciences Unit (University of Cambridge)
%     http://imaging.mrc-cbu.cam.ac.uk/imaging/PhysNoise 
%     http://imaging.mrc-cbu.cam.ac.uk/imaging/PhysNoise?action=AttachFile&do=get&target=phys_noise.m

function [Y, T, data] = phys_noise_JL(varargin)
%PHYS_NOISE Summary of this function goes here
%   Detailed explanation goes here

ppg_file   = varargin{1};
info_file  = varargin{2};
dic_folder = varargin{3};
TR         = varargin{4};
Hz         = varargin{5};

%% start
P.TR        = TR; % Acquisition time (for working around dummy scans)
P.EndClip   = 1000; % Always seems to take a second and change to finish
P.StartClip = 0;    % Some systems have to have a pause between record onset and gradient onset
P.Hz        = Hz;   % Sampling Frequency
P.StartJump = 0;    % Number of shift samples to adjust at the beginning of the series
P.FirstScan = 1;
P.TickTime  = 1000/P.Hz;
[Y,T,P, data] = read_ppg(ppg_file, info_file, P);

% get the start and end of scanning scanning
% start - first series folder and first dicom
% end   - last series fold and last dicom

% % get series folders
%ser_folder = dic_folder(1:strfind(dic_folder, '/Series'));
% dce_f = dir(ser_folder);
% dc1_folder = fullfile(ser_folder,dce_f(3).name);
% dce_folder = fullfile(ser_folder,dce_f(end).name);
% 
% dc1_files = fullfl(dc1_folder,'1.3','dcm');
% dce_files = fullfl(dce_folder,'1.3','dcm');
% dc1_file  = dc1_files{1};
% dce_file  = dce_files{end};
% dc1_hdr   = spm_dicom_headers(dc1_file, false);
% dce_hdr   = spm_dicom_headers(dce_file, false);
% mri_start = dc1_hdr{1}.AcquisitionTime;
% mri_stop  = dce_hdr{1}.ContentTime;
% mri_dur   = mri_stop - mri_start;
% disp(['--- Total MRI time ---']);
% disp(['MRI start        : ' secs2hms(mri_start)]);
% disp(['MRI stop         : ' secs2hms(mri_stop)]);
% disp(['MRI duration     : ' secs2hms(mri_dur)]);
% 
% % load dicom headers
% dicoms = fullfl(dic_folder,'1.3','dcm');
% dicoms = dicoms(P.FirstScan:end);
% d_1st = 1;
% d_end = size(dicoms,1);
% hdr_1 = spm_dicom_headers(dicoms{d_1st}, false);
% hdr_e = spm_dicom_headers(dicoms{d_end}, false);
% 
% % dicom start (in seconds since midnight)
% d_start = hdr_1{1}.AcquisitionTime;
% d_stop  = hdr_e{1}.ContentTime;
% 
% % ppg start
% p_start = P.MRIStartTime/1000; % secs
% p_stop  = P.MRIStopTime/1000;  % secs
% p_dur   = p_stop - p_start;    % secs
% offset  = d_start - p_start;   % secs
% i_s = ceil(offset*P.Hz);       % secs * Hz = indices
% i_e = ceil((offset+(d_end*P.TR/1000))*P.Hz);
% %i_e = ceil((d_stop - p_start)*P.Hz);
% 
% disp(['Phys log start   : ' secs2hms(p_start)]);
% disp(['Phys log stop    : ' secs2hms(p_stop)]);
% disp(['Phys log duration: ' secs2hms(p_dur)]);
% disp(['Freq estimate    : ' num2str(size(Y,2)/p_dur) 'hz, Freq used: ' num2str(P.Hz) 'hz'])
% disp(['--- Run information ---']);
% disp(['Number of dicoms : ' num2str(d_end)]);
% disp(['First dicom      : ' secs2hms(d_start)]);
% disp(['Last dicom       : ' secs2hms(d_stop)]);
% disp(['Phys-dicom       : ' secs2hms(d_start-p_start)]);
% disp(['Phys indices     : ' num2str(i_s) ' ... ' num2str(i_e)]);
% 
% % extract ppg data pertaining to this run
% Y = Y(i_s:i_e);
% T = T(i_s:i_e);
% T = T - T(1);
% disp(['Run duration from dicoms      : ' num2str((size(Y,2).*P.TickTime)/1000/60) ' mins']);
% disp(['Phys. data for this run       : ' num2str((T(end)-T(1))/60) ' mins']);
% 
% % check if there is any data in Y
% if numel(unique(Y))==1
%     disp('-- ERROR --');
%     disp('The physiological signal extracted for this time period is completely flat!');
%     disp(['Check the file ' ppg_file]);
%     disp('for the time period reported above');
%     disp('Noise regressor cannot be estimated.');
%     disp('This might be because your recording equipment was detached during this scanning period');
%     error('Exiting ..');    
% 
% % saving
%[fpath fname fext] = fileparts(ppg_file);
% save(fullfile(fpath,[fname '.mat']),'Y','T','P');
end

function [signal,time,P, data] = read_ppg(pulsefile,info_file,P)

fclose('all');

fid=fopen(pulsefile);
% edit by James Lynch -10/18/17
ignore=textscan(fid,'%s %s %s',5);  % Ignore first 5 rows of header information (3 columns) 
ignore2=textscan(fid, '%s %s %s %s', 1);  % Ignore 6th row - column names (4 columns)

data   = textscan(fid,'%d %s %d %s'); %Read data until end of u16 data.
footer = textscan(fid,'%s');   %Read in remaining data (time stamps and statistics)

fid2=fopen(info_file);
ignore_info=textscan(fid2, '%s %s %s', 6);  % Ignore first 6 rows of header information (3 columns)
ignore_info2=textscan(fid2, '%s %s %s %s', 1);  % Ignore 7th row - column names (4 columns)

data_info   = textscan(fid2,'%d %d %d %d'); %Read data until end of u16 data.
footer_info = textscan(fid2,'%s');   %Read in remaining data (time stamps and statistics)

% ----- Added by James Lynch - 10/18/17 -----
% Our Siemens log files do not have footer information
% To get the start and end times of the log file and functional run, use the 1st and last time stamps from 'data'
firstTic_phys = data{1,1}(1,1);  
firstTic_info = data_info{1,3}(1,1);

lastTic_phys = data{1,1}(end,1);
lastTic_info = data_info{1,4}(end,1);

% Each tic is 2.5 ms, so multiply the tics by 2.5 to convert to ms
P.LogStartTime = firstTic_phys * 2.5;
P.LogStopTime = lastTic_phys * 2.5;
P.MRIStartTime = firstTic_info * 2.5;
P.MRIStopTime = lastTic_info * 2.5;


% if ~size(footer{1},1)
%     disp('No footer found in the log file!');
%     disp('Hence, the start of the logging is unknown');
%     disp('and match with dicoms cannot be calculated.');
%     disp('This might indicate that the equipment was switched off unexpectedly.');
%     error('Exiting..');
% end

% %Get time stamps from footer:
% for n=1:size(footer{1},1)
%     if strcmp(footer{1}(n),'LogStartMDHTime:')  %log start time
%         P.LogStartTime=str2num(footer{1}{n+1});
%     end
%     if strcmp(footer{1}(n),'LogStopMDHTime:')   %log stop time
%         P.LogStopTime=str2num(footer{1}{n+1});
%     end
%     if strcmp(footer{1}(n),'LogStartMPCUTime:') %scan start time
%         P.MRIStartTime=str2num(footer{1}{n+1});
%     end
%     if strcmp(footer{1}(n),'LogStopMPCUTime:')  %scan stop time
%         P.MRIStopTime=str2num(footer{1}{n+1});
%     end
% end

% Remove the systems own evaluation of triggers.
t_on  = find(data{1} == 5000);  % System uses identifier 5000 as trigger ON
t_off = find(data{1} == 5003);  % System uses identifier 5003 as trigger OFF

% Filter the trigger markers from the data
Hz = P.Hz;
data_t=transpose(1:length(data{1}));
indx = setdiff(data_t(:),union(t_on,t_off)); %Note: depending on when the scan ends, the last size(t_off)~=size(t_on).
signal = data{1}(indx);
time = (1:length(signal))./Hz;

% Clip the time series to begin and end with the scan.
if P.TR < 1000
  start = P.StartClip*(Hz/1000);
  stop  = length(signal) - Hz*floor([length(signal) - EndClip*(Hz/1000)]/Hz);
else
  start = P.StartClip*(Hz/1000)+Hz/2;
  stop  = P.EndClip*(Hz/1000)+mod(length(signal(1:end)),Hz)+ Hz/2;
  %stop  = length(signal) - Hz*floor([length(signal) - P.EndClip*(Hz/1000)]/Hz) + Hz/2;
end;

% Reset vectors;
signal = double(signal(start+1-P.StartJump:end-stop-P.StartJump)');
time = time(start+1-P.StartJump:end-stop-P.StartJump);
time = time(:)-time(1);
end

function [Y,T,P] = read_pulsef(pulsefile)

[fid, emsg]  = fopen(pulsefile);

if fid==-1
    error(emsg);
end

tline = fgetl(fid);
C = [];

while ischar(tline)
    C{end+1,1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

Y = str2num(C{1,1});
p = Y(1:5); % first 5 values are parameters
Y = Y(6:end); % resize Y, starting with the 6th value of the first line are the voltage values themselves that reflect the pulse ox.

P.TickTime   = p(3); %time (in msec) between subsequent measurements
P.StartMDH   = get_timeval(C{11,1});
P.StopMDH    = get_timeval(C{12,1});
P.StartMPCU  = get_timeval(C{13,1});
P.StopMPCU   = get_timeval(C{14,1});
end

function timeval = get_timeval(ccell)
[~, timeval] = strtok(ccell);
timeval = str2num(rmwhite(timeval));
end
