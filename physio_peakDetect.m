% James Lynch 10/31/17
% This function combines aspects from the 'phys_noise' function and PhLEM toolbox/peakdet function
%   'phys_noise': MRC Cognition & Brain Sciences Unit (University of Cambridge)
%       http://imaging.mrc-cbu.cam.ac.uk/imaging/PhysNoise 
%       http://imaging.mrc-cbu.cam.ac.uk/imaging/PhysNoise?action=AttachFile&do=get&target=phys_noise.m
%   PhLEM: from Geoffrey Aguirre Lab (UPenn)
%       https://cfn.upenn.edu/aguirre/wiki/public:pulse-oximetry_during_fmri_scanning
%       https://sites.google.com/site/phlemtoolbox/
%   'peakdet': Marco Burges
%       https://www.mathworks.com/matlabcentral/fileexchange/47264-peakdet

function [pks,dep,pid,did,peak_times] = physio_peakDetect(varargin)

% Taken from 'phys_noise'
data_file   = varargin{1};
info_file  = varargin{2};
TR         = varargin{3};
Hz         = varargin{4};
th         = varargin{5};  % Should be 4000 for PULSE, 2300 for RESP (RESP threshold still needs some fine tuning)


P.TR        = TR; % Acquisition time (for working around dummy scans)
P.EndClip   = 1000; % Always seems to take a second and change to finish
P.StartClip = 0;    % Some systems have to have a pause between record onset and gradient onset
P.Hz        = Hz;   % Sampling Frequency
P.StartJump = 0;    % Number of shift samples to adjust at the beginning of the series
P.FirstScan = 1;
P.TickTime  = 1000/P.Hz;
[Y,T,P, data] = read_ppg(data_file, info_file, P);


% Add path for 'peakdet.m' - this is the version from Marco Burges, NOT the version in the PhLEM toolbox
addpath('peakdet');

% ---- Run 'peakdet' -----
% Inputs: pulse data vector, threshold for peak detection
signal = data{1,3}; % Extract pulse/resp info from 'data'
%th = 4000; % Pulse data peaks at 4095

[pks,dep,pid,did] = peakdet(signal, th);

time_stamps = data{1,1};
peak_times = time_stamps(pid);

% ----- Plot first few peaks to visually check the peak detection -----
figure;
plot(signal(1:10000));
y = zeros(30,1);
y(:,1) = th;
x = pid(1:30,1);
hold on;
plot(x,y,'--or');

end


function [signal,time,P, data] = read_ppg(datafile, info_file, P)
% Taken from 'phys_noise'

fclose('all');

fid=fopen(datafile);
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

