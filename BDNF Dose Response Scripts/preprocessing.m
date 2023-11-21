%% Generates time series of spiking activity for specified experiment

%% Load Experiment and Set some parameters
dur      = 300;               % Recording length (s)
fs       = 20000;             % Sampling rate (Hz)
if ~exist('bin_size','var')
    bin_size = 0.002;         % Bin size (s)
end
% exp##_MEAinfo(##).div##.spikeMatrix specifies a particular recording
% Modify as necessary to select a different recording
% Spike Times are in msec
if ~exist('spkMat','var')
    load('MEAdata_BDNFdoseResp.mat');
    spkMat = (1e-3)*exp1_MEAinfo(1).div07.spikeMatrix;
end

%% Time Series
spks = zeros(64, dur/bin_size);
for ii=1:64
    if ismember(ii, [1,5,8,57,64])
        continue;
    end
    
    tspk = spkMat(:,ii); tspk(tspk==0)=[];
    ntmp = zeros(1,dur/bin_size);
    for jj=1:dur/bin_size
        ntmp(jj) = length(find(and(tspk>(jj-1)*bin_size,tspk<=jj*bin_size)));
    end
    spks(ii,:) = ntmp;
end
% spks = double(spks>0);
