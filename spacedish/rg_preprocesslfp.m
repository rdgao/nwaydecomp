%%% 
% preprocess LFP data following SPACE-FSP tutorial
%%%

%% load data
datafile = '~/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_110416/LFP_Sp.mat';
load(datafile)
%%
infofile = '~/Documents/Data/Muotri/Pri_Corticoids/aggregate.mat';
load(infofile,'kerTimes')

%% manually setting up FT structure
day = 18;
well = 5;

% data info
data_raw = LFP{well}';
fs_ds = 1000;
[numchan,datalen]=size(data_raw);
ev_t = kerTimes{day}{well};

% setting up data struct
data = [];
data.label=cell(numchan,1);
for chan=1:numchan
    data.label{chan} = sprintf('%i',chan);
end
data.fsample=fs_ds;
data.trial={data_raw};
data.time={t_ds'};
data.sampleinfo=[1 datalen];

% trial info
cfg=[];
cfg.trl = [ev_t-0.5 ev_t+2.5 zeros(length(ev_t),1)]*1000;
data_tr = ft_redefinetrial(cfg, data);
%% preprocessing
data_tr = ft_preprocessing(cfg, data_tr);

%% computing SPACE input
cfg = []; 
cfg.foi = 5:0.5:10;      
cfg.t_ftimwin = 3./cfg.foi; % sliding windows of 5 cycles wide
cfg.toi = 0.5:0.01:2.5; % time-window slides from 0 to 1.5s in steps of 10ms
cfg.method = 'mtmconvol';  %wavelet convolution
cfg.taper= 'hanning'; %window function for wavelets
cfg.keeptrials = 'yes';
cfg.output = 'fourier';
cfg.pad = 10; % pad the data out to 2s, using spectral interpolation to obtain 
                % Fourier coefficients at integer frequencies
freqdata = ft_freqanalysis(cfg, data_tr);

%% extract rhythmic components
cfg = [];
cfg.model = 'spacefsp';
cfg.datparam = 'fourierspctrm';
cfg.Dmode = 'identity';
cfg.ncompest = 'splitrel';
cfg.ncompeststart = 3;
cfg.ncompeststep = 2;
cfg.ncompestend = 20;
cfg.ncompestsrdatparam = 'oddeven';
cfg.ncompestsrcriticalval = 0.7;
cfg.numiter = 1000;
cfg.convcrit = 1e-6;
cfg.randstart=10;
cfg.ncompestrandstart = 10;
cfg.fsample = data_tr.fsample;
nwaydecomp = nd_nwaydecomposition(cfg,freqdata);


%% --- alternate way of getting input
%% computing SPACE input
freqoi = 5:0.5:20;
% timwin = ones(size(freqoi))*0.5;
timwin = 5./freqoi
timeoi = 0.5:0.01:2.5;
taper = 'hanning';
FCoef = rmr_SPACEinput_electrophysiology(data_tr, freqoi, timeoi, timwin, taper);
FData = freqdata;
FData.fourierspctrm = FCoef;
FData.dimord = 'chan_freq_epoch_tap';
%%
cfg = [];
cfg.model = 'spacefsp';
cfg.datparam = 'fourierspctrm';
cfg.Dmode = 'identity';
cfg.ncompest = 'splitrel';
cfg.ncompeststart = 3;
cfg.ncompeststep = 2;
cfg.ncompestend = 20;
cfg.ncompestsrdatparam = 'oddeven';
cfg.ncompestsrcriticalval = 0.7;
cfg.numiter = 1000;
cfg.convcrit = 1e-6;
cfg.randstart=10;
cfg.ncompestrandstart = 10;
cfg.fsample = data_tr.fsample;
nwaydecomp = nd_nwaydecomposition(cfg,FData);

%% -------------
%% equal partition %%%
%% ------------
% partition parameters
num_part = 12;
part_len = 40; %seconds
well = 5;
fs_ds = 1000;

%% load data
datafile = '~/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_110416/LFP_Sp.mat';
load(datafile)
%% set up 'epoch' LFP data
data = LFP{well};
[len,numchan] = size(data);
part_inds = round(linspace(1, len-fs_ds*part_len, num_part)); % start index of windows

% setting up data struct
data = [];
data.label=cell(numchan,1);
for chan=1:numchan
    data.label{chan} = sprintf('%i',chan);
end
data.fsample=fs_ds;
data.trial={LFP{well}'};
data.time={t_ds'};
data.sampleinfo=[1 len];

% trial info
cfg=[];
cfg.trl = [part_inds' part_inds'+part_len*fs_ds-1 zeros(num_part,1)];
data_tr = ft_redefinetrial(cfg, data);
data_tr = ft_preprocessing(cfg, data_tr);

%% computing SPACE input
freqoi = 0.5:0.5:20;
timwin = 7./freqoi; % 7 cycle filter
timeoi = 1:0.001:part_len-1; % cut off 1 second from beginning and end
taper = 'hanning';
FCoef = rmr_SPACEinput_electrophysiology(data_tr, freqoi, timeoi, timwin, taper);
disp(size(FCoef))
%% save out to disk
clear fourier
datafolder = '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/space_inputs/';
% odd trials
fourier = FCoef(:,:,1:2:num_part,:);
save([datafolder 'fourprt1.mat'], 'fourier', '-v7.3');
% even trials
fourier = FCoef(:,:,2:2:num_part,:);
save([datafolder 'fourprt2.mat'], 'fourier', '-v7.3');
% full trials
fourier = FCoef;
save([datafolder 'fourfull.mat'], 'fourier', '-v7.3');
% config
labels = data.label;
save([datafolder 'space_configs.mat'], 'labels','freqoi', '-v7.3');

%%







%%
%% DO NOT USE BECAUSE FOURIER COEFFICIENTS COME OUT THE WRONG ORDER 
% AND FIELDTRIP DOES A MYSTERY INTERNAL NORMALIZATION THING THAT WE SHOULD
% UNDO
%% compute Fourier coefficients
cfg = []; 
cfg.foi = 0.5:0.5:20;
cfg.t_ftimwin = 7./cfg.foi; % sliding windows of 5 cycles wide
cfg.toi = 1:0.01:part_len-1; % time-window slides from 1 to end-1 in steps of 10ms
cfg.method = 'mtmconvol';  %wavelet convolution
cfg.taper= 'hanning'; %window function for wavelets
cfg.keeptrials = 'yes';
cfg.output = 'fourier';
cfg.pad = 40; % pad the data out to 2s, using spectral interpolation to obtain 
                % Fourier coefficients at integer frequencies
freqdata = ft_freqanalysis(cfg, data_tr);
%% extract rhythmic components
cfg = [];
cfg.model = 'spacefsp';
cfg.datparam = 'fourierspctrm';
cfg.Dmode = 'identity';
cfg.ncompest = 'splitrel';
cfg.ncompeststart = 3;
cfg.ncompeststep = 2;
cfg.ncompestend = 20;
cfg.ncompestsrdatparam = 'oddeven';
cfg.ncompestsrcriticalval = 0.7; [0.7 0.7 0 0.7 0] %space, freq, trial, time, matrixD
cfg.numiter = 1000;
cfg.convcrit = 1e-6;
cfg.randstart=10;
cfg.ncompestrandstart = 10;
cfg.fsample = data_tr.fsample;
nwaydecomp = nd_nwaydecomposition(cfg,freqdata);


%%


%%
    
    



%%

FData = freqdata;
FData.fourierspctrm = FCoef;
FCoef = rmr_SPACEinput_electrophysiology(data_tr, freqoi, timeoi, timwin, taper);
% fake the FT struct
FData.labels = data_tr.label;
FData.freq = freqoi;
FData.time=timeoi;
FData.fourierspctrm = FCoef;
FData.dimord = 'chan_freq_epoch_tap';

% all the stats for all the component tries
% nwaydecomp.splitrelstat.allrandomstatfull

% .cong{rall, .., ..} 
% congrall : congruence over initializations
% check to see that minimum is the same between initializations
% congglobmin : congruence over only the thresholded global mins
% congrcumul : does the same as the above two but for progressing number of
%   pairs, ordered by max explained variance
% --> NEVERMIND: pairs that always include the best one from 2:n
%   look at this to see if we find crap components

% split-trial reliability may or may not be an appropriate metric to use
% with this, we ensure that the components found are definitely there, but
% not necessarily capture all the components





