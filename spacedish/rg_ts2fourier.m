function FourierCoef = rg_ts2fourier(rawdata, fs, tr_info, fr_config)
% rg_ts2fourier(rawdata, fs, tr_info, fr_config)
%   rawdata: time series, [channel x samples]
%   fs: sampling frequency
%   tr_info: trial info in the required fieldtrip format 
%           dim = n_trial x 3 : [trial_start, trial_end, offset]
%   fr_config: fourier parameters:
%       - freqoi: frequency of interest (Hz)
%       - timeoi: time window of interest (within trial, in seconds)
%       - timwin: taper length, in seconds
%       - taper: taper type ('hanning')

[numchan, len] = size(rawdata);

data = [];
% provide data labels
data.label=cell(numchan,1);
for chan=1:numchan
    data.label{chan} = sprintf('%i',chan);
end
% filling in ft data struct parameters
data.fsample=fs;
data.trial={rawdata};
data.time={(0:len-1)/fs};
data.sampleinfo=[1 len];

% trial info
%[part_inds' part_inds'+part_len*fs_ds-1 zeros(num_part,1)];
cfg=[];
cfg.trl = tr_info;

data_tr = ft_redefinetrial(cfg, data);
data_tr = ft_preprocessing(cfg, data_tr);
FourierCoef = rmr_SPACEinput_electrophysiology(data_tr, fr_config.freqoi, fr_config.timeoi, fr_config.timwin, fr_config.taper);
