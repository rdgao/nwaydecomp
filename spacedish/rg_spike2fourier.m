function FourierCoef = rg_spike2fourier(spikedata, fs, t_end, spikeintime, tr_info, fr_config)
%   spike data: spike time stamps in 1xn cells, where n is the # of cells
%   fs: sampling rate of spikes
%   t_end: end time of recording in seconds, if empty, default to max spike
%           time + 0.5 seconds
%   spikeintime: binary flag for whether spiketimes are in time or index
%   tr_info: trial info in fieldtrip format, used to cut epochs manually
%           dim = n_trial x 3 : [trial_start, trial_end, offset]
%   fr_config: fourier parameters:
%       - freqoi: 1xNfreq vector, frequencies at which to perform complex exponential convolution
%       - timwin: scalar, length in seconds of the complex exponentials
%       - method: 'sparseconv', 'segmentsparseconv' or 'segmentfft', method to use for convolving spike trains 
%           - sparseconv: use convolution of sparse matrices in their original size
%           - segmentsparseconv: same as sparseconv, but less memory intensive 
%               (cuts epochs into segments, taking care of overlap)
%           - segmentfft: convolution via multiplication in the frequency domain, more memory
%               efficient when using a long time window (e.g. >100ms)

% take the spiketime vectors and create sparse matrix representation
num_cell = length(spikedata);
cell_id = cell(size(spikedata));
for c = 1:num_cell
    % make a cell array of cell ids
    cell_id{c} = c*ones(length(spikedata{c}),1);
end


spike_id = cell2mat(cell_id');
spiketimes = cell2mat(spikedata');
if spikeintime
    % spike times recorded in time, convert to index
    spiketimes = round(spiketimes*fs);
end

if isempty(t_end)
    % no end time given, default to max time + 0.5 second
    t_end = max(spiketimes)+round(fs*0.5);
else
    % convert to index
    t_end = ceil(t_end*fs);
end

% create sparse matrix: 
%   (spike index, spike time, array of 1s, # of cell, total time)
sparse_spikes = sparse(spike_id, spiketimes, ones(size(spiketimes)), num_cell, t_end);

% make a time vector
t_vec = (1:t_end)/fs;

if ~isempty(tr_info)
    % now cut into trials
    num_epochs = size(tr_info,1);
    spikes_epoched = cell(1, num_epochs);
    for tr = 1:num_epochs
        % assume no offset in the trial info matrix       
        cur_idx = tr_info(tr,1):tr_info(tr,2)-1;
        spikes_epoched{tr} = sparse_spikes(:,cur_idx);
    end    
else
    spikes_epoched = {sparse_spikes};
end

% get fourier coefficients
FourierCoef = rmr_SPACEinput_spiketrain(spikes_epoched,fs,fr_config.freqoi,fr_config.timwin);
