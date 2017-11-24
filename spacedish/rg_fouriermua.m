% compute fourier coefficients for MUA of all recordings
% two ways to go about doing this:
% 1. grab the spikes and convolve with a wide gaussian taper
% 2. filter the LFP for high gamma range and compute the hilbert amplitude

fs = 1000;

% partition parameters
num_part = 24;
part_len = 20; %seconds

% fourier parameters
fr_cfg = [];
fr_cfg.freqoi = 0.5:0.25:10;
fr_cfg.timeoi = 0.001:0.001:part_len;
fr_cfg.timwin = 5./fr_cfg.freqoi; % 5 cycle filter
fr_cfg.taper = 'hanning';

%%
% handle file pathing
cd '/Volumes/My Passport for Mac/Dish/CTC/';
F = dir('CTC_*');
%reorder folders because 2017 recordings are in front
for f=1:length(F)
    if strcmp(F(f).name,'CTC_073116')
        offset = f;
    end
end
F = [F(offset:end);F(1:offset-1)];

%% traversing and compute
mua_method = 'conv_spikes';
wells= [5];
num_well = length(wells);
conv_win = gausswin(51);
conv_win = conv_win/sum(conv_win);

% record days with missing data since fourier calculation is skipped
good_data = zeros(length(F), num_well);

% initialize cells
FC_all = cell(1,12);

for f=1:length(F)
    cd(F(f).name)
    disp(F(f).name)
    
    switch mua_method
        case 'conv_spikes'
            % spike convolution case
            load LFP_Sp.mat spikes t_s
            bsp = binarize_spikes(t_s(end), 12500, spikes, fs);
            [~, numchan, len] = size(bsp);
            
            % grab time indices for cutting trials
            part_inds = round(linspace(1, len-fs*part_len, num_part)); % start index of windows
            tr_info = [part_inds' part_inds'+part_len*fs-1 zeros(num_part,1)];
            
            disp('Smoothing spikes & computing Fourier coefs...')
            for well=wells
                disp(well)
                data = zeros(size(squeeze(bsp(well, :, :))));
                for chan=1:numchan
                    data(chan,:) = conv(squeeze(bsp(well,chan,:)), conv_win, 'same');                    
                end                
                FC = rg_ts2fourier(data, fs, tr_info, fr_cfg, []);
                FC_padded = zeros(numchan, length(fr_cfg.freqoi), num_part, numchan);
                nonzeros = size(FC,4);
                FC_padded(:,:,:,1:nonzeros) = FC;
                FC_all{well} = cat(3, FC_all{well}, FC_padded);                
            end            
            
        case 'filt_lfp'
            %load LFP_Sp.mat LFP spikes spike_cnt t_s t_ds
            for well=5:5
                disp(sprintf('Well %i', well));
                
                %skip well if there's no LFP
                if isempty(LFP{well})
                    disp(sprintf('Well %i skipped, no data.', well))
                    continue
                end
                
                %skip well if there's missing channels
                adj = sum(1-all(LFP{well}==0))/64;
                if adj~=1
                    disp(sprintf('Well %i missing channels.',well))
                    continue
                end
                
                data = LFP{well}';
                [numchan, len] = size(data);
                
                
                part_inds = round(linspace(1, len-fs*part_len, num_part)); % start index of windows
                tr_info = [part_inds' part_inds'+part_len*fs-1 zeros(num_part,1)];
                
                FC = rg_ts2fourier(data, fs, tr_info, fr_cfg);
                FC_all{well-4} = cat(3, FC_all{well-4}, FC);
                good_data(f, well) = 1;
                clear FC
            end
            disp(size(FC_all{1}))
            disp(good_data)
    end
    cd ..
end

%% save out
clear fourier
cd '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/SPACE/MUA/'
for well = wells
    fourier = FC_all{well};
    save(sprintf('fourfull_%i.mat', well), 'fourier', '-v7.3');    
end

%% grab only the good trials
load '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/data_qual.mat';
fc_all = load('fourfull_5.mat');
%%
num_part = 24;
good_recs = find(quality_5==1);
tr_include = zeros(1,size(fc_all.fourier,3));
for rec = 1:length(good_recs)
    tr_include((good_recs(rec)-1)*num_part+(1:num_part))=1;
end
disp(all(unique(ceil(find(tr_include)/num_part))==good_recs')); %should be 1
fourier_clean = fc_all.fourier(:,:,find(tr_include),:); % grab good indices
%%
fourier = fourier_clean;
save('fourfull.mat', 'fourier', '-v7.3');
fourier = fourier_clean(:,:,1:2:end,:);
save('fourprt1.mat', 'fourier', '-v7.3');
fourier = fourier_clean(:,:,2:2:end,:);
save('fourprt2.mat', 'fourier', '-v7.3');

%% grab all spike counts
spike_cnt_all = zeros(length(F),12,64);
for f=1:length(F)
    cd(F(f).name)
    disp(F(f).name)
    
    load LFP_Sp.mat spikes t_s spike_cnt
    spike_cnt_all(f,:,:) = spike_cnt;
    cd ..
end
%% 
save('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/SPACE/MUA/spkcnt.mat', 'spike_cnt_all');