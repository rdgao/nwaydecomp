% compute fourier coefficients for all recordings
fs = 1000;

% partition parameters
num_part = 24;
part_len = 20; %seconds
filt_param = [0 30]; %lowpass filter

% fourier parameters
fr_cfg = [];
fr_cfg.freqoi = 0.5:0.25:10;
fr_cfg.timeoi = 0.001:0.001:part_len; % cut off 1 second from beginning and end
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
num_well = 1;

% record days with missing data since fourier calculation is skipped
good_data = zeros(length(F), num_well);

% initialize cells
FC_all = cell(1,num_well);
tr_all = cell(length(F),1);

for f=1:length(F)
    cd(F(f).name)
    disp(F(f).name)
    load LFP_Sp.mat LFP t_s t_ds
    
    % cut trials
    len = length(t_ds);
    part_inds = round(linspace(1, len-fs*part_len, num_part)); % start index of windows
    tr_info = [part_inds' part_inds'+part_len*fs-1 zeros(num_part,1)];        
    tr_all{f} = tr_info;
    
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
        FC = rg_ts2fourier(data, fs, tr_info, fr_cfg, filt_param);
        FC_all{well-4} = cat(3, FC_all{well-4}, FC);
        good_data(f, well) = 1;
        clear FC
    end
    disp(size(FC_all{1}))
    %disp(good_data)
    cd ..
end

%% save out
clear fourier
cd /Users/rgao/Documents/Data/Muotri/Pri_Corticoids/space_inputs/all_lfp_FC
for well = 1:num_well
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
%% grabbing indices
num_well = 1;

% record days with missing data since fourier calculation is skipped
good_data = zeros(length(F), num_well);

% initialize cells
tr_all = cell(length(F),1);

for f=1:length(F)
    cd(F(f).name)
    disp(F(f).name)
    load LFP_Sp.mat t_s t_ds
    len = length(t_ds);
    part_inds = round(linspace(1, len-fs*part_len, num_part)); % start index of windows
    tr_info = [part_inds' part_inds'+part_len*fs-1 zeros(num_part,1)];
    tr_all{f} = tr_info;
    cd ..
end
%%
save trial_indices.mat tr_all