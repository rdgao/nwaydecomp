% compute fourier coefficients for all recordings
fs = 1000;

% partition parameters
num_part = 16;
part_len = 30; %seconds

% fourier parameters
fr_cfg = [];
fr_cfg.freqoi = 0.5:0.5:12;
fr_cfg.timeoi = 1:0.001:part_len-1; % cut off 1 second from beginning and end
fr_cfg.timwin = 7./fr_cfg.freqoi; % 7 cycle filter
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

for f=15:length(F)
    cd(F(f).name)
    disp(F(f).name)
    load LFP_Sp.mat LFP spikes spike_cnt t_s t_ds

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
        
        keyboard
        
        part_inds = round(linspace(1, len-fs*part_len, num_part)); % start index of windows
        tr_info = [part_inds' part_inds'+part_len*fs-1 zeros(num_part,1)];

        FC = rg_ts2fourier(data, fs, tr_info, fr_cfg);
        FC_all{well-4} = cat(3, FC_all{well-4}, FC);
        good_data(f, well) = 1;
        clear FC
    end
    disp(size(FC_all{1}))
    disp(good_data)
    cd ..
end

%% save out
clear fourier
cd /Users/rgao/Documents/Data/Muotri/Pri_Corticoids/space_inputs/all_lfp_FC
for well = 1:num_well
    fourier = FC_all{well};
    save(sprintf('fourfull_%i.mat', well), 'fourier', '-v7.3');    
end