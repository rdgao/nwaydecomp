% compute fourier coefficients for spikes
fs = 12500;

% partition parameters
num_part = 40;
part_len = 10; %seconds

% fourier parameters
fr_cfg = [];
fr_cfg.timwin = 0.02; % 20 ms window
fr_cfg.freqoi = (1:20)./fr_cfg.timwin;

%% 
% define files to run on
datafolder = '/Volumes/My Passport for Mac/Dish/CTC/';
files = {'CTC_121716','CTC_121716_Drugs'};
%% 
% run it
wells = [6 8 10 12];
num_well = length(wells);
FC_all = cell(1,num_well);
for file = 1:length(files)
    filename = [datafolder, files{file}];
    disp(filename)
    load([filename '/LFP_Sp.mat'], 'spikes', 't_s')

    len = length(t_s);
    part_inds = round(linspace(1, len-fs*part_len, num_part)); % start index of windows
    tr_info = [part_inds' part_inds'+part_len*fs-1 zeros(num_part,1)];
        
    for well = wells
        disp(well)
        FC = rg_spike2fourier(spikes(well,:), 12500, ceil(t_s(end)), 0, tr_info, fr_cfg);
        FC_all{find(wells==well)} = cat(3, FC_all{find(wells==well)}, FC);
    end
end
%%
clear fourier
cd '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/SPACE/Spikes'

for well = 1:num_well
    fourier = FC_all{well};
    save(sprintf('fourfull_%i.mat', wells(well)), 'fourier', '-v7.3');    
end

%% code snippet for visualizing total power at a frequency across channels
freq = 1;
CS = zeros(64,64);
for tr = 1:80
    idx = find(~isnan(squeeze(fourier(1,freq,tr,:))));
    CS = CS + squeeze(fourier(:,freq,tr,idx))*squeeze(fourier(:,freq,tr,idx))';
end
figure
imagesc(abs(CS))
colorbar

