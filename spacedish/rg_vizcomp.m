cd '/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/space_inputs/comps'
load nwaydecomp.mat
%%
dish_dim = [8 8];
comp = 6;
max_chan = 9;

%close all
figure

% amplitude profile
subplot(1,3,1)
imagesc(reshape(nwaycomp.comp{comp}{1}, dish_dim))
title('Amplitude Profile')
colorbar

% frequency profile
subplot(1,3,2)
plot(nwaycomp.freq, nwaycomp.comp{comp}{2})
xlabel('Frequency (Hz)')
ylabel('Component Strength')
title('Frequency Profile')

% phase/time profile
subplot(1,3,3)
[temp max_freq] = max(nwaycomp.comp{comp}{2}); % get max freq
%imagesc(reshape(nwaycomp.comp{comp}{4}(:,max_freq), dish_dim)) %spacefsp
imagesc(reshape(nwaycomp.comp{comp}{4}-nwaycomp.comp{comp}{4}(max_chan), dish_dim)) %spacetime
colorbar
caxis([-0.1 0.1])
title('Phase Profile at Max Frequency')

%% normalization snippet
% Normalize sigma's for unicity, and inversly normalize D or Pkl
% normalize sigma such that the magnitude weighted mean per component is sigmacirc/2
A = comp{1};
S = comp{4};
% set radian conversion factor
radconv = ((2*pi)./sigmacirc);
% convert to complex domain
Scomp = exp(1i*S.*radconv);
% get mangitude weighted mean and add sigmacirc/2
magwmeanScomp = mean(Scomp.*A,1) ./ exp(1i*(sigmacirc/2)*radconv);
% normalize sigma
Scomp = Scomp ./ repmat(magwmeanScomp,[smodey(1) 1]);
% move sigma back to the real domain
S      = angle(Scomp);
S(S<0) = S(S<0) + (2*pi);
S      = S ./ radconv;

