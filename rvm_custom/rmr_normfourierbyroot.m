function [fourier, scaling] = rmr_normfourierbyroot(fourier, k)
% function [fourier, scaling] = rmr_normfourierbyroot(fourier, k)
%   scale fourier such that the average power, of each channel,
%   over epochs/freq of the CSD, is equal to its kth root
%
% fourier: 4D fourier data, [chan,freqs,trials,tapers]
% k: kth root (use power of 2 is fine, 8, 16, 32, etc)

power = zeros(size(fourier,1),1);
for iepoch = 1:size(fourier,3)
    currfour = double(squeeze(fourier(:,:,iepoch,:)));
    power = power + nansum(nansum(abs(currfour).^2,3),2);
end
power = power ./ (size(fourier,3) .* size(fourier,2));
% scale power by its kth root
scaling = power .^ (1/k);
for iepoch = 1:size(fourier,3)
    for ifreq = 1:size(fourier,2)
        currfour  = double(squeeze(fourier(:,ifreq,iepoch,:)));
        nonnanind = ~isnan(currfour(1,:));
        currfour  = currfour(:,nonnanind);
        currfour  = bsxfun(@rdivide,currfour,sqrt(power));
        currfour  = bsxfun(@times,currfour,sqrt(scaling));
        % save in fourier
        fourier(:,ifreq,iepoch,nonnanind) = single(currfour);
    end
end