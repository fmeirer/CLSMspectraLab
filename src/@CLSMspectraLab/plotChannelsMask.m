function varargout = plotChannelsMask(obj,ha)
% PLOTCHANNELSMASK plots the average intensity per channels of the
% foreground and background mask.
%
%   Usage:
%   plotChannels(obj,ha,mask) plots the average intensity per channel over
%   the foreground and background of the mask.
%
%   hps = obj(ii).plotChannels(...) returns the handles to the
%   plot.
%
%   [hps, ha] = obj(ii).plotChannels(...) returns also the handle to
%   the axes handle of the plot.


if nargin < 2 || isempty(ha)
    ha = gca;
end

hps = gobjects(obj.nInput,2);
c = lines(obj.nInput);
for ii = 1:obj.nInput
    channelInt = getImageProcessed(obj,ii,'reshape','c','x','y','z','t'); % in order to have the mask applied, we cannot directly request 'c'
    channelIntVec = nanmean(channelInt,2:5);
    hps(ii,1) = plot(ha,(1:numel(channelIntVec)),channelIntVec,'.-','Color',c(ii,:));
    hold on
    
    obj.mask(ii) = obj.mask(ii).invertImage; % invert to get background
    
    channelInt = getImageProcessed(obj,ii,'reshape','c','x','y','z','t'); % in order to have the mask applied, we cannot directly request 'c'
    channelIntVec = nanmean(channelInt,2:5);
    hps(ii,2) = plot(ha,(1:numel(channelIntVec)),channelIntVec,'.--','Color',c(ii,:));
    obj.mask(ii) = obj.mask(ii).invertImage; % invert back
end

xlabel(ha,'Channel')
ylabel(ha,'Mean intensity')
l = legend(hps(:,1),string(1:size(hps,1)));
title(l,'data')
title('-: foreground, --: background')

% Output
if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

