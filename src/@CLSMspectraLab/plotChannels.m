function varargout = plotChannels(obj,ha,indices)
% PLOTCHANNELS plots the average intensity per channels.
%
%   Usage:
%   plotChannels(obj,ha) plots the average intensity per channel over
%   the whole sample.
%
%   %   plotChannels(obj,ha,indices) plots the average intensity per channel over
%   the whole sample for the specified indices
%
%   hps = obj(ii).plotChannels(...) returns the handles to the
%   plot.
%
%   [hps, ha] = obj(ii).plotChannels(...) returns also the handle to
%   the axes handle of the plot.


if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 3 || isempty(indices)
    indices = 1:obj.nInput;
end

if max(indices) > obj.nInput
    error('Not all specified indices do exist.')
end

hps = gobjects(numel(indices),1);
for ii = indices
    channelInt = getImageProcessed(obj,ii,'mean','c','x','y','z','t'); % in order to have the mask applied, we cannot directly request 'c'
    channelIntVec = nanmean(channelInt,2:5);
    hps(ii) = plot(ha,(1:numel(channelIntVec)),channelIntVec,'.-');
    hold on
end

xlabel(ha,'Channel')
ylabel(ha,'Mean intensity')
legend(string(indices))

% Output
if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

