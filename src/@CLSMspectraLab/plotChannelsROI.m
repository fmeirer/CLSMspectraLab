function varargout = plotChannelsROI(obj,ha,ROItable,meanFlag)
% PLOTCHANNELS plots the average intensity per channels.
%
%   Usage:
%   plotChannelsROI(obj,ROItable,meanFlag) 
%   plots the average intensity per ROIdefined in the ROItable. 
%   channelInt is a column cell array with the spectrum per ROI.
%   The ROItable contains at least the fields 
%   'data', 'z', 't', 'coords', 'label'
%   which are respectively an integer with the id of the data set in the 
%   obj; integer with the number of the z-slice;  integer with the number 
%   of the t-slice, an array with polygon coordinates, and an integer id
%   of the ROI.
%
%   plotChannelsROI(obj,ROItable,meanFlag)
%   if meanFlag is true, all ROIs are averaged weigthed by the number of
%   pixels per ROI.
%
%   hps = obj(ii).plotChannelsROI(...) returns the handles to the
%   plot.
%
%   [hps, ha] = obj(ii).plotChannelsROI(...) returns also the handle to
%   the axes handle of the plot.


if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 4
    meanFlag = false;
end

[Ichannel,legNames] = obj.getChannelsROI(ROItable,meanFlag);

if meanFlag
    hps = plot(ha,(1:numel(Ichannel)),Ichannel,'.-');
else
    hps = gobjects(numel(Ichannel),1);
    for ii = 1:numel(Ichannel)
        hps(ii) = plot(ha,(1:numel(Ichannel{ii})),Ichannel{ii},'.-');
        hold on
    end
    legend(legNames{:})
end
xlabel(ha,'Channel')
ylabel(ha,'Mean intensity')

% Output
if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

