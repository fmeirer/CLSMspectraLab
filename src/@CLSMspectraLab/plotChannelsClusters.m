function varargout = plotChannelsClusters(obj,ha,indices)
% PLOTCHANNELSCLUSTERS plots the average intensity per channel in a cluster.
%
%   Usage:
%   plotChannels(obj,ha) plots the average intensity per channel over
%   the whole sample.
%
%   plotChannels(obj,ha,indices) plots the average intensity per channel over
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


channelVals = obj.getChannelsClusters(indices);
hps = [];
kk = 1;
legNames = {};
for ii = 1:numel(indices)
    hps = [hps; gobjects(obj.clustering(indices(ii)).nClusters,1)];
    for jj = 1:obj.clustering(indices(ii)).nClusters
        hps(kk) = plot(ha,(1:size(channelVals{ii},1)),channelVals{ii}(:,jj),'.-');
        if numel(indices) == 1
            legNames{kk} = sprintf('Cluster %s',num2str(jj));
        else
            legNames{kk} = sprintf('Data %s, cluster %s',num2str(indices(ii)),num2str(jj));
        end
        kk = kk + 1;
        hold on
    end
end

xlabel(ha,'Channel')
ylabel(ha,'Mean intensity')
legend(legNames{:})

% Output
if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

