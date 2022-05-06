function [Ichannel,legNames] = getChannelsROI(obj,ROItable,meanFlag)
% GETCHANNELSROI gets the average intensity per channels.
%
%   Usage:
%   [channelInt,legNames] = getChannelsROI(obj,ROItable,meanFlag) 
%   outputs the average intensity per ROIdefined in the ROItable. 
%   channelInt is a column cell array with the spectrum per ROI. The
%   variable legNames is a column cell array containing the channel names
%   as defined per their label. The ROItable contains at least the fields 
%   'data', 'z', 't', 'coords', 'label'
%   which are respectively an integer with the id of the data set in the 
%   obj; integer with the number of the z-slice;  integer with the number 
%   of the t-slice, an array with polygon coordinates, and an integer id
%   of the ROI.
%
%   [channelInt,legNames] = getChannelsROI(obj,ROItable,meanFlag)
%   if meanFlag is true, all ROIs are averaged weigthed by the number of
%   pixels per ROI. channelInt is a vector and legNames are empty.

if isempty(ROItable)
    Ichannel = [];
    return
end

if nargin < 3
    meanFlag = false;
end

uniqueData = unique(ROItable.data);
Ichannel = cell(size(ROItable,1),1);
legNames = cell(size(ROItable,1),1);
nPixels = nan(size(ROItable,1),1);
kk = 1;

for ii = 1:numel(uniqueData) % loop first over the unique data set and open them one by one
    ROItableThisData = ROItable(ROItable.data == uniqueData(ii),:);
    Ir = getImageProcessed(obj,uniqueData(ii),'reshaped','x','y','z','t','c');
    
    for jj = 1:size(ROItableThisData,1) % loop over entries for ii-th data set
        % get idx and coords from table
        z = ROItableThisData.z(jj);
        t = ROItableThisData.t(jj);
        polycoords = ROItableThisData.coords{jj};
        % get ROI slice 
        I = squeeze(Ir(:,:,z,t,:));
        temp_mask = poly2mask(polycoords(:,1),polycoords(:,2),...
            obj.input(uniqueData(ii)).getDim('x'),obj.input(uniqueData(ii)).getDim('y'));
        nPixels(kk) = sum(temp_mask,'all');
        Ic = bsxfun(@times, I, cast(temp_mask, 'like', I));
        % average over x and y, Ichannel{ii} is a column vector with the mean values per channel of ROI 
        Ichannel{kk} = squeeze(nanmean(Ic,[1 2])); 
        legNames{kk} = sprintf('Region %d',ROItableThisData.label(jj));
        kk = kk + 1; % current roi counter
    end
end


if meanFlag
    nChannels = cellfun(@numel,Ichannel);
    if numel(unique(nChannels)) ~= 1
        error('The variable ''meanFlag'' can only be set to true if all ROIs have the same number of channels.')
    end
    nPixels = nPixels./sum(nPixels); % weight for mean
    Ichannel = horzcat(Ichannel{:}) .* nPixels(:)'; % multiply ROI spectrum by weight
    Ichannel = sum(Ichannel,2);
    legNames = [];
end

end

