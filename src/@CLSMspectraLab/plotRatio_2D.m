function varargout = plotRatio_2D(obj,ha,id,Zidx,Tidx,ch1,ch2,normalizeIntFlag)
% PLOTRATIO_2D plots the intensity ratio of two different channels
%
% The ratio is Ch1/Ch2.
%
% Usage:
%   plotRatio_2D(obj,ha,id,Zidx,Tidx,ch1,ch2) with ch1 the indices of 
%   channel 1, ch2 the indices of channel 2, id the number
%   of the data set (index imageStack array), Zidx the index of the
%   Z-slices, and Tidx the index of the time slice.
%   The figure is plotted in the axis specified by the axis handle ha, or 
%   the current axis is taken if left empty.
%   The brightness in the image corresponds to the intensty of the channels
%   combined.
%
%   plotRatio(obj,ha,ch1,ch2,Zidx,normalizeIntFlag) id true, 
%   normalizeIntFlag keeps the intensity of the pixel in the output image
%   constant regardless of the intensty of the channels combined.
%   
%   hps = obj(ii).plotRatio(...) returns the handles to the
%   plot.
%
%   [hps, ha] = obj(ii).plotRatio(...) returns also the handle to
%   the axes handle of the plot.

if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 2 || isempty(ha)
    ha = gca;
end

if nargin < 3 || isempty(id)
    id = 1;
end

if nargin < 4 || isempty(Zidx)
    Zidx = 1;
end

if nargin < 5 || isempty(Tidx)
    Tidx = 1;
end

if nargin < 8
    normalizeIntFlag = false;
end

if isempty(ch1)
    warning('Channel 1 is empty.')
end

if isempty(ch2)
    warning('Channel 2 is empty.')
end

[Ich1,Ich2] = getIbinsProcessed(obj,id,'sum','c',{ch1,ch2},'x','y','z','t','c');

Ich1 = squeeze(Ich1(:,:,Zidx,Tidx,:)); % take z,t-slice
Ich2 = squeeze(Ich2(:,:,Zidx,Tidx,:));

Iratio = sum(Ich1,3:5)./sum(Ich2,3:5);

if normalizeIntFlag
    hps = imagesc(Iratio);
else
    hps = imagesc(Iratio,"AlphaData",rescale(sum(cat(3,Ich1,Ich2),3)));
end
axis image
colormap jet
hbar = colorbar;
title(hbar, 'ch1/ch2')

% Output
if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end