function plotRGBcolormap2D(obj,ha,id,Zidx,Tidx,ch1,ch2,ch3)
% PLOTRBGCOLORMAP2D plots the intensity of three different channels in red,
% green,and blue.
%
% Channel 1 is plotted in red, channel 2 is plotted in green, and channel 3
% is plotted in blue.
%
% Usage:
%
% plotRGBcolormap(obj,ha,id,Zidx,Tidx,ch1,ch2,ch3) with ch1 the indices of channel 1,
%  ch2 the indices of channel 2, ch3 the indices of channel 3, id the number
%  of the data set (index imageStack array), Zidx the index of the
%  Z-slices, and Tidx the index of the time slice.

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

if isempty(ch1)
    warning('Channel 1 is empty.')
end

if isempty(ch2)
    warning('Channel 2 is empty.')
end

if isempty(ch3)
    warning('Channel 3 is empty.')
end

[Ich1,Ich2,Ich3] = getIbinsProcessed(obj,id,'sum','c',{ch1,ch2,ch3},'x','y','z','t','c');

Irgb = zeros(size(Ich1,1),size(Ich1,2),3);

Ich1 = squeeze(Ich1(:,:,Zidx,Tidx,:)); % take z,t-slice
Ich2 = squeeze(Ich2(:,:,Zidx,Tidx,:));
Ich3 = squeeze(Ich3(:,:,Zidx,Tidx,:));

% red channel
Irgb(:,:,1) = rescale(sum(Ich1,3:5));

% green
Irgb(:,:,2) = rescale(sum(Ich2,3:5));

% blue
Irgb(:,:,3) = rescale(sum(Ich3,3:5));

imshow(Irgb,'Parent',ha)
axis image

end
