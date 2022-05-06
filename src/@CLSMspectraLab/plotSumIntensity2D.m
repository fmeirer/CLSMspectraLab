function varargout = plotSumIntensity2D(obj,ha,id,Zidx,Tidx)
% PLOTSUMINTENSITY2D plots the summed intensity over all channels.
%
%   Usage:
%   plotSumIntensity2D(obj) plots the summed intensity of the first slice in
%   with no binning in the current axis.
%
%   plotSumIntensity2D(obj,ha) plots in the axis specified by the axis handle
%   ha.
%
%   plotSumIntensity2D(obj,ha,id) plots the id-th image. If is
%   empty, the first image is taken.
%
%   plotSumIntensity2D(obj,ha,id,Zidx) plots the Z-slice with index Zidx. If is
%   empty, the first slice is taken.
%
%   plotSumIntensity2D(obj,ha,id,Zidx,Tidx) plots the T-slice with index Tidx. If is
%   empty, the first slice is taken.
%
%   hps = obj(ii).plotSumIntensity2D(...) returns the handles to the
%   plot.
%
%   [hps, ha] = obj(ii).plotSumIntensity2D(...) returns also the handle to
%   the axes handle of the plot.


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

I = obj.getImageProcessed(id,'sum','x','y','z','t');

imagesc(ha,sum(I(:,:,Zidx,Tidx),3:4))
axis image

% Output
if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

