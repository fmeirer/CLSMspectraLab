function varargout = plotMask(obj,ha)
% PLOTMASK plots the mask.
%
%   Usage:
%   plotMask(obj) plots the mask.
%
%   plotMask(obj,ha) plots in the axis specified by the axis handle
%   ha.
%
%   hps = obj(ii).plotMask(...) returns the handles to the
%   plot.
%
%   [hps, ha] = obj(ii).plotMask(...) returns also the handle to
%   the axes handle of the plot.


if nargin < 2 || isempty(ha)
    ha = gca;
end

z = obj.stackSize;
[nR,nC] = subplotDimensions(z);

if ~obj.isVolume
    hps = gobjects(z,1);
    for ii=1:z
       subplot(nR,nC,ii); 
       hps(ii) = imagesc(obj.Imask(:,:,ii)); 
       axis image;
    end
else
    voxelSize = obj.getVoxelSize;
    voxelSize = voxelSize./voxelSize(1); % normalize to X
    voxelSize(isnan(voxelSize)) = 1; % set al unknown values to 1
    input('Displaying volume. Ready to continue? Press Enter... (it is a good idea to first close the volume viewer)','s');
    volumeViewer(obj.Imask,'ScaleFactors',voxelSize);
    hps = [];
end

% Output
if nargout > 0
    varargout{1} = hps;
    if nargout > 1
        varargout{2} = ha;
    end
end

end

