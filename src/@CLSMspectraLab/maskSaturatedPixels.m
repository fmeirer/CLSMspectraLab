function obj = maskSaturatedPixels(obj,threshold)
% MASKSATURATEDPIXELS sets pixels in mask higher than threshold to 0
%
%   Usage:
%   obj = maskSaturatedPixels(obj,threshold) with values in image equal or
%   larger than threshold are masked out by setting mask == false for these
%   voxels.
    

for ii = 1:numel(obj)
    for jj = 1:numel(obj(ii).input)
        % get image with maximum value in masking dimensions
        I = obj(ii).input(jj).getReducedImage('mean',obj(ii).mask(jj).dimLabel{:});
        % do masking
        obj(ii).mask(jj) = obj(ii).mask(jj).setValue(I >= threshold,false);
    end
end

end