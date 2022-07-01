function objOut = getBinnedObj(obj,binFactor)
%GETBINNEDOBJ Bins input image.
%   Function bins the input images and outputs an empty CLSMspectraLab
%   object with the binned images.
%
%   Usage:
%   objOut = getBinnedObj(obj,binFactor): bins the imageStack(s) in obj
%   with a binFactor and outputs an empty CLSMspectraLab object with the
%   binned image in X and Y.

labels = {'x','y','z','t','c'};
iStack_out = imageStack.empty(obj.nInput,0);
for ii = 1:obj.nInput
    iStack_out(ii) = obj.input(ii).binImageFirst2Dims(binFactor,labels{:});
end

objOut = CLSMspectraLab(iStack_out,obj.binFactor*binFactor);

end

