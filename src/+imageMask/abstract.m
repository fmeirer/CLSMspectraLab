classdef (Abstract) abstract < imageStack 
    %IMAGEMASK Abstract class to compute an image mask
    
    properties (Abstract,Constant)
        name; % name of image mask computation method
    end
    
    methods (Abstract)
        obj = compute(obj,imageStack); % compute the mask
    end
    
    methods

    end
    
end

