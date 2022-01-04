classdef (Abstract) abstract < imageStack 
    %IMAGEMASK Abstract class to compute a background-corrected image
    
    properties (Abstract,Constant)
        name; % name of image mask computation method
    end
    
    methods (Abstract)
        obj = compute(obj,imageStack); % compute the background correction
    end
    
    methods

    end
    
end

