function [numRow,numCol] = subplotDimensions(numPlots)
% [numRow,numCol] = subplotDimensions(numPlots)
%--------------------------------------------------------------------------
% Description
%
% Gives the preferred number of rows and coloumns used in a subplot based
% on the number of plots.
%
%--------------------------------------------------------------------------
% Necessary Inputs (name/data type)
% 
% numPlots/integer: number of plots
%
%--------------------------------------------------------------------------
% Outputs (name/data type)
%
% numRow/integer: number of rows for subplot
% numCol/integer: number of columns for subplot
%
%--------------------------------------------------------------------------
% Variable Inputs (flag/ data type /(default)):
% 
% 
%--------------------------------------------------------------------------
% Dependencies
%
%
%--------------------------------------------------------------------------
% Erik Maris
% 11/2018
% MATLAB 2018a
% calculates the the dimensions for a subplot with variable numel(input)

if numPlots == 1
    numRow = 1;
    numCol = 1;
elseif numPlots == 6
    numRow = 2;
    numCol = 3;
elseif numPlots == 8
    numRow = 2;
    numCol = 4;
else
    goldenRatio = 1.618;
    numCol = round(goldenRatio*sqrt(numPlots));
    numRow = ceil(numPlots/numCol);
end

end