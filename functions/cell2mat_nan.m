function [Out] = cell2mat_nan(C)

maxLength = max(cellfun(@(x)numel(x),C));
Out = cell2mat(cellfun(@(x)cat(2,x,nan(1,maxLength-length(x))),C,'UniformOutput',false));
