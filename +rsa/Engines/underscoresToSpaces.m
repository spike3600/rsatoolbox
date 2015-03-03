function stringORstringInCell=underscoresToSpaces(stringORstringInCell)
%
%  underscoresToSpaces is a function based on deunderscore.m
%  It takes an incoming string and replaces all underscores (which
%  are sometimes interpreted by MATLAB as subscript indicators in figures)
%  as spaces (which aren't).
%
%  Cai Wingfield 11-2009
%__________________________________________________________________________
% Copyright (C) 2009 Medical Research Council

import rsa.*

if iscell(stringORstringInCell)
    for lineI=1:numel(stringORstringInCell)
        line=stringORstringInCell{lineI};
        line(line==95)=' ';
        stringORstringInCell{lineI}=line;
    end
else
    stringORstringInCell(stringORstringInCell==95)=' ';
end  