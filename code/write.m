function [ success ] = write(fid,cellarray,color)
%WRITE_CELLS writes a cell array of strings to the file identified by fid
%   Copyright (c) 2013-2015 Sylvain Quoilin & Jessica Schrouff. 
%   All rights reserved.

if nargin == 2
    color = 'black';
elseif nargin == 3
    
else
    disp('Function write requires two or three arguments')
end
        

for i = 1:size(cellarray,1)
    fprintf(fid,[cellarray{i} char(10)]);
    %disp(cellarray{i})
    cprintf(color,cellarray{i})
    cprintf('\n')
end

success = true;

end