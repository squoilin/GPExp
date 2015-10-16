function [ success ] = write(fid,data,color)
%WRITE_CELLS writes a cell array of strings to the file identified by fid
%   Copyright (c) 2013-2015 Sylvain Quoilin & Jessica Schrouff. 
%   All rights reserved.

if nargin == 2
    color = 'black';
elseif nargin == 3
    
else
    disp('Function write requires two or three arguments')
end
        
if iscell(data)
    for i = 1:size(data,1)
        fprintf(fid,[data{i} char(10)]);
        %disp(cellarray{i})
        cprintf(color,data{i})
        cprintf('\n')
    end
elseif ischar(data)
    for i = 1:size(data,1)
        fprintf(fid,[data(i,:) char(10)]);
        cprintf(color,data(i,:))
        cprintf('\n')
    end
else
    error('The input of the write function must be a cell array or a char array')
end

success = true;

end