function [ success ] = write(fid,cellarray)
%WRITE_CELLS writes a cell array of strings to the file identified by fid
%   Copyright (c) 2013-2015 Sylvain Quoilin & Jessica Schrouff. 
%   All rights reserved.

for i = 1:size(cellarray,1)
    fprintf(fid,[cellarray{i} char(10)]);
    disp(cellarray{i})
end

success = true;

end