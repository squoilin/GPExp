function [in,errors,warnings] = inputs_sanity_check(in)
% This function performs a sanity check on the input structure variable of
% GPExp. It returns error messages to assist the user in the model
% parametrization process.
% 
% University of Liège, December 2014
% 
% Copyright (c) 2013-2015 Sylvain Quoilin & Jessica Schrouff. 
% All rights reserved.

errors = cell(0);
warnings = cell(0);

if nargin<1
    error('Please enter a valid input structure')
end

if ~isfield(in,'data')
    errors{end+1} = 'The "in" structure should at least contain the "data" field';
end

if strcmp(in.filename(end-2:end),'csv') || strcmp(in.filename(end-2:end),'xls') 
    in.name = in.filename(1:end-4);
elseif strcmp(in.filename(end-3:end),'xlsx')
    in.name = in.filename(1:end-5);
else
    in.name = in.filename;
end


if ~isfield(in,'headers') % variable names
    in.headers=cell(size(in.data,2),1);
    for i=1:size(in.data,2)
        in.headers{i} = ['Column ',num2str(i)];
    end
else
    if length(in.headers) ~= size(in.data,2)
        errors{end+1} = ['The number of headers (' num2str(length(in.headers)) ') does not match with the number of data columns (' num2str(size(in.data,2)) ')'];
    end
end

list_headers = [];
if size(in.data,2) == length(in.headers)
    for i=1:size(in.data,2)
        list_headers = [list_headers char(10) in.headers{i}];
    end
end

% Checking output and output name:
if ~isfield(in,'considered_output')
    errors{end+1} = ['Please specify the output to be considered in the "considered_output" variable. Possible values are: ' list_headers];
else
    if iscell(in.considered_output) && size(in.considered_output,1)>1
        errors{end+1} = 'Predicting one target (output) at a time, please execute multiple times for multiple targets';
    end
    idx_output = find(strcmp(in.headers,in.considered_output));
    if isempty(idx_output)
        errors{end+1} = ['The specified output (' in.considered_output ') could no be found in the list of headers:' list_headers];
    end
    in.y = in.data(:,idx_output);
end
        
% Checking inputs and inputs names:
if ~isfield(in,'considered_inputs')
    errors{end+1} = ['Please specify the inputs to be considered in the "considered_inputs" variable. Possible values are: ' list_headers];
else
    idx_inputs = [];
    for i = 1:size(in.considered_inputs,1)
        idx_input = find(strcmp(in.headers,in.considered_inputs{i}));
        if isempty(idx_input)
            errors{end+1} = ['The specified input (' in.considered_inputs{i} ') could no be found in the list of headers:' list_headers];
        end
        idx_inputs = [idx_inputs idx_input];
    end
    in.x = in.data(:,idx_inputs);
end

if isfield(in,'considered_inputs') && isfield(in,'considered_output')
    idx = find(strcmp(in.considered_output,in.considered_inputs));
    if ~isempty(idx)
         errors{end+1} = ['The specified output (' in.considered_output{1} ') is the same as one of the selected inputs'];
    end
end

if ~isfield(in,'consider_temporality')
    in.consider_temporality = false;
end

if in.consider_temporality && sum(strcmp('idx_{sample}',in.considered_inputs))<1
    in.x = horzcat(in.x,(1:size(in.x,1))');
    in.considered_inputs{end+1,1}='idx_{sample}';
end
    
Nrows = size(in.x,1);
Nvars = size(in.x,2);


list_considered_inputs = [];
for i=1:size(in.considered_inputs,1)
    list_considered_inputs = [list_considered_inputs char(10) in.considered_inputs{i}];
end

if isfield(in,'plot_xy')
    if length(in.plot_xy) ~= 2
        errors{end+1} = 'plot_xy should contain exactly two cells (one for the x axis and one for the y axis)';
    else
        in.idx_xy = [];
        for i = 1:2
            idx = find(strcmp(in.considered_inputs,in.plot_xy{i}));
            if isempty(idx)
                errors{end+1} = ['The specified axis input (' in.plot_xy{i} ') could no be found in the list of selected input names:' list_considered_inputs];
            end
            in.idx_xy = [in.idx_xy idx];
        end
    end
end
            

if ~isfield(in,'covfunction') %covariance function
    in.covfunction={'covSEard'};
end

if ~isfield(in,'kfolds') %number of folds for cross-validation
    in.kfolds=0;
end

if ~isfield(in,'perm') %number of permutations for significance
    in.perm=0;
end

if in.kfolds < 0 && in.perm > 0
    errors{end+1} = 'Permutations only work with cross validation. Please set the number of folds to a non-negative value or set the number of permutations to zero';
end

if ~isfield(in,'name') %name of the simulation for plotting and for saving results
    in.name = 'default';
end

if isfield(in,'description')
    if ~iscellstr(in.description) && ~ischar(in.description)
        errors{end+1} = 'The "description" field should be a cell array of strings or a char array';
    end
else
    in.description = {' '};
end

if ~isfield(in,'Ngrid') %number of grid points per dimension
    Ngrid_max = floor(1E5^(1/Nvars));
    in.Ngrid = min(20,Ngrid_max);
end

Ngridpoints = in.Ngrid^Nvars;
if Ngridpoints > 1E6
    warning('You have many grid points (probably due to a high number of considered input variables). This might cause memory issues and can be solved by reducing the Ngrid parameter')
end
    
if ~isfield(in,'hyp')
    in.hyp.cov = zeros(Nvars+1,1);
    in.hyp.lik = 0;
end

if isempty(errors)
    in.checked = true;
end

end

