function [in] = inputs_sanity_check(in)
% This function performs a sanity check on the input structure variable of
% GPExp. It returns error messages to assist the user in the model
% parametrization process.
% 
% University of Liège, December 2014
% 
% Copyright (c) 2013-2015 Sylvain Quoilin & Jessica Schrouff. 
% All rights reserved.

if nargin<1
    error('Please enter a valid input structure')
end

if ~isfield(in,'inputs') || ~isfield(in,'outputs')
    error('The "in" structure should at least contain the "inputs" and "outputs" fields')
end

% Checking inputs and input names:
if size(in.inputs,1) ~= size(in.outputs,1)
    error('The data matrix (inputs) and the targets (outputs) should have the same number of samples')
end

if ~isfield(in,'inputnames') %name of the input variables
    in.inputnames=cell(size(in.inputs,2),1);
    for i=1:size(in.inputs,2)
        in.inputnames{i} = ['Input ',num2str(i)];
    end
else
    if size(in.inputnames,1) ~= size(in.inputs,2)
        error(['The number of specified input names (' num2str(size(in.inputnames,1)) ') does not match with the number of columns of the input array (' num2str(size(in.inputs,2)) ')'])
    end
end

if ~isfield(in,'outputnames') %name of the output variables
    in.outputnames=cell(size(in.outputs,2),1);
    for i=1:size(in.outputs,2)
        in.outputnames{i} = ['output ',num2str(i)];
    end
else
    if size(in.outputnames,1) ~= size(in.outputs,2)
        error(['The number of specified output names (' num2str(size(in.outputnames,1)) ') does not match with the number of columns of the output array (' num2str(size(in.outputs,2)) ')'])
    end
end

list_inputs = [];
for i=1:size(in.inputs,2)
    list_inputs = [list_inputs char(10) in.inputnames{i}];
end

list_outputs = [];
for i=1:size(in.outputs,2)
    list_outputs = [list_outputs char(10) in.outputnames{i}];
end


% Checking output and output name:
if ~isfield(in,'considered_output') && size(in.outputs,2) ~= 1
    error(['Please specify the output to be considered in the "considered_output" variable. Possible values are: ' list_outputs])
elseif ~isfield(in,'considered_output') && size(in.output,2) == 1
    in.considered_output = in.outputnames{1};
else
    if iscell(in.considered_output) && size(in.considered_output,1)>1
        error('Predicting one target (output) at a time, please execute multiple times for multiple targets')
    end
    idx_output = find(strcmp(in.outputnames,in.considered_output));
    if isempty(idx_output)
        error(['The specified output (' in.considered_output ') could no be found in the list of output names:' list_outputs])
    end
    in.y = in.outputs(:,idx_output);
end
        
% Checking inputs and inputs names:
if ~isfield(in,'considered_inputs') && size(in.inputs,2) ~= 1
    error(['Please specify the inputs to be considered in the "considered_inputs" variable. Possible values are: ' list_inputs])
elseif ~isfield(in,'considered_inputs') && size(in.input,2) == 1
    in.considered_inputs = in.inputnames{1};
else
    idx_inputs = [];
    for i = 1:size(in.considered_inputs,1)
        idx_input = find(strcmp(in.inputnames,in.considered_inputs{i}));
        if isempty(idx_input)
            error(['The specified input (' in.considered_inputs{i} ') could no be found in the list of input names:' list_inputs])
        end
        idx_inputs = [idx_inputs idx_input];
    end
    in.x = in.inputs(:,idx_inputs);
end
    
Nrows = size(in.x,1);
Nvars = size(in.x,2);


list_considered_inputs = [];
for i=1:size(in.considered_inputs,1)
    list_considered_inputs = [list_considered_inputs char(10) in.considered_inputs{i}];
end

if isfield(in,'plot_xy')
    if size(in.plot_xy,1) ~= 2
        error('plot_xy should contain exactly two cells (one for the x axis and one for the y axis)')
    else
        in.idx_xy = [];
        for i = 1:2
            idx = find(strcmp(in.considered_inputs,in.plot_xy{i}));
            if isempty(idx)
                error(['The specified axis input (' in.plot_xy{i} ') could no be found in the list of selected input names:' list_considered_inputs])
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
    error('Permutations only work with cross validation. Please set the number of folds to a non-negative value or set the number of permutations to zero')
end

if ~isfield(in,'name') %name of the simulation for plotting and for saving results
    in.name = 'default';
end

if isfield(in,'description')
    if ~iscellstr(in.description)
        error('The "description" field should be a cell array of strings')
    end
else
    in.description = {' '};
end

if ~isfield(in,'Ngrid') %number of grid points per dimension
    Ngrid_max = floor(1E5^(1/Nvars));
    in.Ngrid = max(20,Ngrid_max);
end

Ngridpoints = in.Ngrid^Nvars;
if Ngridpoints > 1E6
    warning('You have many grid points (probably due to a high number of considered input variables). This might cause memory issues and can be solved by reducing the Ngrid parameter')
end
    
if ~isfield(in,'hyp')
    in.hyp.cov = zeros(Nvars+1,1);
    in.hyp.lik = 0;
end

in.checked = true;

end

