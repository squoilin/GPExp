function [y,s2] = GP_prediction(in,out,varargin)
% This function uses a previous GPEXp analysis to perform a prediction of
% the output as a function of the selected inputs for that analysis.
%
% The function inputs are the following:
% in      Inputs to the previous GPExp analysis
% out     Output of the previous GPExp analysis
% inputs  For each input of the anlysis, the value must provided in the
%         form 'inputname=0'. 
%
% As an example, if the analysis consists of two inputs the synthax is:
% GP_prediction(in,out,'inputname1',value1,'inputname2',value2)
% 
% 
% Copyright (c) 2013-2015 Sylvain Quoilin & Jessica Schrouff. 
% All rights reserved.

% Adding the GPML toolbox to the path:

gpmlpath = 'gpml3.5/';   % path to the gpml library


% Check that the path exists:

if exist(gpmlpath,'dir') > 0
    % do nothing
elseif exist(['code/' gpmlpath],'dir') > 0
    gpmlpath = ['code/' gpmlpath];
elseif exist(['../code/' gpmlpath],'dir') > 0
    gpmlpath = ['../code/' gpmlpath];
else
    error('Could not find the GPML library. Check the gpmlpath variable in main_model.m')
end

addpath(gpmlpath)
addpath([gpmlpath 'cov/'])
addpath([gpmlpath 'inf/'])
addpath([gpmlpath 'lik/'])
addpath([gpmlpath 'mean/'])
addpath([gpmlpath 'util/'])


N = length(in.considered_inputs);
N_inputs = (nargin - 2)/2;

InputValues = zeros(N,1);

% Building a string for debugging the code. This might slow it down and
% might be deleted:
list_inputs = [];
for i=1:length(in.considered_inputs)
    list_inputs = [list_inputs char(10) in.considered_inputs{i}];
end

if N ~= N_inputs
    error(['The number of provided inputs (' num2str(N_inputs) ') should be equal to the number of inputs of the analysis (' num2str(N) ')']);
else
    temp = in.considered_inputs;
    for i = 1:N
        strvar = varargin(1+(i-1)*2);
        strvar = strvar{1};
        value = varargin(2+(i-1)*2);
        value = value{1};
        if ischar(strvar) && isnumeric(value)
            idx = find(strcmp(strvar,temp));
            if ~isempty(idx)   
                InputValues(i) = value;
                if iscell(temp) && length(temp)>1
                    temp = temp{[1:idx(1)-1 idx(1)+1:end]};
                else
                    temp = [];
                end
            else     % Could not find the specified input
                error(['Could not find ' strvar ' in the list of considered inputs: ' list_inputs])
            end
        else
            error('The function calling format is not correct. It should look like: GP_prediction(in,out,''inputname1'',value1,''inputname2'',value2)')
        end
    end
end

if ~isempty(temp)
    error(['Value for input ' temp{1} ' does not seem to have been provided']);
end

% Checking the ranges:
for i = 1:N
    if InputValues(i) > max(in.x(:,i)) || InputValues(i) < min(in.x(:,i))
        warning(['The value of variable ' in.considered_inputs{i} '(' num2str(InputValues(i)) ') is outside of the definition domain of the Gaussian Process regression ([' num2str(min(in.x(:,i))) ',' num2str(max(in.x(:,i))) '])'])
    end
end

meanfunc=[];
likfunc = @likGauss;

xx = [InputValues];

[y,s2] = gp(out.hyp, @infExact, meanfunc, in.covfunction, likfunc, in.x, in.y, xx');

end