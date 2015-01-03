%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                %%%%%%%
%%%%   GPExp - A Gaussian Process Experimental Data Analysis Tool   %%%%%%%
%%%%                                                                %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the standard template of the input form of GPExp.

% The tool aims at assessing the
% quality of experimental results with regards to the provided inputs. It
% can help understanding which input variables are the most relevant, what
% the repeatability of the test, which data points are likely to be
% outliers, etc. 

% GPExp is based on the GPML Matlab Library version 3.5 developped by Carl 
% Edward Rasmussen and Christopher K. I. Williams. 

% It runs on Matlab® 7.x and later. 

% The code is released under the FreeBSD License.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, 
%    this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY SYLVAIN QUOILIN AND JESSICA SCHROUFF ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Copyright (c) 2013-2015 Sylvain Quoilin & Jessica Schrouff. 
% All rights reserved.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%   USER INPUTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This section gathers the main data and parameters inputs for GPExp. For
% each user input, a short description is provided with the required input
% format. Some inputs are compulsory, while some others can be left blank
% for the default value.
clear in out

%% 1. Description and origin of the input data. 
% Please provide any relevant information regarding the orgin of the data
% to analyse: citation to publications, type of data, description of the 
% experimental apparatus, etc.
% The 'name' variable is a string used for the file name of the results
% The 'description' variable is a cell array of strings for each line.

in.name = 'absorption';

in.description = strcat({
'Steady-state data relative to a Low-capacity (10 kW) absorption chiller'
'Model: Pink chilli PSC 12'
'Fluid: NH3-H2O'
''
'Data from:'
'Jerko Labus, Modelling of small capacity absorption chillers driven by '
'solar thermal energy or waste heat, PhD Thesis, Universitat Rovira I '
'I Virgili, Tarragona, 2011   '
});


%% 2. Data file
% The data is passed to GPExp to a .mat file contaning two array variables:
% 'inputs' : each colunm corresponds to a different input variable
%            each row is a steady-state measurement
% 'y' :      each column corresponds to an output variable
%            each row is a steady-state measurement
% The file should be placed in the current working directory, or the full
% (relative path must be provided)

load dataset_absorption.mat
in.inputs=inputs;
in.outputs=y;

n = size(inputs,1);
%inputs = inputs([1:10:n],:);
%y = y([1:10:n],:);


%% 3. Variable names
% The variable names (both inputs and outputs) are entered as a cell array
% containing string values. Two variables must be defined: inputnames and 
% outputnames. The number of cells for each variables must correspond to
% the number of columns of the "inputs" and "y" variables of the data file.
% The variables names are further used to refer to the right variable and
% for plotting.

in.inputnames = {
    'T_{chw,in} [C]'
    'T_{cw,in} [C]'
    'T_{hw,in} [C]'
    'Vdot_{chw} [m³/h]'
    'Vdot_{cw} [m³/h]'
    'Vdot_{hw} [m³/h]'
    };

in.outputnames = {
    'Qdot_{eva} [kW]'
    'Qdot_{ac} [kW]'
    'Qdot_{gen} [kW]'
    'COP_{th}'
    'Wdot [kW]'
    'COP_{overall}'
    'EER'
    };

%% 4. Selection of inputs and outputs 
% This defines the subset of the inputs to be considered for the
% processing, as well as the output to be considered.
% Each element of the "considered_inputs" variable should be defined in the
% "inputnames" variable. The "considered_output" is a single cell array,
% whose value must be defined in the "outputnames variables".

in.considered_inputs = {
    'T_{chw,in} [C]'
    'T_{cw,in} [C]'
    'T_{hw,in} [C]'
%    'Vdot_{chw} [m³/h]'
%    'Vdot_{cw} [m³/h]'
%    'Vdot_{hw} [m³/h]'
    };

in.considered_output = {
    'Qdot_{eva} [kW]'
%    'Qdot_{ac} [kW]'
%    'Qdot_{gen} [kW]'
%    'COP_{th}'
%    'Wdot [kW]'
%    'COP_{overall}'
    };


%% 5. Number of folds
% GPExp results are all cross-validation-based: all the performance
% indicators are computed on data points that were not used to train de
% model. The number of data points to leave outside of the training set can
% be user-defined through the "kfolds" parameter. 

% kfolds > 1 : 
% The data is randomly split into k folds. The model is trained using k-1
% folds and the error indicators are computed using the last fold. This is
% done iteratively for each fold.

% kfolds == 1 :
% Split in half: the data samples are split into two folds. The model is
% trained on one of them and its accuracy is checked using the second one.

% kfolds == 0 :
% Leave-one-out: the model is trained on all samples minus one (used for
% the cross-validation). This is equivalent to kfolds == Nsamples

% kfolds == -1 : 
% Ignore kfold cross-validation (not recommended). All data samples are 
% used to train the model.

% Default value: 0

in.kfolds =-1;


%% 6. Number of permutations
% To assess the significance of the regression performance compared to the 
% null hypothesis, permutations are used: 
% the input rows are randomly permuted to obtain a “baseline” model 
% performance (i.e. the performance of the regression for a totally random
% and non-correlated set of inputs/output). The “true” model performance is 
% then compared to the baseline level.
% Typically, a model performance is significant (i.e. there is a significant 
% relationship between inputs and output) if the performance obtained by 
% chance does not exceed or equal the true model performance more than 
% 5% of the time (p < 0.05).

% It should be noted that permutations significantly increase the
% computational time. If the proof of a relationship between the output and
% the inputs is not critical for the user, permutations should not be used.

% The number of permutation can be set to zero to disable this option
% If used, it should be set to a number sufficient to ensure statistical
% relevance.

% NB: the permutation method uses cross-validation, kfolds must therefore
% be set to a positive value to activate it.

% Default value: 0

in.perm = 0;


%% 7. 3D plotting options
% Ngrid allows setting the number of grid points in one direction (the
% final number is Ngrid²).
% If not specified, Ngrid is set by default to 20
%
% plot_inputs is a 2-cells array with the names of the inputs variables to
% be plotted. 
% If not specificed, the two most relevant variables (i.e. the ones with
% the lowest lengthscales) are plotted

%in.Ngrid = 15;

% in.plot_xy = {
%     'T_{chw,in} [C]'
%     'T_{cw,in} [C]'
%     'T_{hw,in} [C]'
%     'Vdot_{chw} [m³/h]'
%     'Vdot_{cw} [m³/h]'
%     'Vdot_{hw} [m³/h]'
%     };


%% 8. Definition of the kernel
% These parameters should be modified by experienced users only. 
% 
% Definition of the covariance function. Default is square exponential with
% automatic relevant detection (SEard)

%in.covfunction='covSEard';

% Definition of the hyperparameters. This avoids the use of the
% optimization algorithm and can be activated if the latter has been run at
% least once before. The hyperparameters resulting of the optimization can
% be obtained in the results variable ('results.hyp').

% In general three sets of hyperparameters must be defined, corresponding
% to the mean function, the covariance function, and the likelihood
% function. In GPExp, the mean function is set to 0 => there is ne need for
% hyperparameters. If an ARD kernel is used (most recommended case), the
% hyperparameter of the covariance function are the logarithmic lengthscales 
% for each dimension and the log standard deviation. The likelihood function of
% the prior is specified to be gaussian in GPExp. It corresponds to the
% expected noise level and is defined by the log of its standard deviation.

% If undefined, the parameters are obtained by minimizing the marginal
% likelihood. If defined, the hyperparameters are imposed for the "train"
% regression (i.e. with all data samples) and used as guesses for the
% cross-validations.

%in.hyp.cov= log([3.7088; 2.9898; 1.1500; 0.5153]);     % 3 lengthscales and one standard deviation 
%in.hyp.lik=log(0.0113);                                % Standard deviation of the likelihood

%in = inputs_sanity_check(in);
results = main_model(in)

