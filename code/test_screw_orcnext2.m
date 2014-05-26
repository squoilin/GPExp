% dataset_screw_orcnext
% 
% Nicolas' experimental campaign on the howest screw expander.  All input variables are non-dimensional, except the last one.
% 
% Input variables: 
% inputs(:,1): r_p
% inputs(:,2): p
% inputs(:,3): rpm

% 
% Outputs: 
% y(:,1): epsilon_s [-]

% Load dataset:
% clear all
load dataset_screw_orcnext2.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'rp','p [pa]','rpm'};
outputnames = {'epsilon_s','Mdot'};

% Definition of the kernel to be used
in.covfunction='covSEard';

% Indexes of the inputs/output to be processed:
idx_inputs=[1 3];
idx_output=2;

in.xy = [1 2];
in.Ngrid = 20;

in.perm = 0;
in.kfolds = 5;

%in.hypguess = [1*ones(length(idx_inputs),1);0];

%in.covhyp=[1*ones(length(idx_inputs),1);0];
% < 0.5: overfit sur rp (ondulation dans la map)
% >1 : overfit disparait mais AIC BIC légèrement inférieurs!
% 10: underfit, only rp taken into account

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)