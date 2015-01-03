% dataset_screw_orcnext
% 
% Adriano's experimental campaign on the howest screw expander.  All input variables are non-dimensional, except the last one.
% 


% Load dataset:
% clear all
load dataset_screw_orcnext3.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'p_ex','p_su [pa]','Tsu','rpm','rp'};
outputnames = {'epsilon_s','Mdot','Wdot','FFVs','Tex'};

% Definition of the kernel to be used
in.covfunction='covSEard';

% Indexes of the inputs/output to be processed:
idx_inputs=[1 5];
idx_output=1;

%in.xy = [1 2];
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