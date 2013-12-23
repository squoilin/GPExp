% dataset_screw_orcnext
% 
% Nicolas' experimental campaign on the howest screw expander.  All input variables are non-dimensional, except the last one.
% 
% Input variables: 
% inputs(:,1): p_star
% inputs(:,2): rho_star
% inputs(:,3): rpm_star
% inputs(:,4): r_p_star
% inputs(:,5): M_dot (kg/s)
% 
% Outputs: 
% y(:,1): epsilon_s [-]
% y(:,2): phi*Vs

% Load dataset:
% clear all
load dataset_screw_orcnext.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'p*','rho*','rpm*','rp*','M_dot [kg/s]'};

outputnames = {'epsilon_s','phiVs'};

% Definition of the kernel to be used
in.covfunction='covSEard';

% Indexes of the inputs/output to be processed:
idx_inputs=[1 3 4];
idx_output=2;

in.xy = [1 2];
in.Ngrid = 10;

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