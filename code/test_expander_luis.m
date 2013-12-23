% dataset_expander_Luis
% 
% Input variables: 
% inputs(:,1): r_p [-]
% inputs(:,2): p_su [bar]
% inputs(:,3): N_rot [rpm]
% inputs(:,4): T_su [°C]
% 
% Outputs: 
% y(:,1): epsilon_s [-]
% y(:,2): M_dot

load dataset_expander_luis.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'r_p [-]','p_{su} [bar]','N_{rot} [rpm]','T_{su} [C]'};
outputnames = {'epsilon_s','Mdot'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 2 3 4];
idx_output=1;

in.xy = [2 3];
in.Ngrid = 10;

in.perm = 0;
in.kfolds = 5;

in.covhyp=[5*ones(length(idx_inputs),1);0];
% 1: ok, pas d'overfit
% 5: ok memem chose que 1
% 10: underfit


in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)