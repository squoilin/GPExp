% Dataset_expander_emeline
% 
% Expander tested by Emeline for her master thesis.
% 
% Input variables: 
% inputs(:,1): r_p [-]
% inputs(:,2): p_su [bar]
% inputs(:,3): N_rot [rpm]
% inputs(:,4): T_su [°C]
% 
% Outputs: 
% y(:,1): epsilon_s [-]
% y(:,2): phi [-]
% y(:,3): T_ex [°C]
% y(:,4): W_dot [W]

clear
load dataset_expander_emeline.mat



% Definition of the variables names for plotting and of post-processing:
inputnames = {'r_p [-]','p_{su} [bar]','N_{rot} [rpm]','T_{su} [C]'};
outputnames = {'epsilon_s [-]','phi [-]', 'T_ex [°C]', 'W_dot [W]'};

% Indexes of the inputs/output to be processed:
idx_inputs=[1 2 3];
idx_output=1;

% Definition of the kernel to be used
in.covfunction='covSEard';
in.covhyp=[4*ones(length(idx_inputs),1);0];
% right range for the start values: -2 to 3


%in.hypguess = [1*ones(length(idx_inputs),1);0];
%in.Ngrid = 5;

in.perm = 0;
in.kfolds = -1;

in.xy = [1 3];

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);


results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)

