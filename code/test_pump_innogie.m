% dataset_expander_Luis
% 
% Input variables: 
% inputs(:,1): r_p [-]
% inputs(:,2): f [Hz]
% inputs(:,3): DELTAT [K]
% inputs(:,4): p_su [bar]
% inputs(:,5): DELTAp [bar]
%
% Outputs: 
% y(:,1): W_dot_sh [W]
% y(:,2): W_dot_el [W]
% y(:,3): epsilon_s [-]
% y(:,4): epsilon_v [-]
% y(:,5): M_dot [kg/s]

load dataset_pump_innogie.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'r_p [-]','f [Hz]','DELTAT [K]','p_{su} [bar]','DELTAp [bar]'};
outputnames = {'W_{sh} [W]','W_{el} [W]','epsilon_s [-]','epsilon_v [-]','M [kg/s]'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Select variable with high volumetric effectiveness to avoid points with
% caviation:
idx_rows = y(:,4) > 0.53;

% Indexes of the inputs/output to be processed:
idx_inputs=[5 2 3];
idx_output=4;

%in.xy = [2 1];
in.Ngrid = 10;

in.perm = 0;
in.kfolds = 5;

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

results = GP_model_data(inputs(idx_rows,idx_inputs),y(idx_rows,idx_output),in)