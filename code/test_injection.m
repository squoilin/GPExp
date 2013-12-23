% dataset_injection:
% 
% Injection Danfoss compressor, tested with R407c in december/january 2008.
% The points with inproper flow rate measurement (see report) have been removed.
% 
% Input variables: 
% inputs(:,1): p_su [bar]
% inputs(:,2): p_inj [bar]
% inputs(:,3): P_ex [bar]
% inputs(:,4): DELTAT_su [K]
% onputs(:,5): DELTAT_inj [K]
% 
% Outputs: 
% y(:,1): M_dot_su [kg/s]
% y(:,2): M_dot_inj [kg/s]
% y(:,3): W_dot [W]
% y(:,4): T_ex [°C]


clear
load dataset_injection.mat


% Definition of the variables names for plotting and of post-processing:
inputnames = {'p_{su} [bar]','p_{inj} [bar]','P_{ex} [bar]','DELTAT_{su} [K]','DELTAT_{inj} [K]'};
outputnames = {'Mdot_{su} [kg/s]','Mdot_{inj} [kg/s]', 'Wdot [W]', 'T_{ex} [C]'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 2 3];
idx_output=3;

in.Ngrid = 5;

in.perm = 0;
in.kfolds = 5;

in.covhyp=[10*ones(length(idx_inputs),1);0];

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)