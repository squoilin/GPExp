% dataset_injection_ion:
% 
% Injection Emerson compressor, tested by Ion.
% All the points have been taken from the "mesuratori" exell sheet
% 
% Input variables: 
% inputs(:,1): p_su [bar]
% inputs(:,2): p_inj [bar]
% inputs(:,3): P_ex [bar]
% inputs(:,4): T_su [°C]
% onputs(:,5): T_inj [°C]
% 
% Outputs: 
% y(:,1): M_dot_ex [kg/s]
% y(:,2): M_dot_inj [kg/s]
% y(:,3): W_dot [W]
% y(:,4): T_ex [°C]


clear
load dataset_injection_ion.mat


% Definition of the variables names for plotting and of post-processing:
inputnames = {'p_{su} [bar]','p_{inj} [bar]','P_{ex} [bar]','T_{su} [C]','T_{inj} [C]'};
outputnames = {'Mdot_{ex} [kg/s]','Mdot_{inj} [kg/s]', 'Wdot [kW]', 'T_{ex} [C]'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 2 3 4 5];
idx_output=3;

in.Ngrid = 10;

in.perm = 0;
in.kfolds = 5;

in.covhyp=[0*ones(length(idx_inputs),1);0];
% two ok optimums: one at 0 (takes T_su into account), one at 4 (neglects
% Tsu)

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)