% dataset_innogie:
% 



clear
load dataset_innogie_oli.mat


% Definition of the variables names for plotting and of post-processing:
inputnames = {'Mdot_r [g/s]','Mdot_w [l/s]','dp [bar]','p_{su} [bar]','T_{su} [C]','T_{ex} [C]'};
outputnames = {'pinch_{cd} [K]'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 2];
idx_output=1;

in.Ngrid = 20;

%in.xy = [1 3];


in.perm = 0;
in.kfolds = 5;

in.covhyp=[0*ones(length(idx_inputs),1);0];
% two ok optimums: one at 0 (takes T_su into account), one at 4 (neglects
% Tsu)

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)