% dataset_hermetic2:
% 
% Open-drive expander, tested by S.Quoilin & V. Lemort
% 
% Lemort, V, Declaye, S, & Quoilin, S. (2011). Experimental characterization of a hermetic scroll expander for use in a micro-scale Rankine cycle. Proceedings of the Institution of Mechanical Engineers, Part A, Journal of Power and Energy. 
% 
% Input variables: 
% inputs(:,1): r_p [-]
% inputs(:,2): p_su [bar]
% inputs(:,3): T_su [°C]
% inputs(:,4): X_oil [-]
% 
% Outputs: 
% y(:,1): epsilon_s [-]
% y(:,2): phi [-]
%
% same as hermetic, but 53 data points, among which some were rejected because they were suspected invalid (too large error on the heat balance)



load dataset_hermetic2.mat


% Definition of the variables names for plotting and of post-processing:
inputnames = {'r_p [-]','p_{su} [bar]','T_{su} [C]','X_{oil} [-]'};
outputnames = {'epsilon_s','phi'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 2 3];
idx_output=1;

in.Ngrid = 10;

in.perm = 0;
in.kfolds = 5;

in.covhyp=[10*ones(length(idx_inputs),1);0];



in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

in.name = 'hermetic2';

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)