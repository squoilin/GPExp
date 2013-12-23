% dataset_opendrive: 
% Open-drive expander, tested by S.Declaye and L. Grégoire.
% 
% Sébastien Declaye, Sylvain Quoilin, Ludovic Guillaume and Vincent Lemort, Experimental study on an open-drive scroll expander integrated into an ORC system working with R245fa, Energy, 2013
% 
% Input variables: 
% inputs(:,1): r_p [-]
% inputs(:,2): N_rot [rpm]
% inputs(:,3): p_su [bar]
% 
% Outputs: 
% y1: epsilon_s [-]
% y2: phi [-]

clear
load dataset_opendrive.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'r_p [-]','N_{rot} [rpm]','p_{su} [bar]'};
outputnames = {'epsilon_s','phi'};

% Definition of the kernel to be used
in.covfunction='covSEard';

% Indexes of the inputs/output to be processed:
idx_inputs=[1 2 3];
idx_output=1;

in.perm = 0;
in.kfolds = -1;

% in.hyp.cov = [-2;-0.251562787897010];
% in.hyp.lik = -2.992;
% in.hyp.cov = [-1;-4];
% in.hyp.lik = -2.992;

%in.hypguess = [1*ones(length(idx_inputs),1);0];

in.xy = [1 2];
in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)