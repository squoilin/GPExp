% dataset_pump_orcnext
% 
% First experimental campaign on the ORCNext setup by Nicolas Melotte
% The pump consumption was measured with a clamp meter, which later turned
% out to be very inaccurate in the low power range, leading to a very poor
% dataset.
% The flow rate was measured by a coriolis flow meter, which should be more
% accurate
%
% Nicolas Melotte, Experimental study and dynamic modeling of a Waste Heat 
% Recovery Organic Rankine Cycle, Master Thesis, University of Liège, 2012

clear
load dataset_pump_orcnext.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'\delta_p* [-]','f* [-]','r_p* [-]','rho* [-]','p* [-]'};
outputnames = {'Mdot','epsilon_s'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 2 5];
idx_output=2;

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

in.name = 'pump_orcnext';

in.covhyp=[0*ones(length(idx_inputs),1);0];
%bad dataset. overfits at -2. At 0, ok

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)