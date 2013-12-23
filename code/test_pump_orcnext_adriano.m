% dataset_pump_orcnext
% 
% First experimental campaign on the ORCNext setup by Nicolas Melotte
% The pump consumption was measured with a clamp meter, which later turned
% out to be very inaccurate in the low power range, leading to a very poor
% dataset.
% The flow rate was measured by a coriolis flow meter, which should be more
% accurate
%

clear
load dataset_pump_orcnext_adriano.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'\delta_p [pa]','f [Hz]','r_p [-]','p_{su} [pa]','p_{ex} [pa]','T_{su} [C]','p_{ng,su} [pa]','DELTAT_{su} [K]'};
outputnames = {'Mdot','epsilon_s','Wdot'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[2 5 7];
idx_output=1;

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

in.name = 'pump_orcnext_adriano';

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)