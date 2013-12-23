% dataset_pump_orcnext
% 
% Prediction of the ORCNext pump flow rate. Tests by N. Melotte and
% A.desideri
%

clear
load dataset_pump_orcnext_all.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'\delta_p [pa]','f [Hz]','r_p [-]','p_{su} [pa]','p_{ex} [pa]','T_{su} [C]','p_{ng,su} [pa]','DELTAT_{su} [K]','N_{exp} [rpm]'};
outputnames = {'Mdot'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
% in.covhypp=[0,0];
% in.hypguess=[0;0;0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 2 9];
idx_output=1;

%in.xy = [1 3];

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

in.name = 'pump_orcnext_all';

results = GP_model_data_featselect(inputs(:,idx_inputs),y(:,idx_output),in)