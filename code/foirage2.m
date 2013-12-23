% map is flat except two bumps


% dataset_microchannel_hari
% 
% Dataset relative to the measurement of heat tranfer coefficients in micochannels heat exchangers. 
% For this analysis, the input set "inputs_raw" is loaded. This corresponds
% to non-processed raw data (i.e. the variables that were modified in the
% experiments. Constants inputs are not taken into account.
%
% symbol      unit        Description
% w           [m]         microchannel width
% w_f         [m]         fin width
% d           [m]         microchannel depth
% P_sat       [bar]       saturation pressure (at exit)
% T_fl_in     [C]         fluid temperature at microchannel inlet
% G           [kg/m²-s]	  mass flux
%
% Source:
% Harirchian, T. and S. V. Garimella, 2008, “Microchannel size effects on local flow boiling heat transfer to a dielectric fluid,” Int. J. Heat Mass Transfer (available online). 
% Cited by: 
% Bertsch, S.S., Groll, E.A., and Garimella, S.V., 2008, "Review and comparative analysis of studies on saturated flow boiling in small channels," Nanoscale and Microscale Thermophysical Engineering, (available online).

clear
load dataset_microchannel_hari.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'w [m]','w_f [m]','d [m]','p_{sat} [bar]','T_{su} [C]','G [kg/m2-s]'};
outputnames = {'h [W/m2-K]'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 4 6];
idx_output=1;

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

in.name = 'foirage';

results = GP_model_data(inputs_raw(:,idx_inputs),y(:,idx_output),in)