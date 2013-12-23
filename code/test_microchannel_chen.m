% dataset_microchannel_hari
% 
% Dataset relative to the measurement of heat tranfer coefficients in micochannels heat exchangers. 
% For this analysis, the input set "inputs" is loaded. This corresponds
% to non-processed raw data (i.e. the variables that were modified in the
% experiments. Constants inputs are not taken into account.
%
% symbol      unit        Description
% T_fl_in     [C]         fluid temperature at microchannel inlet
% T_surf      [C]         wall temperature
% G           [kg/m²-s]	  mass flux
% qdot        [W/cm²]     heat flux
%
% Only the first 64 points were taken into account, corresponding to the
% first microchannel geometry. There is therefore no influence of any
% geometrical parameter.
%
% Source:
% Chen, T.  and S. V. Garimella, 2007, “Flow boiling heat transfer to a dielectric coolant in a microchannel heat sink,” IEEE Transactions of Components and Packaging Technologies, 30 (1), pp. 24 - 31.
% Cited by: 
% Bertsch, S.S., Groll, E.A., and Garimella, S.V., 2008, "Review and comparative analysis of studies on saturated flow boiling in small channels," Nanoscale and Microscale Thermophysical Engineering, (available online).

% the 3D plot confirms the following statement: 
% The heat transfer coefficients, as well as the critical heat flux (CHF), were found to increase with flow rate

clear
load dataset_microchannel_chen.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'T_{su} [C]','T_{surf} [C]','G [kg/m2-s]','qdot [W/cm2]'};
outputnames = {'h [W/m2-K]'};

% Definition of the kernel to be used
%in.covfunction='covSEiso';
in.covfunction='covSEard';
%in.covhyp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 3 4];
idx_output=1;

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);

%in.xy = [1 3];
in.perm = 0;
in.kfolds = 5;

%in.hypguess=[0*ones(length(idx_inputs),1);0];
%5: underfits T_su, although its weight is not negligible
% start overfitting below -1


in.name = 'microchannel_hari_chen';

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)