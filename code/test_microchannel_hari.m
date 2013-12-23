% dataset_microchannel_hari
% 
% Dataset relative to the measurement of heat tranfer coefficients in micochannels heat exchangers. 
%
% symbol      unit        Description
% fluid       -           fluid used for the experiment
% D_h         [m]         hydraulic diameter
% aspect      -           aspect ratio: channel height/width
% L_inl       [m]         distance of testpoint from inlet of the microchannel
% P_sat       [bar]       saturation pressure (at exit)
% T_fl_in     [C]         fluid temperature at microchannel inlet
% T_surf      [C]         surface (wall) temperature of the microchannels
% G           [kg/m²-s]	  mass flux
% q_flux      [W/m²]      heat flux (with respect to heated wetted area)
% x            -          thermodynamic quality
% h           [W/m²-K]	  heat transfer coefficient
%
% Source:
% Harirchian, T. and S. V. Garimella, 2008, “Microchannel size effects on local flow boiling heat transfer to a dielectric fluid,�? Int. J. Heat Mass Transfer (available online). 
% Cited by: 
% Bertsch, S.S., Groll, E.A., and Garimella, S.V., 2008, "Review and comparative analysis of studies on saturated flow boiling in small channels," Nanoscale and Microscale Thermophysical Engineering, (available online).

clear all
load dataset_microchannel_hari.mat

% Definition of the variables names for plotting and of post-processing:
inputnames = {'D_h [m]','apect [-]','p_{sat} [bar]','T_{su} [C]','T_{surf} [C]', 'G [kg/m2-s]', 'q_{flux} [W/cm2]','x [-]'};
outputnames = {'h [W/m2-K]'};

% Definition of the kernel to be used
%in.covfunction='covSEard';
%in.covhypp=[0,0];

% Indexes of the inputs/output to be processed:
idx_inputs=[1 8];
idx_output=1;

in.inputnames = inputnames(idx_inputs);
in.outputname = outputnames(idx_output);
in.ngrid=5;
in.name = 'microchannel_hari';

results = GP_model_data(inputs(:,idx_inputs),y(:,idx_output),in)