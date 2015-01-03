close all
clearvars


addpath('../gpml-matlab-v3.2-2013-01-15');
addpath('../gpml-matlab-v3.2-2013-01-15/util');

type_fitt='poly2';
numparam=3;
in.perm=0;
in.kfolds=0;

in.axis = [0 1 -1 3];
in.showboundaries=1;
in.generateMat=0;

% Avoid the minimize function to keep the same hyperparameters throughout
% the whole animation
%in.hyp=log([0.348057609586948;0.824797145620358]);
%in.hyp = [-1.06780572496403;-0.408706489908971];
in.hyp.cov = [-1.2865;   -0.6105];
in.hyp.lik = -3.16279883606355;



%% Animated illustration of Gaussian process


% Definition of the dataset:
% x=rand(ndp,1);
% y=sin(x*12)+1.2;
% y_std=y+0.1*gpml_randn(0.2,ndp, 1);
% 
%  save inputs.mat x y y_std

load inputs.mat

x_in = x(1:15,:);
y_in = y(1:15,:);
y_std_in = y_std(1:15,:);

%Defining the noise level:

delta = (y_in - y_std_in);
noise = 5;                  % Defining the noise level
yo_in = y_std_in + noise*delta;

%Vary the number of points
ndp=size(x_in,1);

Nsteps = 25;

ymax = max(y_std);

for i=1:ndp             % iterating on the leave out points
    
x = x_in([1:i-1 i+1:end],:);
y = y_in([1:i-1 i+1:end],:);
yo = yo_in([1:i-1 i+1:end],:);
%yo = y_std_in;
%x = x_in;
%y = y_in;

    in.name=['sim_CVgoodfit_',num2str(i)];
    in.x_leavout=x_in(i);
    in.y_leavout=yo_in(i);
    if ~exist('out')
        out=GP_model_data(x,yo,in);
    else
        out(i)=GP_model_data(x,yo,in);
    end    
end





%close all