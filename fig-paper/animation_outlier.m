close all
clearvars


addpath('../gpml-matlab-v3.2-2013-01-15');
addpath('../gpml-matlab-v3.2-2013-01-15/util');

type_fitt='poly2';
numparam=3;
in.perm=0;
in.kfolds=-1;

in.axis = [0 1 -0.1 2.5];
in.showboundaries=1;
in.generateMat=0;

% Avoid the minimize function to keep the same hyperparameters throughout
% the whole animation
%in.hyp=log([0.348057609586948;0.824797145620358]);
%in.hyp = [-1.06780572496403;-0.408706489908971];
in.hyp.cov = [-1.2865;   -0.6105];
in.hyp.lik = -3.16279883606355;



%% Animated illustration of Gaussian process


%Vary the number of points
ndp=40;
mae_fit_ndp=zeros(length(ndp),1);
mae_ndp=zeros(length(ndp),1); 
maecv_ndp=zeros(length(ndp),1);
func_ndp=zeros(length(ndp),numparam);
out_ndp=struct();

% Definition of the dataset:
% x=rand(ndp,1);
% y=sin(x*12)+1.2;
% y_std=y+0.1*gpml_randn(0.2,ndp, 1);
% 
%  save inputs.mat x y y_std

load inputs.mat



Nsteps = 25;
factor = (logspace(0,1,Nsteps)-1)/9;

ymax = max(y_std);

yo = y_std;

% Find two points situated at the upper part of the sine wave
gap = abs(y_std-2.2);
[gap,idx] = sort(gap);

for i=1:Nsteps
    
%Defining an outlier:
yo(idx(1:2)) = y_std(idx(1:2))-1.5*(i-1)/Nsteps;
    
    in.name=['sim_outlier_',num2str(i)];
    if ~exist('out')
        out=GP_model_data(x,yo,in);
    else
        out(i)=GP_model_data(x,yo,in);
    end
end

close all