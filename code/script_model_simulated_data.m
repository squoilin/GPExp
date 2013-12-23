
% y=sin(x/2)+4;
type_fitt='poly2';
numparam=3;
in.perm=100;

%Vary the number of points
ndp=3:3:40;
mae_fit_ndp=zeros(length(ndp),1);
mae_ndp=zeros(length(ndp),1);
maecv_ndp=zeros(length(ndp),1);
func_ndp=zeros(length(ndp),numparam);
out_ndp=struct();
for i=1:length(ndp)
    x=rand(ndp(i),1);
    in.kfolds=3;
    y=-1.5*x.^2+3*x+4;
    in.name=['sim_ndp_',num2str(ndp(i))];
    yo=y+0.1*gpml_randn(0.2,ndp(i), 1);
    if i==1
        tmp=GP_model_data(x,yo,in);
        out_ndp=tmp;
    else
        out_ndp(i)=GP_model_data(x,yo,in);
    end
    mae_ndp(i)=out_ndp(i).train.mae;
    maecv_ndp(i)=out_ndp(i).CV.mae;
    yy=yo-min(yo);
    yy=yy/max(yy);
    xx=x-repmat(mean(x),size(x,1),1);
    xx=xx./repmat(max(abs(x)),size(x,1),1);
    [cfun,gof,outp] = fit(xx,yy,type_fitt);
    mae_fit_ndp(i)=mean(abs(outp.residuals));
    [cfun,gof,outp] = fit(xx,yo,type_fitt);
    aa=struct(cfun);
    bb=aa.coeffValues;
    func_ndp(i,:)=[bb{:}];
%     hold on;
%     plot(cfun,'r')
%     legend('95% boundaries','fitted GP','data','fitted polynom','Location','NorthEastOutside')
    close gcf
end

%covariance function (default squared exponential)
% in.covfunction='covPeriodic';
% in.covhyp=[0 0 0]; %log[ length-scale (ell), period (om), bias (sf)]

ndp=20;
x=rand(ndp,1);
% indout=find(2>x);
% indin=find(7<x);
% x=[x(indout);x(indin)];
ndp=length(x);
in.kfolds=5;

y=-1.5*x.^2+3*x+4;

%Vary the level of noise in data
noise=0.02:0.1:1.12;
mae_fit_noise=zeros(length(noise),1);
mae_noise=zeros(length(noise),1);
maecv_noise=zeros(length(noise),1);
func_noise=zeros(length(noise),numparam);
out_noise=struct([]);
for i=1:length(noise)
    in.name=['sim_noise_',num2str(i)];
    yp=y+noise(i)*gpml_randn(0.2,ndp, 1);
    if i==1
        tmp=GP_model_data(x,yp,in);
        out_noise=tmp;
    else
        out_noise(i)=GP_model_data(x,yp,in);
    end
    mae_noise(i)=out_noise(i).train.mae;
    maecv_noise(i)=out_noise(i).CV.mae;
    yy=yp-min(yp);
    yy=yy/max(yy);
    xx=x-repmat(mean(x),size(x,1),1);
    xx=xx./repmat(max(abs(x)),size(x,1),1);
    [cfun,gof,outp] = fit(xx,yy,type_fitt);
    mae_fit_noise(i)=mean(abs(outp.residuals));
    [cfun] = fit(xx,yp,type_fitt);
    aa=struct(cfun);
    bb=aa.coeffValues;
    func_noise(i,:)=[bb{:}];
%     hold on;
%     plot(cfun,'r')
%     legend('95% boundaries','fitted GP','data','fitted polynom','Location','NorthEastOutside')
    close gcf
end


%Add random outliers

nout=1:10;
mae_fit_outliers=zeros(length(nout),1);
mae_outliers=zeros(length(nout),1);
maecv_outliers=zeros(length(nout),1);
func_outliers=zeros(length(nout),numparam);
indp=randperm(length(x));
out_outliers=struct();
for i=1:length(nout)
    in.name=['sim_outliers_',num2str(nout(i))];
    yo=y+0.1*gpml_randn(0.2,ndp, 1);
    yo(indp(1:nout(i)))=yo(indp(1:nout(i)))+mean(y)+3*rand(1)*std(y);
    if i==1
        tmp=GP_model_data(x,yo,in);
        out_outliers=tmp;
    else
        out_outliers(i)=GP_model_data(x/max(x),yo,in);
    end
    mae_outliers(i)=out_outliers(i).train.mae;
    maecv_outliers(i)=out_outliers(i).CV.mae;
    yy=yo-min(yo);
    yy=yy/max(yy);
    xx=x-repmat(mean(x),size(x,1),1);
    xx=xx./repmat(max(abs(x)),size(x,1),1);
    [cfun,gof,outp] = fit(xx,yy,type_fitt);
    mae_fit_outliers(i)=mean(abs(outp.residuals));
    [cfun] = fit(xx,yo,type_fitt);
    aa=struct(cfun);
    bb=aa.coeffValues;
    func_outliers(i,:)=[bb{:}];
    hold on;
    plot(cfun,'r')
    legend('95% boundaries','fitted GP','data','fitted polynom','Location','NorthEastOutside')
    close gcf
end


%Add a pack of outliers
indo=find(0.25<x<0.4);
% indo=find(2<x<5);
yo=y+0.1*gpml_randn(0.2,ndp, 1);
yo(indo)=mean(y)+4*rand(1)*std(y);
yy=yo-min(yo);
yy=yy/max(yy);
xx=x-repmat(mean(x),size(x,1),1);
xx=xx./repmat(max(abs(x)),size(x,1),1);
in.name=['sim_packoutliers_',num2str(length(indo))];
out_packoutliers=GP_model_data(x,yo,in);
% close gcf
[cfun,gof,outp] = fit(xx,yy,type_fitt);
mae_fit_packoutliers=mean(abs(outp.residuals));
[cfun,gof,outp] = fit(xx,yo,type_fitt);
aa=struct(cfun);
bb=aa.coeffValues;
func_packoutliers=[bb{:}];

%save all results
out_GP=struct('noise',out_noise,'outliers',out_outliers,'packoutliers',...
    out_packoutliers,'ndp',out_ndp);
fit_noise=struct('mae',mae_fit_noise,'func',func_noise);
fit_ndp=struct('mae',mae_fit_ndp,'func',func_ndp);
fit_outliers=struct('mae',mae_fit_outliers,'func',func_outliers);
fit_packoutliers=struct('mae',mae_fit_packoutliers,'func',func_packoutliers);
out_fit=struct('noise',fit_noise,'outliers',fit_outliers,'packoutliers',fit_packoutliers);
save('results_simulated_data.mat','out_GP','out_fit','mae_noise',...
    'maecv_noise','mae_outliers','maecv_outliers','mae_ndp','maecv_ndp')

