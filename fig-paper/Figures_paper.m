clearvars
close all

addpath('../gpml-matlab-v3.2-2013-01-15');
addpath('../gpml-matlab-v3.2-2013-01-15/util');

type_fitt='poly2';
numparam=3;
in.perm=0;
in.generateMat=0;
in.showboundaries=1;
in.kfolds=1;

%% 1-D plot with sine function and one outlier


%Vary the number of points
ndp=40:40;
mae_ndp=zeros(length(ndp),1); 
maecv_ndp=zeros(length(ndp),1);
func_ndp=zeros(length(ndp),numparam);
out_ndp=struct();

% x=rand(ndp,1);
% y=sin(x*12)+4;
% y_std=y+0.5*gpml_randn(0.2,ndp, 1);
% 
% save inputs_single.mat

load inputs_single.mat

y_std(20) = y_std(20)+1;

for i=1:length(ndp)
    in.name=['sim_ndp_',num2str(ndp(i))];
    % See the impact of a reduced number of data points:
    yo=y_std;
%     idx = 1:2:ndp;
%     yo = yo(idx);
%     x = x(idx);
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
    [cfun,gof,outp] = fit(x,yo,type_fitt);
    aa=struct(cfun);
    bb=aa.coeffValues;
    func_ndp(i,:)=[bb{:}];
    
    hold on;
    plot(cfun,'r')
    hold off;
    % legend('95% boundaries','fitted GP','data','fitted polynom','Location','NorthEastOutside')
    %close gcf
end




%% 1-D plot with polynomial function and one outlier


%Vary the number of points
ndp=40:40;
mae_fit_ndp=zeros(length(ndp),1);
mae_ndp=zeros(length(ndp),1);
maecv_ndp=zeros(length(ndp),1);
func_ndp=zeros(length(ndp),numparam);
out_ndp=struct();
for i=1:length(ndp)
    x=rand(ndp(i),1);
    in.kfolds=1;
    y=-1.5*x.^2+3*x+4;
    %y=sin(x*12)+4;
    %Defining an outlier:
    y(20) = y(20)+1;
    in.name=['sim_ndp_',num2str(ndp(i))];
    yo=y+0.25*gpml_randn(0.2,ndp(i), 1);
    if i==1
        tmp=GP_model_data(x,yo,in);
        out_ndp=tmp;
    else
        out_ndp(i)=GP_model_data(x,yo,in);
    end
    mae_ndp(i)=out_ndp(i).train.mae;
    maecv_ndp(i)=out_ndp(i).CV.mae;
    [cfun,gof,outp] = fit(x,yo,type_fitt);
    aa=struct(cfun);
    bb=aa.coeffValues;
    func_ndp(i,:)=[bb{:}];
    
     hold on;
     plot(cfun,'r')
     hold off;
    % legend('95% boundaries','fitted GP','data','fitted polynom','Location','NorthEastOutside')
    %close gcf
end



for i=1:length(ndp)
    x1=2*(rand(ndp(i),1)-0.5);
    x2=2*(rand(ndp(i),1)-0.5);
    in.kfolds=1;
    y=sqrt(x1.^2+x2.^2)/1.414 ;
    %Defining an outlier:
    y(20) = y(20)+0.5;
    in.name=['sim_ndp_',num2str(ndp(i))];
    yo=y+0.05*gpml_randn(0.2,ndp(i), 1);
    if i==1
        tmp=GP_model_data([x1 x2],yo,in);
        out_ndp=tmp;
    else
        out_ndp(i)=GP_model_data([x1 x2],yo,in);
    end
    

    % hold on;
    % plot(cfun,'r')
    % legend('95% boundaries','fitted GP','data','fitted polynom','Location','NorthEastOutside')
    %close gcf
end