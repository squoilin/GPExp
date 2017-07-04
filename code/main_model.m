function [out]=main_model(in)

% Gaussian Processes modelling of the data x (n_samples * n_variables) to
% the regression target y (n_samples). All values have to be continuous.
% When 1 or 2 variables only, plots are generated. Other inputs can be
% specified, such as the:
%      - in.covfunction: type of covariance function/kernel to use. By
%      defaults, 'covSEard'. See usageCov.m in the GPML toolbox for other
%      options and specification of hyperparameters.
%      - in.covhyp: specifies the initialization of the covariance
%      hyperparameter.
%      - in.kfolds: number of folds (i.e. dataset divisions) to perform
%      cross-validation.
%      - in.perm: number of permutations for non-parametric testing of the
%      performance of the model.
%      - in.name: name or number associated with the dataset
% save_fig and view_fig are two plot_dims, set to 1 by default to view the
% figures if 1D or 2D input variables.
% The data will then be modelled, using the whole dataset and then using
% cross-validation. Different parameters assessing the quality of the data
% and model will then be computed and output:
%      - out.AIC= Akaike Information Criterion (corrected)
%      - out.BIC= Bayesian Information Criterion
%      - out.hypcov = optimized hyperparameters for the covariance function
%      - out.train = structure comprising the mean absolute error (MAE) and
%      the correlation between the predictions and the targets, as well as
%      their p-values, when considering the whole dataset.
%      - out.CV = same structure as for train but when considering k-folds
%      cross-validation.
%      - out.outliers = vector comprising the indexes of the outliers, i.e.
%      the datapoints outside the 95% confidence intervals defined by the
%      model.
%--------------------------------------------------------------------------
% Written by J. Schrouff and S. Quoilin, University of Liège, 2013-2014
% Using the GPML toolbox written by C.E. Rasmussen and C.K. Williams (2006)

% Adding the GPML toolbox to the path:

gpmlpath = 'gpml3.5/';   % path to the gpml library


% Check that the path exists:

if exist(gpmlpath,'dir') > 0
    % do nothing
elseif exist(['code/' gpmlpath],'dir') > 0
    gpmlpath = ['code/' gpmlpath];
elseif exist(['../code/' gpmlpath],'dir') > 0
    gpmlpath = ['../code/' gpmlpath];
else
    error('Could not find the GPML library. Check the gpmlpath variable in main_model.m')
end

addpath(gpmlpath)
addpath([gpmlpath 'cov/'])
addpath([gpmlpath 'inf/'])
addpath([gpmlpath 'lik/'])
addpath([gpmlpath 'mean/'])
addpath([gpmlpath 'util/'])

%Sanity check and initialize
%----------------------------
tic

if ~isfield(in,'checked')
    in = inputs_sanity_check(in);
end

x = in.x;
y = in.y;
covfunction=in.covfunction;
k=in.kfolds;
nperm=in.perm;
name = in.name;
Ngrid = in.Ngrid;
hyp = in.hyp;

meanfunc=[];

%initialize the likelihood (Gaussian, with exact inference)
likfunc = @likGauss;
%for weights (linear kernel)
hyp_lin=hyp;
hyp_lin.cov=0;
%number of function evaluation
%nfe=-100;
nfe=-200;


%Prepare inputs
%----------------

%Normalize the variables in order to obtain values ranging from 0 to 1 for
%each of them.
ns=size(x,1);
nv=size(x,2);
% Save the original value of x:
xx = x;
moyx = repmat(mean(x),ns,1);
x=x-moyx;
maxx = repmat(max(abs(x)),ns,1);
x=x./maxx;

%normalize the output variable
yy=y;
moyy = min(y);
y = y-moyy;
maxy = max(y);
y=y./maxy;

%Are there dependent variables? i.e. is the matrix full rank?
s = svd(x);
ss=s/max(s);
% tol = size(x,2)*eps(max(s)); %matlab form, very conservative
%r = sum(s > tol);
tol=0.001;
r=length(find(ss> tol));
% if r<size(x,2)
%     warning('MATLAB:Rankdeff','Rank defficient input matrix: dependent variables in input')
%     ind=find(ss<=tol);
%     for i=1:length(n)
%         sprintf('Variables number %d \n',ind(i))
%     end
% end

bi=min(x);
bs=max(x);
if nv==1  % create 1D grid
    z=linspace(bi,bs,1000)';
    plot_dim=1;
else
    if nv==2 %create 2D grid, rectangular
        plot_dim=2;
        %[z1,z2]=meshgrid(bi(1):0.1:bs(1),bi(2):0.1:bs(2));
        [z1,z2]=meshgrid(linspace(bi(1),bs(1),Ngrid),linspace(bi(2),bs(2),Ngrid));
        z=[z1(:),z2(:)];
    elseif nv==3 %create 3D grid, parallelepiped
        plot_dim=3;
        [z1,z2,z3]=ndgrid(linspace(bi(1),bs(1),Ngrid),linspace(bi(2),bs(2),Ngrid),linspace(bi(3),bs(3),Ngrid));
        z=[z1(:),z2(:),z3(:)];
    else        %create nd grid, hypercubic
        out=[]; out2=[];
        for i=1:nv
            out=[out,[',z',num2str(i)]];
            out2 = [out2, [',linspace(bi(',num2str(i),'),bs(',num2str(i),'),Ngrid)']];
        end
        eval([ '[' out(2:end) '] =ndgrid(' out2(2:end) ');'])
        z=[];
        for i=1:nv
            eval([ 'z = [z, z' num2str(i) '(:)];'])            
        end
        plot_dim=3;
    end 
end


%GP modelling
%---------------
pcvmae=0;
out=struct();
for p=0:nperm %for each permutation of the targets + the "true" targets
    if p>0
        ind=randperm(length(y));
        yp=y(ind);
    else
        yp=y;
    end
    
    %Train on the whole dataset
    %---------------------------  
    
    if p==0
        %first, optimize the hyperparameters
        hyp = in.hyp;
        hyp = minimize(hyp, @gp, nfe, @infExact, meanfunc, covfunction, likfunc, x, yp);
        %test on the dataset
        [m, s2] = gp(hyp, @infExact, meanfunc, covfunction, likfunc, x, yp, x);
        out.train.y_pred = m*maxy + moyy;
        out.train.mae=mean(abs(m-yp));
        [out.train.rsquare out.train.rmse] = rsquare(out.train.y_pred,yy);
        out.hypcov = exp(hyp.cov);
        out.hyp = hyp;
        out.outliers=(y-m)./sqrt(s2);  %detect outliers
        %compute model evidence, AIC and BIC
        nlml = gp(hyp, @infExact, meanfunc, covfunction, likfunc, x, yp);
        out.AIC=2*nlml+2*(length(hyp.cov));             % Akaike information criterion
        out.BIC=2*nlml+(length(hyp.cov))*log(ns);       % Bayesian information criterion
        disp('Training completed')
        
        %Linear kernel: determine weights for each variable
        hyp_lin = minimize(hyp_lin, @gp, nfe, @infExact, meanfunc, {@covLINone}, likfunc, x, yp);
        [d1,d2,d3,d4,d5,post] = gp(hyp_lin, @infExact, meanfunc, {@covLINone}, likfunc, x, yp, x);
        alpha=post.alpha;
        weights=x'*alpha;
        out.weights=abs(weights/norm(weights));
        disp('computation of weights completed (linear kernel)')
        
        if plot_dim ~= 0
            %compute the predicted values and their variances over the whole hypercube:
            [m s2] = gp(hyp, @infExact, meanfunc, covfunction, likfunc, x, yp, z);
            m = m*maxy + moyy;
            s2 = s2*(maxy^2);
            % Restore the original value of the inputs by unscaling z:
            zz = zeros(size(z,1),nv);
            for i = 1:nv
                zz(:,i) = z(:,i)*maxx(1,i)+moyx(1,i);
            end
            out.ndgrid.y_gp=m;
            out.ndgrid.s2=s2;
            out.ndgrid.z=zz;
        disp('Computation of the ndgrid completed')    
        end
    end
    
    % Cross-validation
    %-------------------
    
    %compute number of data points per folds
    if k>1
        if k>ns
            error('GP_model_data:largekfolds',['Number of folds too large',...
                'compared to number of samples, please correct'])
        end
        nsf=floor(ns/k);
        mns=mod(ns,k);
        dk=nsf*ones(1,k);           % Vector with the number of samples for each fold
        dk(end)=dk(end)+mns;
        inds=1;
        sk=[];
        for ii=1:length(dk)
            sk=[sk,inds*ones(1,dk(ii))];
            inds=inds+1;
        end
    elseif k==0 %Leave-One-Sample-Out
        sk=1:ns;
    elseif k==1 %split in half
        ins=floor(ns/2);
        sk=[ones(1,ins),2*ones(1,ns-ins)];
    elseif k==-1
        disp('Ignoring kfold cross-validation')
    else
        error('GP_model_data:negativekfolds',['Number of folds negative',...
            ', please correct'])
    end
    
    if k>=0
        pt=zeros(ns,1);
        for i=1:max(sk)
            xtr=x(sk~=i,:);     % All inputs that don't belong to fold i
            xte=x(sk==i,:);     % All inputs that belong to fold i
            ytr=yp(sk~=i,:);    % All outputs that don't belong to fold i
            yte=yp(sk==i,:);    % All outputs that belong to fold i
            %optimize hyperparameters
            hyp = in.hyp;
            %hyp = out.hyp;             % Take the results of the train as guess for the CV optimization (!! can be assimilated to double-dipping)
            hyp = minimize(hyp, @gp,nfe, @infExact, meanfunc, covfunction, likfunc, xtr, ytr);
            %test model on predictions
            [m] = gp(hyp, @infExact, meanfunc, covfunction, likfunc, xtr, ytr, xte,yte);
            pt(sk==i)=m;
        end
        %compute mae
        tempm=mean(abs(pt-yp));
        if p==0                     % If we are not permuted, compute metrics for the goodness of fit
            out.CV.mae=tempm;
            out.CV.y_pred=pt*maxy + moyy;
            [out.CV.rsquare out.CV.rmse] = rsquare(out.CV.y_pred,yy);
        else                        % If we are permuted, calculate the average mae for each permutation
            if tempm<=out.CV.mae
                pcvmae=pcvmae+1;
            end
        end
        disp('Computation of the cross validation completed')
    end
end

out.ss = ss;

if k>=0
    out.CV.pmae=(pcvmae+1)/(nperm+1);
end

toc

