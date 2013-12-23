function [out]=GP_model_data_featselect(x,y,in)

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
% save_fig and view_fig are two flags, set to 1 by default to view the
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
% Written by J. Schrouff, Cyclotron Research Centre, 31/01/2013
% Using the GPML toolbox written by C.E. Rasmussen and C.K. Williams (2006)

% todo: save the hypercube and its vectors
% optimize the minimization

%recursive feature addition/elimination

%Sanity check and initialize
%----------------------------
if nargin<2
    error('Please enter the data matrix and the regression targets as two inputs')
end

if size(x,1) ~= size(y,1)
    error('The data matrix and the targets should have the same number of samples')
end

if size(y,2)>1
    error('Predicting one target at a time, please execute multiple times for multiple targets')
end

if nargin==3 || nargin == 2
    if isfield(in,'covfunction') %covariance function
        covfunction={in.covfunction};
    else
        covfunction={'covSEard'};
    end
    if isfield(in,'kfolds') %number of folds for cross-validation
        k=in.kfolds;
    else
        k=5;
    end
    if isfield(in,'perm') %number of permutations for significance
        nperm=in.perm;
    else
        nperm=10;
    end
    if isfield(in,'name') %name of the simulation for plotting and for saving results
        name = in.name;
    else
        name = 'default';
    end   
    if isfield(in,'inputnames') %name of the input variables
        inputnames = in.inputnames;
    else
        inputnames=cell(size(x,2),1);
        for i=1:size(x,2)
            inputnames{i} = ['Input ',num2str(i)];
        end
    end   
    if isfield(in,'outputname') %name of the output variable
        outputname = in.outputname;
    else
        outputname = {'Output'};
    end     
    if isfield(in,'Ngrid') %number of grid points per dimension
        Ngrid = in.Ngrid;
    else
        Ngrid = 20;
    end         
end

% meanfunc = @meanConst;
% hyp.mean = 0;
meanfunc=[];

%initialize the likelihood (Gaussian, with exact inference)
likfunc = @likGauss;
hyp.lik=0;
%for weights
hyp2=hyp;
hyp2.cov=0;
%number of function evaluation
% nfe=-100;
nfe=-200;

%tolerance for feature selection
tolfeat=0.005;

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
if r<size(x,2)
    warning('MATLAB:Rankdeff','Rank defficient input matrix: dependent variables in input')
    ind=find(ss<=tol);
    for i=1:length(n)
        sprintf('Variables number %d \n',ind(i))
    end
end

%compute number of data points per folds
if k>1
    if k>ns
        error('GP_model_data:largekfolds',['Number of folds too large',...
        'compared to number of samples, please correct'])
    end
    nsf=floor(ns/k);
    mns=mod(ns,k);
    dk=nsf*ones(1,k);
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


%Main loop: nested cross-validation
%--------------------------------------------------------------------------
pt=zeros(ns,1);
yp=y;
pval=zeros(k,nperm);
nfsel=zeros(k,1);
for i=1:max(sk)
    xtr=x(sk~=i,:);
    xte=x(sk==i,:);
    ytr=yp(sk~=i,:);
    
    %feature selection step
    %----------------------------------------------------------------------
    
    %compute the number of data points per subfold
    if k>1
        nsf=floor(size(xtr,1)/k);
        mns=mod(size(xtr,1),k);
        dk=nsf*ones(1,k);
        dk(end)=dk(end)+mns;
        inds=1;
        sksub=[];
        for ii=1:length(dk)
            sksub=[sksub,inds*ones(1,dk(ii))];
            inds=inds+1;
        end
    elseif k==0 %Leave-One-Sample-Out
        sksub=1:ns;
    elseif k==1 %split in half
        ins=floor(ns/2);
        sksub=[ones(1,ins),2*ones(1,ns-ins)];
    end
    
    
    %strategy for feature selection: Recursive Feature Addition (my paper)
    % based on Recursive Feature Elimination (RFE, Kriegeskorte, 2006)
    tempms=zeros(size(x,2),1);
    features=1:size(x,2);
    pts=zeros(size(xtr,1),1);
    for nfeat=features %start with all features and decrease number
        for j=1:max(sksub)
            xtrs=xtr(sksub~=j,:);
            xtes=xtr(sksub==j,:);
            ytrs=ytr(sksub~=j,:);
            ytes=ytr(sksub==j,:);
            %optimize hyperparameters
            hyp.cov = zeros(nfeat+1,1);
            hyp = minimize(hyp, @gp,nfe, @infExact, meanfunc, ...
                covfunction, likfunc, xtrs(:,1:nfeat), ytrs);
            %test model on predictions
            [m] = gp(hyp, @infExact, meanfunc, covfunction, likfunc, ...
                xtrs(:,1:nfeat), ytrs, xtes(:,1:nfeat),ytes);
            pts(sksub==j)=m;
        end
        tempms(nfeat)=mean(abs(pts-ytr)); %compute MAE for each feature set
        if nfeat>2
            if tempms(nfeat)>tempms(nfeat-1) && tempms(nfeat)>tempms(nfeat-2) %error starts increasing again
                nfsel(i)=nfeat-2;
                break
            elseif abs(tempms(nfeat)-tempms(nfeat-1))<tolfeat && ...
                    abs(tempms(nfeat)-tempms(nfeat-2))<tolfeat                      %does not vary much anymore
                nfsel(i)=nfeat-2;
                break
            else
               [d1,nfsel(i)]=min(tempms);
            end
        end               
    end
    
    
    %Model data using defined number of features (with permutations)
    %----------------------------------------------------------------------
    
    for p=0:nperm %for each permutation of the targets + the "true" targets
        if p>0
            ind=randperm(length(y));
            yp=y(ind);
        else
            yp=y;
        end 
        ytr=yp(sk~=i,:);
        yte=yp(sk==i,:);
        %optimize hyperparameters
        hyp.cov = zeros(nfsel(i)+1,1);
        hyp = minimize(hyp, @gp,nfe, @infExact, meanfunc, covfunction,...
            likfunc, xtr(:,1:nfsel(i)), ytr);
        %test model on predictions
        [m] = gp(hyp, @infExact, meanfunc, covfunction, likfunc,...
            xtr(:,1:nfsel(i)), ytr, xte(:,1:nfsel(i)),yte);
        if p==0
            pt(sk==i)=m;
        else
            pval(i,p)=mean(abs(yte-m));
        end        
    end       
end
%compute mae
out.CV.mae=mean(abs(pt-yp));
out.CV.pmae=length(find(mean(pval)>out.CV.mae))/nperm;
disp('Computation of the cross validation completed')

nfeat=round(mean(nfsel));
sprintf('Number of features selected: %d',nfeat)
sprintf('Variance of the number of features selected across folds: %d',var(nfsel))

out.featsel.numberperfold=nfsel;
out.featsel.averagenumber=nfeat;
out.featsel.varianceperfold=var(nfsel);

%GP modelling of the whole dataset using selected features
%--------------------------------------------------------------------------

if isfield(in,'hyp')
    hyp = in.hyp;
else %first, optimize the hyperparameters
    hyp.cov = zeros(nfeat+1,1);
    hyp = minimize(hyp, @gp, nfe, @infExact, meanfunc, covfunction, likfunc, x(:,1:nfeat), yp);
    
end

%test on the dataset
[m, s2] = gp(hyp, @infExact, meanfunc, covfunction, likfunc, ...
    x(:,1:nfeat), yp, x(:,1:nfeat));
out.train.mae=mean(abs(m-yp));
out.train.rsquare = rsquare(m,yp);
out.hypcov = exp(hyp.cov);
out.outliers=abs((m-y))./sqrt(s2);  %detect outliers
%compute model evidence, AIC and BIC
nlml = gp(hyp, @infExact, meanfunc, covfunction, likfunc, x(:,1:nfeat), yp);
out.AIC=2*nlml+2*(length(hyp.cov));
out.BIC=2*nlml+(length(hyp.cov))*log(ns);
disp('Training completed')

%Display of the trained model
%--------------------------------------------------------------------------
bi=min(x(:,1:nfeat));
bs=max(x(:,1:nfeat));
if nfeat==1  % create 1D grid
    z=linspace(bi,bs,1000)';
    flag=1;
else
    if nfeat==2 %create 2D grid, rectangular
        flag=2;
        [z1,z2]=meshgrid(linspace(bi(1),bs(1),Ngrid),linspace(bi(2),bs(2),Ngrid));
        z=[z1(:),z2(:)];
    elseif nfeat==3 %create 3D grid, parallelepiped
        flag=3;
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
        for i=1:nfeat
            eval([ 'z = [z, z' num2str(i) '(:)];'])            
        end
        flag=3;
    end 
end

if flag ~= 0
    %compute the predicted values and their variances over the whole hypercube:
    [m s2] = gp(hyp, @infExact, meanfunc, covfunction, likfunc, x(:,1:nfeat), yp, z);
    m = m*maxy + moyy;
    if nfeat > 1
        to=reshape(m, size(z1));
    end
    % Restore the original value of the inputs by unscaling z:
    zz = zeros(size(z,1),nfeat);
    for i = 1:nfeat
        zz(:,i) = z(:,i)*maxx(1,i)+moyx(1,i);
    end
    disp('Computation of the ndgrid completed')
end
if flag==1
    s2 = s2*(maxy^2);
    f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
    figure;
    fill([zz; flipdim(zz,1)], f, [7 7 7]/8)
    hold on; plot(zz, m); plot(xx, yy, '+')
    xlabel(inputnames); ylabel(outputname);
    saveas(gcf,name,'fig')
elseif flag==2
    figure;
    hold on;
    surf(z1*maxx(1,1)+moyx(1,1),z2*maxx(1,2)+moyx(1,2),to);
    plot3(xx(:,1),xx(:,2), yy, '+')
    xlabel(inputnames(1)); ylabel(inputnames(2)) ; zlabel(outputname);
    grid on
    saveas(gcf,name,'fig')
elseif flag ==3
    if isfield(in,'xy') %if the x and y axes are imposed
        idx_sort = [in.xy, 1:(min(in.xy)-1),(min(in.xy)+1):(max(in.xy)-1),(max(in.xy)+1:nfeat)] ;
    else % Sort the variables in terms of their weights or lengthscale in absolute value:
        if ~isempty(strfind(covfunction{:},'ard'))
            [aa idx_sort] = sort(out.hypcov(1:end-1));
        else
            %idx_sort = 1:nv; %otherwise, plot the 2 first variables
            [aa idx_sort] = sort(abs(weights),'descend');
        end
    end
    % Set the nv-2 least relevant variables to their median value:
    med = median(z);
    % Take a slice of the hypercube to plot the two main relevant
    % variables:
    to2 = permute(to,idx_sort);
    med=med(idx_sort);
    itp=[];
    plottext = '';
    for i=1:nfeat
        if i==1 || i==2
            itp=[itp,',:'];
        else
            indi=idx_sort(i);
            vec = linspace(bi(indi),bs(indi),Ngrid);
            [dd,indt] = min(abs(vec-med(i)));
            itp=[itp,[',',num2str(indt)]];
            plottext = strcat(plottext,{' '},inputnames(indi), {' = '}, num2str(vec(indt)*maxx(1,indi)+moyx(1,indi)), {'; '});
        end
    end
    eval([ 'to2 = to2(' itp(2:end) ');'])
    vecx=linspace(bi(idx_sort(1)),bs(idx_sort(1)),Ngrid)'*maxx(1,idx_sort(1))+moyx(1,idx_sort(1));
    vecy=linspace(bi(idx_sort(2)),bs(idx_sort(2)),Ngrid)'*maxx(1,idx_sort(2))+moyx(1,idx_sort(2));
    figure;
    surf(vecy,vecx,to2);
    ylabel(inputnames(idx_sort(1))); xlabel(inputnames(idx_sort(2))) ; 
    zlabel(outputname); title(plottext);
    grid on
    saveas(gcf,name,'fig')
end
out.ss=ss;
save(['GP_modelling_dataset_',num2str(name),'.mat'],'out','out')  

% %Linear kernel: determine weights for each variable
% hyp2 = minimize(hyp2, @gp, nfe, @infExact, meanfunc, {@covLINone}, likfunc, x, yp);
% [d1,d2,d3,d4,d5,post] = gp(hyp2, @infExact, meanfunc, {@covLINone}, likfunc, x, yp, x);
% alpha=post.alpha;
% weights=x'*alpha;
% out.weights=abs(weights/norm(weights));
% disp('computation of weights completed (linear kernel)')   
    

