function [RptgenML_CReport1] = buildReport(in,out)
%Report generation
 
% Create RptgenML.CReport
RptgenML_CReport1 = RptgenML.CReport('Description',...
['Results explanation'],...
'Stylesheet','html-!SingleClearTitleTocLot',...
'FilenameName',[datestr(now,'yymmdd HHMM') ' - ' in.name],...
'FilenameType','other',...
'DirectoryType','pwd');
% setedit(RptgenML_CReport1);
 
% Create rptgen.cfr_titlepage
rptgen_cfr_titlepage1 = rptgen.cfr_titlepage('Copyright_Date','2015',...
'Author','S. Quoilin, A. Legros',...
'Title',['Results analysis of the dataset "' in.name '"']);

rptgen_cfr_text2 = rptgen.cfr_text;
set(rptgen_cfr_titlepage1,'LegalNoticeComp',rptgen_cfr_text2);
rptgen_cfr_image1 = rptgen.cfr_image('FileName','','MaxViewportSize',[7 9]);
set(rptgen_cfr_titlepage1,'ImageComp',rptgen_cfr_image1);
setParent(rptgen_cfr_titlepage1,RptgenML_CReport1);
 
%% Create the first section
rptgen_cfr_section1 = rptgen.cfr_section('SectionTitle',...
['DISCLAMER']);
setParent(rptgen_cfr_section1,RptgenML_CReport1);
rptgen_cfr_text1 = rptgen.cfr_text('Content',...
    ['The information provided below is provided without garantee and should be considered as an illustrative example of how GPExp results can be exploited. '...
    'Most of the qualitative analysis is based on experience and cannot be considered as firm rules. Depending on the dataset, the imposed euristic '...
    'boundaries can vary significantly']);
setParent(rptgen_cfr_text1,rptgen_cfr_section1);

%% Create the second section

rptgen_cfr_section2 = rptgen.cfr_section('SectionTitle','Dataset Information');
setParent(rptgen_cfr_section2,RptgenML_CReport1);

%% First subsection : MAIN RESULTS

rptgen_cfr_subSection21 = rptgen.cfr_section('SectionTitle','Main Results');
setParent(rptgen_cfr_subSection21,rptgen_cfr_section2);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content','These results are stored in the out.CV and out.train variables.');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content','In training mode (with all the data points):');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

Table1 = {'MARE [%]', num2str(out.train.mae*100);
         'R� [%]', num2str(out.train.rsquare*100)
         'RMSE [-]', num2str(out.train.rmse)}'; 

rptgen_cfr_table1 = rptgen.cfr_table('Source',Table1,'TableTitle','Training results','AllAlign','center','ColumnWidths',[2 2 2],'isPgwide',false);
setParent(rptgen_cfr_table1,rptgen_cfr_subSection21);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content','*** Bla bla sur la d�finition du MARE, R� et RMSE.');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

if isfield(out,'CV')
    
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content','In cross-validation:');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

    Table2 = {'MARE [%]', num2str(out.CV.mae*100);
             'R� [%]', num2str(out.CV.rsquare*100)
             'RMSE [-]', num2str(out.CV.rmse)}'; 

    rptgen_cfr_table1 = rptgen.cfr_table('Source',Table2,'TableTitle','Training results','AllAlign','center','ColumnWidths',[2 2 2],'isPgwide',false);
    setParent(rptgen_cfr_table1,rptgen_cfr_subSection21);
else
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content','No cross-validation was Performed');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
end

string = [in.considered_inputs{1}];
for i = 2:length(in.considered_inputs)
    string = [string ', ' in.considered_inputs{i}];
end

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content','These results can be interpreted as follows:');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content',['   - The lowest average error that could be reached by a model predicting ' in.considered_output{1} ...
    ' as a function of ' string ' is ' num2str(out.train.mae*100) ' %'],'Color','blue');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

if isfield(out,'CV')
    
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content',['   - When predicting values that are outside of the data, an average error of ' ...
        num2str(out.CV.mae*100) ' % can be expected'],'Color','blue');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);    
end

rptgen_cfr_image2 = rptgen.cfr_image('Filename',[in.name '_fig2.png'],...
'MaxViewportSize',[7 9],...
'DocHorizAlign','center','isTitle','local',...
'Title',['Predicted vs measured values of ' in.considered_output{1}]);
setParent(rptgen_cfr_image2,rptgen_cfr_subSection21);

%% Second subsection : DETECTION OF OVERFITTING

rptgen_cfr_subSection21 = rptgen.cfr_section('SectionTitle','Detection of overfitting');
setParent(rptgen_cfr_subSection21,rptgen_cfr_section2);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content',['Overfitting can be detected by comparing the error in training mode and in cross-' ...
    'validation: in case of overfitting the train error should remain low, while the '...
    'CV error should increase significantly.']);
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content',['However, there is no firm quantitative rule to detec overfitting. This analysis'...
    'therefore provides a warning, which should be checked visually by plotting the function']);
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

if isfield(out,'CV')
    ratio_error = out.CV.mae/out.train.mae;
    
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content',['The ratio between the error in CV and training is ' num2str(ratio_error,3)],'Color','blue');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
    
    if ratio_error < 2
        rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
        rptgen_cfr_text1 = rptgen.cfr_text('Content','This value relatively low, which tends to indicate that there is no overfitting','Color','blue');
        set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
        setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
    elseif ratio_error < 4
        rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
        rptgen_cfr_text1 = rptgen.cfr_text('Content',['This value is comprised between 2 and 4, which might indicate some overfitting! ' ...
            'If overfitting is visually detected, try reducing the number of inputs, or increasing the ' ...
            'number of data samples'],'Color','blue');
        set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
        setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
    else 
        rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
        rptgen_cfr_text1 = rptgen.cfr_text('Content',['This value is higher than 4, which indiates a high chance of overfitting. '...
            'The 3D plot should be carefully checked: unwanted wiggles and noise in the function are '...
            'related to overfitting. In that case, no good model could be found with the input data. '...
            'Try reducing the number of inputs, or increase the number of data points'],'Color','red');
        set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
        setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
    end
else
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content',['No cross-validation was performed, impossible to detect overfitting other ' ...
        'than visually, by displaying the 3D plot of the function'],'Color','blue');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
end

rptgen_cfr_image2 = rptgen.cfr_image('Filename',[in.name '_fig1.png'],...
'MaxViewportSize',[7 9],...
'DocHorizAlign','center','isTitle','local',...
'Title','Detection of overfitting');
setParent(rptgen_cfr_image2,rptgen_cfr_subSection21);
    
%% Third subsection : RELEVANCE OF THE SELECTED INPUTS (PERMUTATIONS)

rptgen_cfr_subSection21 = rptgen.cfr_section('SectionTitle','Relevance of the selected inputs (Permutations)');
setParent(rptgen_cfr_subSection21,rptgen_cfr_section2);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content','(NB: This result is stored in the out.CV.pmae variable)');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content',['The dependency of the output with the given inputs is checked by comparing '...
'the mean average error (in cross-validation) of the dataset with a random '...
'dataset (the same with randomly permuted outputs. The reported statistics '...
'in the CV.pmae variable, corresponding to the probability of having a random '...
'dataset performing better in terms of mean average error than the actual '...
'dataset.']);
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content',['A pmae lower than 5 percent indicates that there is a significant correlation between '...
'inputs and outputs']);
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

if in.perm < 1
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content','Permutations were not used in the present analysis','Color','blue');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
elseif in.perm < 20
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content',['Warning: The number of permutations is quite low (' num2str(in.perm) '),' ...
        'the provided pmae might not have a sufficient statistical relevance'],'Color','red');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
    
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content',['Computed pmae value: ' num2str(out.CV.pmae*100) ' %'],'Color','blue');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
else
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content',['Computed pmae value: ' num2str(out.CV.pmae*100) ' %'],'Color','blue');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
end

%% Fourth subsection : OUTLIERS

rptgen_cfr_subSection21 = rptgen.cfr_section('SectionTitle','Outliers');
setParent(rptgen_cfr_subSection21,rptgen_cfr_section2);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content','(NB: These results are stored in the out.outliers vector)');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content',['The posterior of the Gaussian Process provides the likelihood function, '...
'with its mean and standard deviation. For each data point, the Grubbs ' ...
'test for outliers can therefore be applied: the error is expressed as the '...
'a multiple of the standard deviation at that particular point. According '...
'to the normal distribution hypothesis, a multiple higher than 1.96 '...
'indicates a significance level lower than 5 percent '...
'Users should also refer to the error plot to visualize the repartition of '...
'the error on the normal distribution']);
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content','The data point with the highest probability to be an outlier is:');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

[val idx] = max(abs(out.outliers)); 

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content',['Data point number ' num2str(idx) ', with an error ' num2str(val) ' times higher than the standard ' ...
'deviation of the gaussian process function at that particular point'],'Color','blue');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
rptgen_cfr_text1 = rptgen.cfr_text('Content',['The following data points present a significance level lower than 5 percent. '...
'They are therefore likely to be outliers:'],'Color','blue');
set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);

idx = find(abs(out.outliers)>=1.96);
if isempty(idx)
    rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
    rptgen_cfr_text1 = rptgen.cfr_text('Content',' - None','Color','blue');
    set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
    setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
else
    for i = 1:length(idx)
        rptgen_cfr_paragraph1 = rptgen.cfr_paragraph;
        rptgen_cfr_text1 = rptgen.cfr_text('Content',[' - Data point number ' num2str(idx(i)) ': ' num2str(abs(out.outliers(idx(i))))],'Color','blue');
        set(rptgen_cfr_paragraph1,'ParaTextComp',rptgen_cfr_text1);
        setParent(rptgen_cfr_paragraph1,rptgen_cfr_subSection21);
    end
end

rptgen_cfr_image2 = rptgen.cfr_image('Filename',[in.name '_fig3.png'],'isTitle','local',...
'Title','Detection of outliers: Normal distribution of the error, with 5% significance intervals',...
'MaxViewportSize',[7 9],...
'DocHorizAlign','center');
setParent(rptgen_cfr_image2,rptgen_cfr_subSection21);

