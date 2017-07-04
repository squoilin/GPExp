function [ success ] = result_analysis(in,out)
%RESULTS_ANALYSIS GPExp result analysis tool
%   This function provides insights and interpretations of the GPexp
%   results. Most of the proposed analysis is qualitative and can be use as
%   an illustrative example of how the data can be exploited.
%   The text results are displayed in the Matlab interface, and saved to
%   the log.txt file in the working directory.
%   Parameters:
%   in : the input structure of GPExp
%   out: the output (results) structure of GPExp
%
%   Copyright (c) 2013-2015 Sylvain Quoilin & Jessica Schrouff. 
%   All rights reserved.

fid = fopen('log.txt','w');

write(fid,{
' '
['GPExp: results analysis of the "' in.name '" dataset']
' '
' '
});

write(fid,{    
'DISCLAMER'
' '
'The information provided below is provided without garantee and should'
'be considered as an illustrative example of how GPExp results can be exploited'
'Most of the qualitative analysis is based on experience and cannot be'
'considered as firm rules. Depending on the dataset, the imposed euristic'
'boundaries can vary significantly'
' '
' '
});

write(fid,{
'DATASET INFORMATION'
' '
});

write(fid,in.description,'blue');

write(fid,{
    ' '
    ' '
    'MAIN RESULTS'
    '(NB: These results are stored in the out.CV and out.train variables)'
    ' '
    'Whole training set (i.e. with all the data points):'
    });

write(fid,{
    ['Normalized mean absolute error: ' num2str(out.train.mae*100) ' %%']
    ['Coefficient of determination (R square): '   num2str(out.train.rsquare*100) ' %%']
    ['Normalized root mean square error (RMSE): ' num2str(out.train.rmse)]    
    },'blue');

 write(fid,{
        ' '
        'In cross-validation:'
    });

if isfield(out,'CV')
write(fid,{
        ['Normalized mean absolute error: ' num2str(out.CV.mae*100) ' %%']
        ['Coefficient of determination (R square): '   num2str(out.CV.rsquare*100) ' %%']
        ['Normalized root mean square error (RMSE): ' num2str(out.CV.rmse)]    
        },'blue');
else
    write(fid,{
        ' '
        'No cross-validation was Performed'
        },'blue');
end

string = [in.considered_inputs{1}];
for i = 2:length(in.considered_inputs)
    string = [string ', ' in.considered_inputs{i}];
end

write(fid,{
        ' '
        ['These results can be interpreted as follows:']
        ' '
        });

write(fid,{
        ['The lowest average error that could be reached by a model predicting ' in.considered_output{1}]
        ['as a function of ' string ' is ' num2str(out.train.mae*100) ' %%']   
        },'blue');

if isfield(out,'CV')
    write(fid,{
        ' '
        ['When predicting values that are outside of the data, an average error of']
        [num2str(out.CV.mae*100) ' %% can be expected']   
        },'blue');    
end


write(fid,{
    ' '
    ' '
    'DETECTION OF OVERFITTING'
    ' '
    'Overfitting can be detected by comparing the error in training mode and in cross-'
    'validation: in case of overfitting the train error should remain low, while the'
    'CV error should increase significantly'
    'However, there is no firm quantitative rule to detec overfitting. This analysis'
    'therefore provides a warning, which should be checked visually by plotting the function'
    ' '
    });


if isfield(out,'CV')
    ratio_error = out.CV.mae/out.train.mae;
    write(fid,{['The ratio between the error in CV and training is ' num2str(ratio_error,3)]},'blue');
    if ratio_error < 2
        write(fid,{'This value relatively low, which tends to indicate that there is no overfitting'},'blue');
    elseif ratio_error < 4
        write(fid,{'This value is comprised between 2 and 4, which might indicate some overfitting!'
            'If overfitting is visually detected, try reducing the number of inputs, or increasing the'
            'number of data samples'},'blue');
    else 
        write(fid,{'This value is higher than 4, which indiates a high chance of overfitting'
            'The 3D plot should be carefully checked: unwanted wiggles and noise in the function are'
            'related to overfitting. In that case, no good model could be found with the input data'
            'Try reducing the number of inputs, or increase the number of data points'},'red');
    end
else
    write(fid,{
        'No cross-validation was Performed, impossible to detect overfitting other'
        'than visually, by displaying the 3D plot of the function'
        },'blue');
end

% To be done: add a criteria with the smallest lengthscale (equivalent to
% visual check?)



write(fid,{
'  '
' '
'RELEVANCE OF THE SELECTED INPUTS (PERMUTATIONS)'
'(NB: This result is stored in the out.CV.pmae variable)'
' '
'The dependency of the output with the given inputs is checked by comparing'
'the mean average error (in cross-validation) of the dataset with a random'
'dataset (the same with randomly permuted outputs. The reported statistics'
'in the CV.pmae variable, corresponding to the probability of having a random'
'dataset performing better in terms of mean average error than the actual'
'dataset.'
'A pmae lower than 5 percent indicates that there is a significant correlation between'
'inputs and outputs'
' '
});

if in.perm < 1
    write(fid,{
        'Permutations were not used in the present analysis'
        },'blue');
elseif in.perm < 100
    write(fid,{
        ['Warning: The number of permutations is quite low (' num2str(in.perm) '),']
        'the provided pmae might not have a sufficient statistical relevance'
        ' '
        },'red');
    write(fid,{
        ['Computed pmae value: ' num2str(out.CV.pmae*100) ' %']
        },'blue');
else
    write(fid,{    
    ['Computed pmae value: ' num2str(out.CV.pmae*100) ' %']
    },'blue');
end
    



write(fid,{
' '
' '
'OUTLIERS'
'(NB: These results are stored in the out.outliers vector)'
' '
'The posterior of the Gaussian Process provides the likelihood function,'
'with its mean and standard deviation. For each data point, the Grubbs' 
'test for outliers can therefore be applied: the error is expressed as the'
'a multiple of the standard deviation at that particular point. According '
'to the normal distribution hypothesis, a multiple higher than 1.96 '
'indicates a significance level lower than 5 percent'
'Users should also refer to the error plot to visualize the repartition of '
'the error on the normal distribution'
' '
'The data point with the highest probability to be an outlier is:'
});

[val idx] = max(abs(out.outliers)); 

write(fid,{
['Data point number ' num2str(idx) ', with an error ' num2str(val) ' times higher than the standard']
'deviation of the gaussian process function at that particular point'
},'blue');

write(fid,{
' '
'The following data points present a significance level lower than 5 percent.'
'They are therefore likely to be outliers:'
});
idx = find(abs(out.outliers)>=1.96);
if isempty(idx)
    write(fid,{'None'},'blue');
else
    for i = 1:length(idx)
         write(fid,{['Data point number ' num2str(idx(i)) ': ' num2str(abs(out.outliers(idx(i))))]},'blue');
    end
end


write(fid,{
' '
' '
'LENGTHSCALES AND SENSITIVITY ANALYSIS'
' '
'Feature selection (i.e. determining the relevance of each input variable'
'can be performed based on the optimal lengthscales of the ARD kernel'
'The lengthscales are a good indicator of the relevance of a particular '
'input to predict the output. High lengthscale values indicate that the'
'output varies slowly in the direction of the selected input. The sensitivity'
'is therefore low. Extermely high lengthscales indicate that the optimization'
'of the marginal likelyhool leads to neglect the influence of the specified variable'
});
write(fid,{' '},'blue');

for i=1:length(in.considered_inputs)
    string = [in.considered_inputs{i},': ',num2str(out.hypcov(i))];
    write(fid,{string},'blue');
end
write(fid,{' '},'blue');
write(fid,{['Prior Likelyhood: ', num2str(out.hypcov(end))]},'blue');

[minval,minpos] = min(out.hypcov(1:end-1));
[maxval,maxpos] = max(out.hypcov(1:end-1));

write(fid,{' '},'blue');

write(fid,{['The most relevant variable seems to be ', in.considered_inputs{minpos}]},'blue');
write(fid,{['The least relevant variable seems to be ', in.considered_inputs{maxpos}]},'blue');

write(fid,{
'   '
' '
' '
},'blue');





write(fid,{
' '
' '
'PREDICTION'
' '
'Once the analysis has been performed and the results are satisfying (i.e. '
'no overfitting, outliers have been removed, accuracy is sufficient, etc),'
'the result files can be used for prediction, i.e. to predict the output ' 
'of a set of inputs outside of the data.'
'This can be done by :'
'1. loading the "in" and "out" structures from the .mat result analysis file:'
});
write(fid,{'   "load data-file.mat"'},'blue');
write(fid,{
    '2. Use the GPExp prediction function by assigning a value to each input:'
});

string = [];
for i=1:length(in.considered_inputs)
    string = [string '''' in.considered_inputs{i} ''', value' num2str(i) ', '];
end
string = string(1:end-2);    % removing end comma
write(fid,{
['   "GP_prediction(in,results, ' string ')"']
' '
' '
},'blue');

end

