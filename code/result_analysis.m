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

write(fid,in.description);

write(fid,{
    ' '
    ' '
    'MAIN RESULTS'
    ' '
    'In training mode (with all the data points):'
    ['Mean relative error: ' num2str(out.train.mae*100) '%']
    ['Coefficient of determination (R square): '   num2str(out.train.rsquare*100) '%']
    ['Root mean square error (RMSE): ' num2str(out.train.rmse)]    
    });
if isfield(out,'CV')
    write(fid,{
        ' '
        'In cross-validation:'
        ['Mean relative error: ' num2str(out.CV.mae*100) '%']
        ['Coefficient of determination (R square): '   num2str(out.CV.rsquare*100) '%']
        ['Root mean square error (RMSE): ' num2str(out.CV.rmse)]    
        });
else
    write(fid,{
        ' '
        'No cross-validation was Performed'
        });
end


write(fid,{
    ' '
    ' '
    'DETECTION OF OVERFITTING'
    ' '
    });
if isfield(out,'CV')
    write(fid,{
        ' '
        'In cross-validation:'
        ['Mean relative error: ' num2str(out.CV.mae*100) '%']
        ['Coefficient of determination (R square): '   num2str(out.CV.rsquare*100) '%']
        ['Root mean square error (RMSE): ' num2str(out.CV.rmse)]    
        });
else
    write(fid,{
        'No cross-validation was Performed, impossible to detect overfitting other'
        'than visually, by looking at the plot'
        });
end




write(fid,{
'  '
' '
'RELEVANCE OF THE SELECTED INPUTS (PERMUTATIONS)'
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
'OUTLIERS'
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
'RE-USE OF THE MODEL RESULTS'
' '
'It is not necessary to perform a whole new optimization to re-use the results'
'e.g. for plotting, or to predict the function output for a given input set'
'Most of the qualitative analysis is based on experience and cannot be'
'considered as firm rules. Depending on the dataset, the imposed euristic'
'boundaries can vary significantly'
' '
' '
});



end

