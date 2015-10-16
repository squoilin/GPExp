function [ success ] = plot_prediction(in,out)
%PLOT_RESULTS Plotting of the GPExp simulation results
%   Required inputs:
%   in : the input structure of GPExp
%   out: the output (results) structure of GPExp

% plot of predicted vs. measured variables

miny = min(in.y) - (max(in.y)-min(in.y))/15;
maxy = max(in.y) + (max(in.y)-min(in.y))/15;

axis([miny,maxy,miny,maxy])
hold on
plot(in.y,out.train.y_pred,'r*')
plot([miny,maxy],[miny,maxy],'k')
if isfield(out,'CV')
    plot(in.y,out.CV.y_pred,'b+')
    legend('Train (i.e. with all data samples)','45 deg','Cross-Validation','Location','NorthWest')
else
    legend('Train (i.e. with all data samples)','45 deg','Location','NorthWest')
end
title(['Pred. vs meas. values of ' in.considered_output{1}])

success = true;
end
