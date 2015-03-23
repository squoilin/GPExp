figure('Visible','off')

axis([miny,maxy,miny,maxy])
hold on
plot(y,handles.results.train.y_pred,'r*')
plot([miny,maxy],[miny,maxy],'k')
if isfield(handles.results,'CV')
    plot(y,handles.results.CV.y_pred,'b+')
    legend('Train (i.e. with all data samples)','45 deg','Cross-Validation','Location','best')
else
    legend('Train (i.e. with all data samples)','45 deg','Location','best')
end
% title(['Predicted vs measured values of ' handles.in.considered_output])
saveas(gcf,[handles.in.name '_fig2'],'png')

delete(gcf)