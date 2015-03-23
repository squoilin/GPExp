set(handles.ErrorDisplay,'ForegroundColor','k')
set(handles.ConsideredInputListbox,'Enable','off')
set(handles.ConsideredInputListbox,'String','')
set(handles.ConsideredOutputListbox,'Enable','off')
set(handles.ConsideredOutputListbox,'String','')
set(handles.TimeVariableCheckbox,'Enable','off')
set(handles.kfoldTextbox,'Enable','off')
set(handles.PermutationTextbox,'Enable','off')
set(handles.RunButton,'Enable','off')
set(handles.AnalysisStatusButton,'String','')
set(handles.ResultTable,'Visible','off');
set(handles.Textbox1,'String','');
set(handles.Textbox2,'String','');

cla(handles.MainPlot)
cla(handles.SecondPlot)
cla(handles.ThirdPlot)