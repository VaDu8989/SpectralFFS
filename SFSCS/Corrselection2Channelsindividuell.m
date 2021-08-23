

function varargout = Corrselection2Channelsindividuell(varargin)
% CORRSELECTION2CHANNELSINDIVIDUELL MATLAB code for Corrselection2Channelsindividuell.fig
%      CORRSELECTION2CHANNELSINDIVIDUELL, by itself, creates a new CORRSELECTION2CHANNELSINDIVIDUELL or raises the existing
%      singleton*.
%
%      H = CORRSELECTION2CHANNELSINDIVIDUELL returns the handle to a new CORRSELECTION2CHANNELSINDIVIDUELL or the handle to
%      the existing singleton*.
%
%      CORRSELECTION2CHANNELSINDIVIDUELL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRSELECTION2CHANNELSINDIVIDUELL.M with the given input arguments.
%
%      CORRSELECTION2CHANNELSINDIVIDUELL('Property','Value',...) creates a new CORRSELECTION2CHANNELSINDIVIDUELL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Corrselection2Channelsindividuell_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Corrselection2Channelsindividuell_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Corrselection2Channelsindividuell

% Last Modified by GUIDE v2.5 17-Jan-2017 10:26:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Corrselection2Channelsindividuell_OpeningFcn, ...
                   'gui_OutputFcn',  @Corrselection2Channelsindividuell_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Corrselection2Channelsindividuell is made visible.
function Corrselection2Channelsindividuell_OpeningFcn(hObject, eventdata, handles, varargin)
global curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  corfit meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Corrselection2Channelsindividuell (see VARARGIN)

% Choose default command line output for Corrselection2Channelsindividuell
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.slider3,'Max',size(correlationcurvesCh1,2)-1);
set(handles.slider3,'SliderStep',[1/(size(correlationcurvesCh1,2)-2) 1/(size(correlationcurvesCh1,2)-2)]);
correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;


for kl=2:size(correlationcurvesCh1,2)
if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
end
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
end
end
meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);

set(handles.checkbox2,'Value',curveincl1(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl2(get(handles.slider3,'Value'),1));
    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1)],'r-')
     ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     hold (handles.axes2,'off')
     
     plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(handles.slider3,'Value')+1))
    hold (handles.axes3,'on')
    plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1)],'r-')
     ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
     hold (handles.axes3,'off')
     semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'Marker','+','LineStyle','none','Color',[171 221 164]/255)
     hold (handles.axes1,'on')
     semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'Marker','+','LineStyle','none','Color',[215 25 28]/255)
     hold (handles.axes1,'on')
     semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'Marker','+','LineStyle','none','Color',[43 131 186]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'LineStyle','-','Color',[171 221 164]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesCh2(:,1),meancorrcurve(:,2),'LineStyle','-','Color',[215 25 28]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesChCC(:,1),meancorrcurveCC,'LineStyle','-','Color',[43 131 186]/255)
 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% For now: Don't plot fits of CFs
%semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
% For now: Don't set axis manually
%axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
hold (handles.axes1,'off')

% UIWAIT makes Corrselection2Channelsindividuell wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Corrselection2Channelsindividuell_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
global curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  corfit meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


if get(hObject,'Value')<size(correlationcurvesCh1,2)
    
    correlationcurves2Ch1=correlationcurvesCh1;
    correlationcurves2Ch2=correlationcurvesCh2;
    correlationcurves2ChCC=correlationcurvesChCC;

for kl=2:size(correlationcurvesCh1,2)
if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
end
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
end
end
meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);

set(handles.checkbox2,'Value',curveincl1(get(hObject,'Value'),1));
set(handles.checkbox3,'Value',curveincl2(get(hObject,'Value'),1));
% get(hObject,'Min')
% get(hObject,'Max')
%      get(hObject,'Value')+1
    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(hObject,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(hObject,'Value')+1),1) mean(ItraceIICh1(:,get(hObject,'Value')+1),1)],'r-')
     ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     hold (handles.axes2,'off')
     
     plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(hObject,'Value')+1))
    hold (handles.axes3,'on')
    plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(hObject,'Value')+1),1) mean(ItraceIICh2(:,get(hObject,'Value')+1),1)],'r-')
     ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
     hold (handles.axes3,'off')

semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(hObject,'Value')+1),'Marker','+','LineStyle','none','Color',[171 221 164]/255)
     hold (handles.axes1,'on')
     semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(hObject,'Value')+1),'Marker','+','LineStyle','none','Color',[215 25 28]/255)
     hold (handles.axes1,'on')
     semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(hObject,'Value')+1),'Marker','+','LineStyle','none','Color',[43 131 186]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'LineStyle','-','Color',[171 221 164]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesCh2(:,1),meancorrcurve(:,2),'LineStyle','-','Color',[215 25 28]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesChCC(:,1),meancorrcurveCC,'LineStyle','-','Color',[43 131 186]/255)





 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% For now: Don't plot fits of CFs
%semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(hObject,'Value')+1),'k-','LineWidth',2)
% For now: Don't set axis manually
%axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(hObject,'Value')+1))) max(max(correlationcurves(1:10,get(hObject,'Value')+1)))]) 
hold (handles.axes1,'off')
end

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)

global curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  corfit meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC

% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



curveincl1(get(handles.slider3,'value'),1)=get(hObject,'Value');
correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

for kl=2:size(correlationcurvesCh1,2)
if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
end
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
end
end
meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);
  
set(hObject,'Value',curveincl1(get(handles.slider3,'Value'),1));
    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1)],'r-')
     ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     hold (handles.axes2,'off')
     
      plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(handles.slider3,'Value')+1))
    hold (handles.axes3,'on')
    plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1)],'r-')
     ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
     hold (handles.axes3,'off')
   
semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'Marker','+','LineStyle','none','Color',[171 221 164]/255)
     hold (handles.axes1,'on')
     semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'Marker','+','LineStyle','none','Color',[215 25 28]/255)
     hold (handles.axes1,'on')
     semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'Marker','+','LineStyle','none','Color',[43 131 186]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'LineStyle','-','Color',[171 221 164]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesCh2(:,1),meancorrcurve(:,2),'LineStyle','-','Color',[215 25 28]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesChCC(:,1),meancorrcurveCC,'LineStyle','-','Color',[43 131 186]/255)




 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% For now: Don't plot fits of CFs
%semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
% For now: Don't set axis manually
%axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
hold (handles.axes1,'off')


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global goon
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
goon=1;
close Corrselection2Channelsindividuell


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
global curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  corfit meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC

% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



curveincl2(get(handles.slider3,'value'),1)=get(hObject,'Value');
correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

for kl=2:size(correlationcurvesCh1,2)
if curveincl1(kl-1)==0
    correlationcurves2Ch1(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
end
if curveincl2(kl-1)==0
correlationcurves2Ch2(:,kl)=NaN;
correlationcurves2ChCC(:,kl)=NaN;
end
end
meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);
  
set(hObject,'Value',curveincl2(get(handles.slider3,'Value'),1));
    plot(handles.axes2,ItraceIICh1(:,1),ItraceIICh1(:,get(handles.slider3,'Value')+1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ ItraceIICh1(1,1) ItraceIICh1(size(ItraceIICh1,1),1) ],[mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh1(:,get(handles.slider3,'Value')+1),1)],'r-')
     ylim(handles.axes2,[min(ItraceIICh1(:,get(handles.slider3,'Value')+1)) max(ItraceIICh1(:,get(handles.slider3,'Value')+1))])
     hold (handles.axes2,'off')
     
      plot(handles.axes3,ItraceIICh2(:,1),ItraceIICh2(:,get(handles.slider3,'Value')+1))
    hold (handles.axes3,'on')
    plot(handles.axes3,[ ItraceIICh2(1,1) ItraceIICh2(size(ItraceIICh2,1),1) ],[mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1) mean(ItraceIICh2(:,get(handles.slider3,'Value')+1),1)],'r-')
     ylim(handles.axes3,[min(ItraceIICh2(:,get(handles.slider3,'Value')+1)) max(ItraceIICh2(:,get(handles.slider3,'Value')+1))])
     hold (handles.axes3,'off')

semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'Marker','+','LineStyle','none','Color',[171 221 164]/255)
     hold (handles.axes1,'on')
     semilogx(handles.axes1,correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'Marker','+','LineStyle','none','Color',[215 25 28]/255)
     hold (handles.axes1,'on')
     semilogx(handles.axes1,correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'Marker','+','LineStyle','none','Color',[43 131 186]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'LineStyle','-','Color',[171 221 164]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesCh2(:,1),meancorrcurve(:,2),'LineStyle','-','Color',[215 25 28]/255)
hold (handles.axes1,'on')
semilogx(handles.axes1,correlationcurvesChCC(:,1),meancorrcurveCC,'LineStyle','-','Color',[43 131 186]/255)


 xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
% For now: Don't plot fits of CFs
%semilogx(handles.axes1,correlationcurves(:,1),corfit(:,get(handles.slider3,'Value')+1),'k-','LineWidth',2)
% For now: Don't set axis manually
%axis(handles.axes1,[min(min(correlationcurves(:,1))) max(max(correlationcurves(:,1))) min(min(correlationcurves(:,get(handles.slider3,'Value')+1))) max(max(correlationcurves(1:10,get(handles.slider3,'Value')+1)))]) 
hold (handles.axes1,'off')
