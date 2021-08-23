function varargout = Corrselection4Channels(varargin)
% CORRSELECTION4CHANNELS MATLAB code for Corrselection4Channels.fig
%      CORRSELECTION4CHANNELS, by itself, creates a new CORRSELECTION4CHANNELS or raises the existing
%      singleton*.
%
%      H = CORRSELECTION4CHANNELS returns the handle to a new CORRSELECTION4CHANNELS or the handle to
%      the existing singleton*.
%
%      CORRSELECTION4CHANNELS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRSELECTION4CHANNELS.M with the given input arguments.
%
%      CORRSELECTION4CHANNELS('Property','Value',...) creates a new CORRSELECTION4CHANNELS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Corrselection4Channels_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Corrselection4Channels_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Corrselection4Channels

% Last Modified by GUIDE v2.5 27-Apr-2020 15:39:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Corrselection4Channels_OpeningFcn, ...
                   'gui_OutputFcn',  @Corrselection4Channels_OutputFcn, ...
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

% --- Executes just before Corrselection4Channels is made visible.
function Corrselection4Channels_OpeningFcn(hObject, eventdata, handles, varargin)
global x0 y0 fixed curveincl_Chs correlationcurves_Chs correlationcurvesCC_Chs ...
sigmascurves_Chs sigmascurvesCC_Chs meancorrcurve meancorrcurveCC ItraceII_Chs...
correlationcurves_Chs_2 correlationcurvesCC_Chs_2 sigmascurves_Chs_2 ...
sigmascurvesCC_Chs_2 yuplimit ylowlimit colormatrix_Chs colormatrixCC_Chs

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Corrselection4Channels (see VARARGIN)

% Choose default command line output for Corrselection4Channels
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.slider3,'Max',size(correlationcurves_Chs,2)-1);
set(handles.slider3,'SliderStep',[1/(size(correlationcurves_Chs,2)-2) 1/(size(correlationcurves_Chs,2)-2)]);
correlationcurves_Chs_2=correlationcurves_Chs;
correlationcurvesCC_Chs_2=correlationcurvesCC_Chs;

sigmascurves_Chs_2=sigmascurves_Chs;
sigmascurvesCC_Chs_2=sigmascurvesCC_Chs;

for kl=2:size(correlationcurves_Chs,2)
    ccindex=1;
    for i=1:4
        if curveincl_Chs(kl-1,i)==0
            correlationcurves_Chs_2(:,kl,i)=NaN;
            sigmascurves_Chs_2(:,kl,i)=NaN;
        end
        if i<4
            for j=i+1:4
                if curveincl_Chs(kl-1,i)==0 || curveincl_Chs(kl-1,j)==0
                correlationcurvesCC_Chs_2(:,kl,ccindex)=NaN;
                sigmascurvesCC_Chs_2(:,kl,ccindex)=NaN;
                end
                ccindex=ccindex+1;
            end
        end
    end
end
meancorrcurve=squeeze(nanmean(correlationcurves_Chs_2(:,2:end,:),2));
meancorrcurveCC=squeeze(nanmean(correlationcurvesCC_Chs_2(:,2:end,:),2));

set(handles.checkbox2,'Value',curveincl_Chs(get(handles.slider3,'Value'),1));
set(handles.checkbox3,'Value',curveincl_Chs(get(handles.slider3,'Value'),2));
set(handles.checkbox4,'Value',curveincl_Chs(get(handles.slider3,'Value'),3));
set(handles.checkbox5,'Value',curveincl_Chs(get(handles.slider3,'Value'),4));

plot(handles.axes2,ItraceII_Chs(:,1,1),ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))
hold (handles.axes2,'on')
plot(handles.axes2,[ItraceII_Chs(1,1,1) ItraceII_Chs(size(ItraceII_Chs,1),1,1) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1),1)],'r-')
ylim(handles.axes2,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))])
hold (handles.axes2,'off')

plot(handles.axes3,ItraceII_Chs(:,1,2),ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))
hold (handles.axes3,'on')
plot(handles.axes3,[ItraceII_Chs(1,1,2) ItraceII_Chs(size(ItraceII_Chs,1),1,2) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1)],'r-')
ylim(handles.axes3,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))])
hold (handles.axes3,'off')

plot(handles.axes6,ItraceII_Chs(:,1,3),ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))
hold (handles.axes6,'on')
plot(handles.axes6,[ItraceII_Chs(1,1,3) ItraceII_Chs(size(ItraceII_Chs,1),1,3)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1)],'r-')
ylim(handles.axes6,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))])
hold (handles.axes6,'off')

plot(handles.axes7,ItraceII_Chs(:,1,4),ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))
hold (handles.axes7,'on')
plot(handles.axes7,[ItraceII_Chs(1,1,4) ItraceII_Chs(size(ItraceII_Chs,1),1,4)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1)],'r-')
ylim(handles.axes7,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))])
hold (handles.axes7,'off')
     
yuplimit=max(max(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));
ylowlimit=min(min(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));

for i=1:4
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),correlationcurves_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),meancorrcurve(:,i),'LineStyle','-','Color',colormatrix_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
end
xlim(handles.axes1,[min(correlationcurves_Chs(:,1,1)) max(correlationcurves_Chs(:,1,1))]);
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')
save('meancorrcurveCC')
for i=1:6
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),correlationcurvesCC_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),meancorrcurveCC(:,i),'LineStyle','-','Color',colormatrixCC_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
end
xlim(handles.axes8,[min(correlationcurvesCC_Chs(:,1,1)) max(correlationcurvesCC_Chs(:,1,1))]);
ylim(handles.axes8,[ylowlimit yuplimit]);
hold (handles.axes8,'off')

% --- Outputs from this function are returned to the command line.
function varargout = Corrselection4Channels_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl_Chs correlationcurves_Chs correlationcurvesCC_Chs ...
sigmascurves_Chs sigmascurvesCC_Chs meancorrcurve meancorrcurveCC ItraceII_Chs...
correlationcurves_Chs_2 correlationcurvesCC_Chs_2 sigmascurves_Chs_2 ...
sigmascurvesCC_Chs_2 yuplimit ylowlimit colormatrix_Chs colormatrixCC_Chs

% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if get(hObject,'Value')<size(correlationcurves_Chs,2)
    correlationcurves_Chs_2=correlationcurves_Chs;
    correlationcurvesCC_Chs_2=correlationcurvesCC_Chs;
    sigmascurves_Chs_2=sigmascurves_Chs;
    sigmascurvesCC_Chs_2=sigmascurvesCC_Chs;

    for kl=2:size(correlationcurves_Chs,2)
        ccindex=1;
        for i=1:4
            if curveincl_Chs(kl-1,i)==0
                correlationcurves_Chs_2(:,kl,i)=NaN;
                sigmascurves_Chs_2(:,kl,i)=NaN;
            end
            if i<4
                for j=i+1:4
                    if curveincl_Chs(kl-1,i)==0 || curveincl_Chs(kl-1,j)==0
                    correlationcurvesCC_Chs_2(:,kl,ccindex)=NaN;
                    sigmascurvesCC_Chs_2(:,kl,ccindex)=NaN;
                    end
                    ccindex=ccindex+1;
                end
            end 
        end
    end

    meancorrcurve=nanmean(correlationcurves_Chs_2(:,2:end,:),2);
    meancorrcurveCC=nanmean(correlationcurvesCC_Chs_2(:,2:end,:),2);

    set(handles.checkbox2,'Value',curveincl_Chs(get(hObject,'Value'),1));
    set(handles.checkbox3,'Value',curveincl_Chs(get(hObject,'Value'),2));
    set(handles.checkbox4,'Value',curveincl_Chs(get(hObject,'Value'),3));
    set(handles.checkbox5,'Value',curveincl_Chs(get(hObject,'Value'),4));

    plot(handles.axes2,ItraceII_Chs(:,1,1),ItraceII_Chs(:,get(hObject,'Value')+1,1))
    hold (handles.axes2,'on')
    plot(handles.axes2,[ItraceII_Chs(1,1,1) ItraceII_Chs(size(ItraceII_Chs,1),1,1) ],[mean(ItraceII_Chs(:,get(hObject,'Value')+1,1),1) mean(ItraceII_Chs(:,get(hObject,'Value')+1,1),1)],'r-')
    ylim(handles.axes2,[min(ItraceII_Chs(:,get(hObject,'Value')+1,1)) max(ItraceII_Chs(:,get(hObject,'Value')+1,1))])
    hold (handles.axes2,'off')

    plot(handles.axes3,ItraceII_Chs(:,1,2),ItraceII_Chs(:,get(hObject,'Value')+1,2))
    hold (handles.axes3,'on')
    plot(handles.axes3,[ItraceII_Chs(1,1,2) ItraceII_Chs(size(ItraceII_Chs,1),1,2) ],[mean(ItraceII_Chs(:,get(hObject,'Value')+1,2),1) mean(ItraceII_Chs(:,get(hObject,'Value')+1,2),1)],'r-')
    ylim(handles.axes3,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2)) max(ItraceII_Chs(:,get(hObject,'Value')+1,2))])
    hold (handles.axes3,'off')

    plot(handles.axes6,ItraceII_Chs(:,1,3),ItraceII_Chs(:,get(hObject,'Value')+1,3))
    hold (handles.axes6,'on')
    plot(handles.axes6,[ItraceII_Chs(1,1,3) ItraceII_Chs(size(ItraceII_Chs,1),1,3)],[mean(ItraceII_Chs(:,get(hObject,'Value')+1,3),1) mean(ItraceII_Chs(:,get(hObject,'Value')+1,3),1)],'r-')
    ylim(handles.axes6,[min(ItraceII_Chs(:,get(hObject,'Value')+1,3)) max(ItraceII_Chs(:,get(hObject,'Value')+1,3))])
    hold (handles.axes6,'off')

    plot(handles.axes7,ItraceII_Chs(:,1,4),ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))
    hold (handles.axes7,'on')
    plot(handles.axes7,[ItraceII_Chs(1,1,4) ItraceII_Chs(size(ItraceII_Chs,1),1,4)],[mean(ItraceII_Chs(:,get(hObject,'Value')+1,4),1) mean(ItraceII_Chs(:,get(hObject,'Value')+1,4),1)],'r-')
    ylim(handles.axes7,[min(ItraceII_Chs(:,get(hObject,'Value')+1,4)) max(ItraceII_Chs(:,get(hObject,'Value')+1,4))])
    hold (handles.axes7,'off')

    yuplimit=max(max(correlationcurves_Chs(:,get(hObject,'Value')+1,:)));
    ylowlimit=min(min(correlationcurves_Chs(:,get(hObject,'Value')+1,:)));

    for i=1:4
        semilogx(handles.axes1,correlationcurves_Chs(:,1,i),correlationcurves_Chs(:,get(hObject,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrix_Chs(i,:))
        hold (handles.axes1,'on')
        semilogx(handles.axes1,correlationcurves_Chs(:,1,i),meancorrcurve(:,i),'LineStyle','-','Color',colormatrix_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrix_Chs(i,:))
        hold (handles.axes1,'on')
    end
    xlim(handles.axes1,[min(correlationcurves_Chs(:,1,1)) max(correlationcurves_Chs(:,1,1))]);
    ylim(handles.axes1,[ylowlimit yuplimit]);
    hold (handles.axes1,'off')
    for i=1:6
        semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),correlationcurvesCC_Chs(:,get(hObject,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrixCC_Chs(i,:))
        hold (handles.axes8,'on')
        semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),meancorrcurveCC(:,i),'LineStyle','-','Color',colormatrixCC_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrixCC_Chs(i,:))
        hold (handles.axes8,'on')
    end
    xlim(handles.axes8,[min(correlationcurvesCC_Chs(:,1,1)) max(correlationcurvesCC_Chs(:,1,1))]);
    ylim(handles.axes8,[ylowlimit yuplimit]);
    hold (handles.axes8,'off')
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
global x0 y0 fixed curveincl_Chs correlationcurves_Chs correlationcurvesCC_Chs ...
sigmascurves_Chs sigmascurvesCC_Chs meancorrcurve meancorrcurveCC ItraceII_Chs...
correlationcurves_Chs_2 correlationcurvesCC_Chs_2 sigmascurves_Chs_2 ...
sigmascurvesCC_Chs_2 yuplimit ylowlimit colormatrix_Chs colormatrixCC_Chs
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

curveincl_Chs(get(handles.slider3,'value'),1)=get(hObject,'Value');
correlationcurves_Chs_2=correlationcurves_Chs;
correlationcurvesCC_Chs_2=correlationcurvesCC_Chs;

sigmascurves_Chs_2=sigmascurves_Chs;
sigmascurvesCC_Chs_2=sigmascurvesCC_Chs;

for kl=2:size(correlationcurves_Chs,2)
    ccindex=1;
    for i=1:4
        if curveincl_Chs(kl-1,i)==0
            correlationcurves_Chs_2(:,kl,i)=NaN;
            sigmascurves_Chs_2(:,kl,i)=NaN;
        end
        if i<4
            for j=i+1:4
                if curveincl_Chs(kl-1,i)==0 || curveincl_Chs(kl-1,j)==0
                correlationcurvesCC_Chs_2(:,kl,ccindex)=NaN;
                sigmascurvesCC_Chs_2(:,kl,ccindex)=NaN;
                end
                ccindex=ccindex+1;
            end
        end
    end
end

meancorrcurve=nanmean(correlationcurves_Chs_2(:,2:end,:),2);
meancorrcurveCC=nanmean(correlationcurvesCC_Chs_2(:,2:end,:),2);
  
set(hObject,'Value',curveincl_Chs(get(handles.slider3,'Value'),1));

plot(handles.axes2,ItraceII_Chs(:,1,1),ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))
hold (handles.axes2,'on')
plot(handles.axes2,[ItraceII_Chs(1,1,1) ItraceII_Chs(size(ItraceII_Chs,1),1,1) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1),1)],'r-')
ylim(handles.axes2,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))])
hold (handles.axes2,'off')

plot(handles.axes3,ItraceII_Chs(:,1,2),ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))
hold (handles.axes3,'on')
plot(handles.axes3,[ItraceII_Chs(1,1,2) ItraceII_Chs(size(ItraceII_Chs,1),1,2) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1)],'r-')
ylim(handles.axes3,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))])
hold (handles.axes3,'off')

plot(handles.axes6,ItraceII_Chs(:,1,3),ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))
hold (handles.axes6,'on')
plot(handles.axes6,[ItraceII_Chs(1,1,3) ItraceII_Chs(size(ItraceII_Chs,1),1,3)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1)],'r-')
ylim(handles.axes6,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))])
hold (handles.axes6,'off')

plot(handles.axes7,ItraceII_Chs(:,1,4),ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))
hold (handles.axes7,'on')
plot(handles.axes7,[ItraceII_Chs(1,1,4) ItraceII_Chs(size(ItraceII_Chs,1),1,4)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1)],'r-')
ylim(handles.axes7,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))])
hold (handles.axes7,'off')
     
yuplimit=max(max(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));
ylowlimit=min(min(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));

for i=1:4
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),correlationcurves_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),meancorrcurve(:,i),'LineStyle','-','Color',colormatrix_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
end
xlim(handles.axes1,[min(correlationcurves_Chs(:,1,1)) max(correlationcurves_Chs(:,1,1))]);
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')
for i=1:6
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),correlationcurvesCC_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),meancorrcurveCC(:,i),'LineStyle','-','Color',colormatrixCC_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
end
xlim(handles.axes8,[min(correlationcurvesCC_Chs(:,1,1)) max(correlationcurvesCC_Chs(:,1,1))]);
ylim(handles.axes8,[ylowlimit yuplimit]);
hold (handles.axes8,'off')

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global goon
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
goon=1;
close Corrselection4Channels


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl_Chs correlationcurves_Chs correlationcurvesCC_Chs ...
sigmascurves_Chs sigmascurvesCC_Chs meancorrcurve meancorrcurveCC ItraceII_Chs...
correlationcurves_Chs_2 correlationcurvesCC_Chs_2 sigmascurves_Chs_2 ...
sigmascurvesCC_Chs_2 yuplimit ylowlimit colormatrix_Chs colormatrixCC_Chs
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3

curveincl_Chs(get(handles.slider3,'value'),2)=get(hObject,'Value');
correlationcurves_Chs_2=correlationcurves_Chs;
correlationcurvesCC_Chs_2=correlationcurvesCC_Chs;

sigmascurves_Chs_2=sigmascurves_Chs;
sigmascurvesCC_Chs_2=sigmascurvesCC_Chs;

for kl=2:size(correlationcurves_Chs,2)
    ccindex=1;
    for i=1:4
        if curveincl_Chs(kl-1,i)==0
            correlationcurves_Chs_2(:,kl,i)=NaN;
            sigmascurves_Chs_2(:,kl,i)=NaN;
        end
        if i<4
            for j=i+1:4
                if curveincl_Chs(kl-1,i)==0 || curveincl_Chs(kl-1,j)==0
                correlationcurvesCC_Chs_2(:,kl,ccindex)=NaN;
                sigmascurvesCC_Chs_2(:,kl,ccindex)=NaN;
                end
                ccindex=ccindex+1;
            end
        end
    end
end

meancorrcurve=nanmean(correlationcurves_Chs_2(:,2:end,:),2);
meancorrcurveCC=nanmean(correlationcurvesCC_Chs_2(:,2:end,:),2);
  
set(hObject,'Value',curveincl_Chs(get(handles.slider3,'Value'),2));

plot(handles.axes2,ItraceII_Chs(:,1,1),ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))
hold (handles.axes2,'on')
plot(handles.axes2,[ItraceII_Chs(1,1,1) ItraceII_Chs(size(ItraceII_Chs,1),1,1) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1),1)],'r-')
ylim(handles.axes2,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))])
hold (handles.axes2,'off')

plot(handles.axes3,ItraceII_Chs(:,1,2),ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))
hold (handles.axes3,'on')
plot(handles.axes3,[ItraceII_Chs(1,1,2) ItraceII_Chs(size(ItraceII_Chs,1),1,2) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1)],'r-')
ylim(handles.axes3,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))])
hold (handles.axes3,'off')

plot(handles.axes6,ItraceII_Chs(:,1,3),ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))
hold (handles.axes6,'on')
plot(handles.axes6,[ItraceII_Chs(1,1,3) ItraceII_Chs(size(ItraceII_Chs,1),1,3)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1)],'r-')
ylim(handles.axes6,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))])
hold (handles.axes6,'off')

plot(handles.axes7,ItraceII_Chs(:,1,4),ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))
hold (handles.axes7,'on')
plot(handles.axes7,[ItraceII_Chs(1,1,4) ItraceII_Chs(size(ItraceII_Chs,1),1,4)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1)],'r-')
ylim(handles.axes7,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))])
hold (handles.axes7,'off')
     
yuplimit=max(max(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));
ylowlimit=min(min(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));

for i=1:4
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),correlationcurves_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),meancorrcurve(:,i),'LineStyle','-','Color',colormatrix_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
end
xlim(handles.axes1,[min(correlationcurves_Chs(:,1,1)) max(correlationcurves_Chs(:,1,1))]);
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')
for i=1:6
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),correlationcurvesCC_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),meancorrcurveCC(:,i),'LineStyle','-','Color',colormatrixCC_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
end
xlim(handles.axes8,[min(correlationcurvesCC_Chs(:,1,1)) max(correlationcurvesCC_Chs(:,1,1))]);
ylim(handles.axes8,[ylowlimit yuplimit]);
hold (handles.axes8,'off')

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl_Chs correlationcurves_Chs correlationcurvesCC_Chs ...
sigmascurves_Chs sigmascurvesCC_Chs meancorrcurve meancorrcurveCC ItraceII_Chs...
correlationcurves_Chs_2 correlationcurvesCC_Chs_2 sigmascurves_Chs_2 ...
sigmascurvesCC_Chs_2 yuplimit ylowlimit colormatrix_Chs colormatrixCC_Chs
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


curveincl_Chs(get(handles.slider3,'value'),3)=get(hObject,'Value');
correlationcurves_Chs_2=correlationcurves_Chs;
correlationcurvesCC_Chs_2=correlationcurvesCC_Chs;

sigmascurves_Chs_2=sigmascurves_Chs;
sigmascurvesCC_Chs_2=sigmascurvesCC_Chs;

for kl=2:size(correlationcurves_Chs,2)
    ccindex=1;
    for i=1:4
        if curveincl_Chs(kl-1,i)==0
            correlationcurves_Chs_2(:,kl,i)=NaN;
            sigmascurves_Chs_2(:,kl,i)=NaN;
        end
        if i<4
            for j=i+1:4
                if curveincl_Chs(kl-1,i)==0 || curveincl_Chs(kl-1,j)==0
                correlationcurvesCC_Chs_2(:,kl,ccindex)=NaN;
                sigmascurvesCC_Chs_2(:,kl,ccindex)=NaN;
                end
                ccindex=ccindex+1;
            end
        end
    end
end

meancorrcurve=nanmean(correlationcurves_Chs_2(:,2:end,:),2);
meancorrcurveCC=nanmean(correlationcurvesCC_Chs_2(:,2:end,:),2);
  
set(hObject,'Value',curveincl_Chs(get(handles.slider3,'Value'),3));

plot(handles.axes2,ItraceII_Chs(:,1,1),ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))
hold (handles.axes2,'on')
plot(handles.axes2,[ItraceII_Chs(1,1,1) ItraceII_Chs(size(ItraceII_Chs,1),1,1) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1),1)],'r-')
ylim(handles.axes2,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))])
hold (handles.axes2,'off')

plot(handles.axes3,ItraceII_Chs(:,1,2),ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))
hold (handles.axes3,'on')
plot(handles.axes3,[ItraceII_Chs(1,1,2) ItraceII_Chs(size(ItraceII_Chs,1),1,2) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1)],'r-')
ylim(handles.axes3,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))])
hold (handles.axes3,'off')

plot(handles.axes6,ItraceII_Chs(:,1,3),ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))
hold (handles.axes6,'on')
plot(handles.axes6,[ItraceII_Chs(1,1,3) ItraceII_Chs(size(ItraceII_Chs,1),1,3)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1)],'r-')
ylim(handles.axes6,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))])
hold (handles.axes6,'off')

plot(handles.axes7,ItraceII_Chs(:,1,4),ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))
hold (handles.axes7,'on')
plot(handles.axes7,[ItraceII_Chs(1,1,4) ItraceII_Chs(size(ItraceII_Chs,1),1,4)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1)],'r-')
ylim(handles.axes7,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))])
hold (handles.axes7,'off')
     
yuplimit=max(max(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));
ylowlimit=min(min(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));

for i=1:4
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),correlationcurves_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),meancorrcurve(:,i),'LineStyle','-','Color',colormatrix_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
end
xlim(handles.axes1,[min(correlationcurves_Chs(:,1,1)) max(correlationcurves_Chs(:,1,1))]);
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')
for i=1:6
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),correlationcurvesCC_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),meancorrcurveCC(:,i),'LineStyle','-','Color',colormatrixCC_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
end
xlim(handles.axes8,[min(correlationcurvesCC_Chs(:,1,1)) max(correlationcurvesCC_Chs(:,1,1))]);
ylim(handles.axes8,[ylowlimit yuplimit]);
hold (handles.axes8,'off')


function edit1_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl_Chs correlationcurves_Chs correlationcurvesCC_Chs ...
sigmascurves_Chs sigmascurvesCC_Chs meancorrcurve meancorrcurveCC ItraceII_Chs...
correlationcurves_Chs_2 correlationcurvesCC_Chs_2 sigmascurves_Chs_2 ...
sigmascurvesCC_Chs_2 yuplimit ylowlimit colormatrix_Chs colormatrixCC_Chs
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
yuplimit=str2double(get(hObject,'String'));
ylowlimit=-str2double(get(hObject,'String'));

for i=1:4
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),correlationcurves_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),meancorrcurve(:,i),'LineStyle','-','Color',colormatrix_Chs(i,:),'Marker','d','MarkerFaceColor',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
end
xlim(handles.axes1,[min(correlationcurves_Chs(:,1,1)) max(correlationcurves_Chs(:,1,1))]);
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')
for i=1:6
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),correlationcurvesCC_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),meancorrcurveCC(:,i),'LineStyle','-','Color',colormatrixCC_Chs(i,:),'Marker','d','MarkerFaceColor',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
end
xlim(handles.axes8,[min(correlationcurvesCC_Chs(:,1,1)) max(correlationcurvesCC_Chs(:,1,1))]);
ylim(handles.axes8,[ylowlimit yuplimit]);
hold (handles.axes8,'off')

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl_Chs correlationcurves_Chs correlationcurvesCC_Chs ...
sigmascurves_Chs sigmascurvesCC_Chs meancorrcurve meancorrcurveCC ItraceII_Chs...
correlationcurves_Chs_2 correlationcurvesCC_Chs_2 sigmascurves_Chs_2 ...
sigmascurvesCC_Chs_2 yuplimit ylowlimit colormatrix_Chs colormatrixCC_Chs
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5

curveincl_Chs(get(handles.slider3,'value'),4)=get(hObject,'Value');
correlationcurves_Chs_2=correlationcurves_Chs;
correlationcurvesCC_Chs_2=correlationcurvesCC_Chs;

sigmascurves_Chs_2=sigmascurves_Chs;
sigmascurvesCC_Chs_2=sigmascurvesCC_Chs;

for kl=2:size(correlationcurves_Chs,2)
    ccindex=1;
    for i=1:4
        if curveincl_Chs(kl-1,i)==0
            correlationcurves_Chs_2(:,kl,i)=NaN;
            sigmascurves_Chs_2(:,kl,i)=NaN;
        end
        for j=i+1:4
            if curveincl_Chs(kl-1,i)==0 || curveincl_Chs(kl-1,j)==0
            correlationcurvesCC_Chs_2(:,kl,ccindex)=NaN;
            sigmascurvesCC_Chs_2(:,kl,ccindex)=NaN;
            end
            ccindex=ccindex+1;
        end
            
    end
end

meancorrcurve=nanmean(correlationcurves_Chs_2(:,2:end,:),2);
meancorrcurveCC=nanmean(correlationcurvesCC_Chs_2(:,2:end,:),2);
  
set(hObject,'Value',curveincl_Chs(get(handles.slider3,'Value'),4));

plot(handles.axes2,ItraceII_Chs(:,1,1),ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))
hold (handles.axes2,'on')
plot(handles.axes2,[ItraceII_Chs(1,1,1) ItraceII_Chs(size(ItraceII_Chs,1),1,1) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1),1)],'r-')
ylim(handles.axes2,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,1))])
hold (handles.axes2,'off')

plot(handles.axes3,ItraceII_Chs(:,1,2),ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))
hold (handles.axes3,'on')
plot(handles.axes3,[ItraceII_Chs(1,1,2) ItraceII_Chs(size(ItraceII_Chs,1),1,2) ],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2),1)],'r-')
ylim(handles.axes3,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,2))])
hold (handles.axes3,'off')

plot(handles.axes6,ItraceII_Chs(:,1,3),ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))
hold (handles.axes6,'on')
plot(handles.axes6,[ItraceII_Chs(1,1,3) ItraceII_Chs(size(ItraceII_Chs,1),1,3)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3),1)],'r-')
ylim(handles.axes6,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,3))])
hold (handles.axes6,'off')

plot(handles.axes7,ItraceII_Chs(:,1,4),ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))
hold (handles.axes7,'on')
plot(handles.axes7,[ItraceII_Chs(1,1,4) ItraceII_Chs(size(ItraceII_Chs,1),1,4)],[mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1) mean(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4),1)],'r-')
ylim(handles.axes7,[min(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4)) max(ItraceII_Chs(:,get(handles.slider3,'Value')+1,4))])
hold (handles.axes7,'off')
     
yuplimit=max(max(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));
ylowlimit=min(min(correlationcurves_Chs(:,get(handles.slider3,'Value')+1,:)));

for i=1:4
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),correlationcurves_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
    semilogx(handles.axes1,correlationcurves_Chs(:,1,i),meancorrcurve(:,i),'LineStyle','-','Color',colormatrix_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrix_Chs(i,:))
    hold (handles.axes1,'on')
end
xlim(handles.axes1,[min(correlationcurves_Chs(:,1,1)) max(correlationcurves_Chs(:,1,1))]);
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')
for i=1:6
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),correlationcurvesCC_Chs(:,get(handles.slider3,'Value')+1,i),'LineStyle','none','Marker','d','Color',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
    semilogx(handles.axes8,correlationcurvesCC_Chs(:,1,i),meancorrcurveCC(:,i),'LineStyle','-','Color',colormatrixCC_Chs(i,:),'Marker','.','MarkerFaceColor',colormatrixCC_Chs(i,:))
    hold (handles.axes8,'on')
end
xlim(handles.axes8,[min(correlationcurvesCC_Chs(:,1,1)) max(correlationcurvesCC_Chs(:,1,1))]);
ylim(handles.axes8,[ylowlimit yuplimit]);
hold (handles.axes8,'off')
