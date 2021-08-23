function varargout = Corrselection2Channelsindividuellpreview(varargin)
% CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW MATLAB code for Corrselection2Channelsindividuellpreview.fig
%      CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW, by itself, creates a new CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW or raises the existing
%      singleton*.
%
%      H = CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW returns the handle to a new CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW or the handle to
%      the existing singleton*.
%
%      CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW.M with the given input arguments.
%
%      CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW('Property','Value',...) creates a new CORRSELECTION2CHANNELSINDIVIDUELLPREVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Corrselection2Channelsindividuellpreview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Corrselection2Channelsindividuellpreview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Corrselection2Channelsindividuellpreview

% Last Modified by GUIDE v2.5 16-Mar-2020 19:13:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Corrselection2Channelsindividuellpreview_OpeningFcn, ...
                   'gui_OutputFcn',  @Corrselection2Channelsindividuellpreview_OutputFcn, ...
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

% --- Executes just before Corrselection2Channelsindividuellpreview is made visible.
function Corrselection2Channelsindividuellpreview_OpeningFcn(hObject, eventdata, handles, varargin)
global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves corfit meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2 yuplimit ylowlimit
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Corrselection2Channelsindividuellpreview (see VARARGIN)

% Choose default command line output for Corrselection2Channelsindividuellpreview
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.slider3,'Max',size(correlationcurvesCh1,2)-1);
set(handles.slider3,'SliderStep',[1/(size(correlationcurvesCh1,2)-2) 1/(size(correlationcurvesCh1,2)-2)]);
correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

sigmas1curves2=sigmas1curves;
sigmas2curves2=sigmas2curves;
sigmascccurves2=sigmascccurves;


for kl=2:size(correlationcurvesCh1,2)
    if curveincl1(kl-1)==0
        correlationcurves2Ch1(:,kl)=NaN;
        correlationcurves2ChCC(:,kl)=NaN;
        sigmas1curves2(:,kl-1)=NaN;
        sigmascccurves2(:,kl-1)=NaN;
    end
    if curveincl2(kl-1)==0
    correlationcurves2Ch2(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
    sigmas2curves2(:,kl-1)=NaN;
    sigmascccurves2(:,kl-1)=NaN;
    end
end
meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);

sigmas1curvesmean=nanmean(sigmas1curves2,2);
sigmas2curvesmean=nanmean(sigmas2curves2,2);
sigmascccurvesmean=nanmean(sigmascccurves2,2);

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

semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+',correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+',...
correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
hold (handles.axes1,'on')

semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])

yuplimit=max(max([correlationcurvesCh1(:,get(handles.slider3,'Value')+1) correlationcurvesCh2(:,get(handles.slider3,'Value')+1)]));
ylowlimit=min(min([correlationcurvesCh1(:,get(handles.slider3,'Value')+1) correlationcurvesCh2(:,get(handles.slider3,'Value')+1)]));
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')

tcorr1segfinalpre=correlationcurvesCh1(:,1);
tcorr2segfinalpre=correlationcurvesCh2(:,1);
tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
lbseg=tcorr1segfinalpre(1);
ubseg=tcorr1segfinalpre(end);
lbccseg=tcorrccsegfinalfitpre(1);
ubccseg=tcorrccsegfinalfitpre(end);

fcorr1segfinalpre=meancorrcurve(:,1);
fcorr2segfinalpre=meancorrcurve(:,2);
fcrosscorrsegfinalfitpre=meancorrcurveCC;

weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);

% Possibility here to fit average of selected segments and show fit
% functions in seperate window (handles.axes4/ handles.axes5) 
% -------> in progress, has not been finishes yet

    % lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
    % flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);
    % 
    % %x0=[100,0.5, S]; % Set in header or GUI manually later!
    % %fixed=[false false true];
    % fixedcc=[false false true false];
    % %flscrossfitfunc=lsautofitfunc;
    % %flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
    % %y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
    % %y0=[2000,0.5,S];
    % 
    % % [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
    % % [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
    % % [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);
    % 
    % [N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
    % [N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiff(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);
    % [Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiff(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
    % 
    % % h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
    % % positionvector1=[0.1 0.35 0.8 0.55];
    % % positionvector2=[0.1 0.1 0.8 0.15];
    % % subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
    % % hold on
    % % subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
    % % xlabel('time')
    % % ylabel('autocorrelation')
    % 
    % %hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
    % %positionvector1=[0.1 0.35 0.8 0.55];
    % %positionvector2=[0.1 0.1 0.8 0.15];
    % semilogx(handles.axes5,tcorr1segfinalpre,fitcurve1segfinalpre,'-g',tcorr2segfinalpre,fitcurve2segfinalpre,'-r',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-b')
    % hold (handles.axes5,'on')
    % semilogx(handles.axes5,tcorr1segfinalpre,fcorr1segfinalpre,'gs',tcorr2segfinalpre,fcorr2segfinalpre,'rd',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'bx')
    % hold (handles.axes5,'off')


% UIWAIT makes Corrselection2Channelsindividuellpreview wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Corrselection2Channelsindividuellpreview_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves corfit meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2 yuplimit ylowlimit
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


if get(hObject,'Value')<size(correlationcurvesCh1,2) 
    correlationcurves2Ch1=correlationcurvesCh1;
    correlationcurves2Ch2=correlationcurvesCh2;
    correlationcurves2ChCC=correlationcurvesChCC;

    sigmas1curves2=sigmas1curves;
    sigmas2curves2=sigmas2curves;
    sigmascccurves2=sigmascccurves;

    for kl=2:size(correlationcurvesCh1,2)
        if curveincl1(kl-1)==0
            correlationcurves2Ch1(:,kl)=NaN;
            correlationcurves2ChCC(:,kl)=NaN;
            sigmas1curves2(:,kl-1)=NaN;
            sigmascccurves2(:,kl-1)=NaN;
        end
        if curveincl2(kl-1)==0
        correlationcurves2Ch2(:,kl)=NaN;
        correlationcurves2ChCC(:,kl)=NaN;
        sigmas2curves2(:,kl-1)=NaN;
        sigmascccurves2(:,kl-1)=NaN;
        end
    end
    meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
    meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
    meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);

    sigmas1curvesmean=nanmean(sigmas1curves2,2);
    sigmas2curvesmean=nanmean(sigmas2curves2,2);
    sigmascccurvesmean=nanmean(sigmascccurves2,2);

    set(handles.checkbox2,'Value',curveincl1(get(hObject,'Value'),1));
    set(handles.checkbox3,'Value',curveincl2(get(hObject,'Value'),1));
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


    semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(hObject,'Value')+1),'g+',correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(hObject,'Value')+1),'r+',...
     correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(hObject,'Value')+1),'b+')
    hold (handles.axes1,'on')

    semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
    xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
    yuplimit=max(max([correlationcurvesCh1(:,get(handles.slider3,'Value')+1) correlationcurvesCh2(:,get(handles.slider3,'Value')+1)]));
    ylowlimit=min(min([correlationcurvesCh1(:,get(handles.slider3,'Value')+1) correlationcurvesCh2(:,get(handles.slider3,'Value')+1)]));
    ylim(handles.axes1,[ylowlimit yuplimit]);
    hold (handles.axes1,'off')

    tcorr1segfinalpre=correlationcurvesCh1(:,1);
    tcorr2segfinalpre=correlationcurvesCh2(:,1);
    tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
    lbseg=tcorr1segfinalpre(1);
    ubseg=tcorr1segfinalpre(end);
    lbccseg=tcorrccsegfinalfitpre(1);
    ubccseg=tcorrccsegfinalfitpre(end);

    fcorr1segfinalpre=meancorrcurve(:,1);
    fcorr2segfinalpre=meancorrcurve(:,2);
    fcrosscorrsegfinalfitpre=meancorrcurveCC;

    weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
    weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
    weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);

   
    % Possibility here to fit average of selected segments and show fit
    % functions in seperate window (handles.axes4/ handles.axes5) 
    % -------> in progress, has not been finishes yet

        % lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
        % flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);
        % 
        % %x0=[100,0.5, S]; % Set in header or GUI manually later!
        % %fixed=[false false true];
        % fixedcc=[false false true false];
        % %flscrossfitfunc=lsautofitfunc;
        % %flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
        % %y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
        % %y0=[2000,0.5,S];
        % 
        % % For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
        % % [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
        % % [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
        % % [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);
        % 
        % [N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
        % [N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiff(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);
        % [Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiff(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
        % 
        % % h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
        % % positionvector1=[0.1 0.35 0.8 0.55];
        % % positionvector2=[0.1 0.1 0.8 0.15];
        % % subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
        % % hold on
        % % subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
        % % xlabel('time')
        % % ylabel('autocorrelation')
        % 
        % %hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
        % %positionvector1=[0.1 0.35 0.8 0.55];
        % %positionvector2=[0.1 0.1 0.8 0.15];
        % semilogx(handles.axes5,tcorr1segfinalpre,fitcurve1segfinalpre,'-g',tcorr2segfinalpre,fitcurve2segfinalpre,'-r',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-b')
        % hold (handles.axes5,'on')
        % semilogx(handles.axes5,tcorr1segfinalpre,fcorr1segfinalpre,'gs',tcorr2segfinalpre,fcorr2segfinalpre,'rd',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'bx')
        % hold (handles.axes5,'off')
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
global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves corfit meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2 yuplimit ylowlimit

% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

curveincl1(get(handles.slider3,'value'),1)=get(hObject,'Value');
correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

sigmas1curves2=sigmas1curves;
sigmas2curves2=sigmas2curves;
sigmascccurves2=sigmascccurves;

for kl=2:size(correlationcurvesCh1,2)
    if curveincl1(kl-1)==0
        correlationcurves2Ch1(:,kl)=NaN;
        correlationcurves2ChCC(:,kl)=NaN;
        sigmas1curves2(:,kl-1)=NaN;
        sigmascccurves2(:,kl-1)=NaN;
    end
    if curveincl2(kl-1)==0
    correlationcurves2Ch2(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
    sigmas2curves2(:,kl-1)=NaN;
    sigmascccurves2(:,kl-1)=NaN;
    end
end
meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);

sigmas1curvesmean=nanmean(sigmas1curves2,2);
sigmas2curvesmean=nanmean(sigmas2curves2,2);
sigmascccurvesmean=nanmean(sigmascccurves2,2);
  
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

semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+',correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+',...
correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
hold (handles.axes1,'on')

semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')

tcorr1segfinalpre=correlationcurvesCh1(:,1);
tcorr2segfinalpre=correlationcurvesCh2(:,1);
tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
lbseg=tcorr1segfinalpre(1);
ubseg=tcorr1segfinalpre(end);
lbccseg=tcorrccsegfinalfitpre(1);
ubccseg=tcorrccsegfinalfitpre(end);

fcorr1segfinalpre=meancorrcurve(:,1);
fcorr2segfinalpre=meancorrcurve(:,2);
fcrosscorrsegfinalfitpre=meancorrcurveCC;

weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);


% Possibility here to fit average of selected segments and show fit
% functions in seperate window (handles.axes4/ handles.axes5) 
% -------> in progress, has not been finishes yet

    % lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
    % flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);
    % 
    % %x0=[100,0.5, S]; % Set in header or GUI manually later!
    % %fixed=[false false true];
    % fixedcc=[false false true false];
    % %flscrossfitfunc=lsautofitfunc;
    % %flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
    % %y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
    % %y0=[2000,0.5,S];
    % 
    % % For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
    % % [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
    % % [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
    % % %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
    % % [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);
    % 
    % [N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
    % [N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiff(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);
    % %[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);
    % [Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiff(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
    % 
    % % h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
    % % positionvector1=[0.1 0.35 0.8 0.55];
    % % positionvector2=[0.1 0.1 0.8 0.15];
    % % subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
    % % hold on
    % % subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
    % % xlabel('time')
    % % ylabel('autocorrelation')
    % 
    % %hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
    % %positionvector1=[0.1 0.35 0.8 0.55];
    % %positionvector2=[0.1 0.1 0.8 0.15];
    % semilogx(handles.axes5,tcorr1segfinalpre,fitcurve1segfinalpre,'-g',tcorr2segfinalpre,fitcurve2segfinalpre,'-r',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-b')
    % hold (handles.axes5,'on')
    % semilogx(handles.axes5,tcorr1segfinalpre,fcorr1segfinalpre,'gs',tcorr2segfinalpre,fcorr2segfinalpre,'rd',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'bx')
    % hold (handles.axes5,'off')

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
global goon
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
goon=1;
close Corrselection2Channelsindividuellpreview


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves corfit meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2 yuplimit ylowlimit

% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

curveincl2(get(handles.slider3,'value'),1)=get(hObject,'Value');
correlationcurves2Ch1=correlationcurvesCh1;
correlationcurves2Ch2=correlationcurvesCh2;
correlationcurves2ChCC=correlationcurvesChCC;

sigmas1curves2=sigmas1curves;
sigmas2curves2=sigmas2curves;
sigmascccurves2=sigmascccurves;

for kl=2:size(correlationcurvesCh1,2)
    if curveincl1(kl-1)==0
        correlationcurves2Ch1(:,kl)=NaN;
        correlationcurves2ChCC(:,kl)=NaN;
        sigmas1curves2(:,kl-1)=NaN;
        sigmascccurves2(:,kl-1)=NaN;
    end
    if curveincl2(kl-1)==0
    correlationcurves2Ch2(:,kl)=NaN;
    correlationcurves2ChCC(:,kl)=NaN;
    sigmas2curves2(:,kl-1)=NaN;
    sigmascccurves2(:,kl-1)=NaN;
    end
end
meancorrcurve(:,1)=nanmean(correlationcurves2Ch1(:,2:end),2);
meancorrcurve(:,2)=nanmean(correlationcurves2Ch2(:,2:end),2);
meancorrcurveCC=nanmean(correlationcurves2ChCC(:,2:end),2);

sigmas1curvesmean=nanmean(sigmas1curves2,2);
sigmas2curvesmean=nanmean(sigmas2curves2,2);
sigmascccurvesmean=nanmean(sigmascccurves2,2);
  
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

semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+',correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+',...
correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
hold (handles.axes1,'on')

semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')

tcorr1segfinalpre=correlationcurvesCh1(:,1);
tcorr2segfinalpre=correlationcurvesCh2(:,1);
tcorrccsegfinalfitpre=correlationcurvesChCC(:,1);
lbseg=tcorr1segfinalpre(1);
ubseg=tcorr1segfinalpre(end);
lbccseg=tcorrccsegfinalfitpre(1);
ubccseg=tcorrccsegfinalfitpre(end);

fcorr1segfinalpre=meancorrcurve(:,1);
fcorr2segfinalpre=meancorrcurve(:,2);
fcrosscorrsegfinalfitpre=meancorrcurveCC;

weights1segfinalpre=abs(fcorr1segfinalpre./sigmas1curvesmean);
weights2segfinalpre=abs(fcorr2segfinalpre./sigmas2curvesmean);
weightsccsegfinalfitpre=abs(fcrosscorrsegfinalfitpre./sigmascccurvesmean);


% Possibility here to fit average of selected segments and show fit
% functions in seperate window (handles.axes4/ handles.axes5) 
% -------> in progress, has not been finishes yet

    % lsautofitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
    % flscrossfitfuncpre=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);
    % 
    % %x0=[100,0.5, S]; % Set in header or GUI manually later!
    % %fixed=[false false true];
    % fixedcc=[false false true false];
    % %flscrossfitfunc=lsautofitfunc;
    % %flscrossfitfunc=@(y,tcc)1/y(1)*((1+4*y(2)*tcc./y(3)^2).^-0.5).*((1+4*y(2)*tcc./(y(3)*y(4))^2).^-0.5).*exp(-d^2./(y(3)^2+4*y(2).*tcc));
    % %y0=[5000,100,S,10^-5]; %[N0,tau0,S], Set in header or GUI manually later!
    % %y0=[2000,0.5,S];
    % 
    % % For now: separate fit. Later: Global fit of all 3 curves (2x auto & cc)
    % % [N1finalpre,taud1finalpre,Sfit1finalpre,CI1finalpre,fitcurve1finalpre,residuals1finalpre] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfuncpre,x0,fixed);
    % % [N2finalpre,taud2finalpre,Sfit2finalpre,CI2finalpre,fitcurve2finalpre,residuals2finalpre] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfuncpre,x0,fixed);
    % % %[Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiffconst(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,flscrossfitfunc,y0,fixedcc);
    % % [Nccfinalpre,tauccfinalpre,Sfitccfinalpre,CIccfinalpre,fitcurveccfinalpre,residualsccfinalpre] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfuncpre,y0,fixed);
    % 
    % [N1segfinalpre,taud1segfinalpre,Sfit1segfinalpre,CI1segfinalpre,fitcurve1segfinalpre,residuals1segfinalpre] = autocorrfit2Ddiff(tcorr1segfinalpre,fcorr1segfinalpre,lbseg,ubseg,weights1segfinalpre,lsautofitfuncpre,x0,fixed);
    % [N2segfinalpre,taud2segfinalpre,Sfit2segfinalpre,CI2segfinalpre,fitcurve2segfinalpre,residuals2segfinalpre] = autocorrfit2Ddiff(tcorr2segfinalpre,fcorr2segfinalpre,lbseg,ubseg,weights2segfinalpre,lsautofitfuncpre,x0,fixed);
    % %[Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiffconst(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',flscrossfitfunc,y0,fixedcc);
    % [Nccsegfinalpre,tauccsegfinalpre,Sfitccsegfinalpre,CIccsegfinalpre,fitcurveccsegfinalpre,residualsccsegfinalpre] =autocorrfit2Ddiff(tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,lbccseg,ubccseg,weightsccsegfinalfitpre,lsautofitfuncpre,y0,fixed);
    % 
    % % h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
    % % positionvector1=[0.1 0.35 0.8 0.55];
    % % positionvector2=[0.1 0.1 0.8 0.15];
    % % subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1finalpre,'-g',tcorr2fit,fitcurve2finalpre,'-r',tcorrccfit,fitcurveccfinalpre,'-b')
    % % hold on
    % % subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
    % % xlabel('time')
    % % ylabel('autocorrelation')
    % 
    % %hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
    % %positionvector1=[0.1 0.35 0.8 0.55];
    % %positionvector2=[0.1 0.1 0.8 0.15];
    % semilogx(handles.axes5,tcorr1segfinalpre,fitcurve1segfinalpre,'-g',tcorr2segfinalpre,fitcurve2segfinalpre,'-r',tcorrccsegfinalfitpre,fitcurveccsegfinalpre,'-b')
    % hold (handles.axes5,'on')
    % semilogx(handles.axes5,tcorr1segfinalpre,fcorr1segfinalpre,'gs',tcorr2segfinalpre,fcorr2segfinalpre,'rd',tcorrccsegfinalfitpre,fcrosscorrsegfinalfitpre,'bx')
    % hold (handles.axes5,'off')

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


function edit1_Callback(hObject, eventdata, handles)
global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves corfit meancorrcurve meancorrcurveCC ItraceIICh1 ItraceIICh2  correlationcurves2Ch1 correlationcurves2Ch2 correlationcurves2ChCC sigmas1curves2 sigmas2curves2 sigmascccurves2 yuplimit ylowlimit
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

yuplimit=str2double(get(hObject,'String'));
ylowlimit=-str2double(get(hObject,'String'));

semilogx(handles.axes1,correlationcurvesCh1(:,1),correlationcurvesCh1(:,get(handles.slider3,'Value')+1),'g+',correlationcurvesCh2(:,1),correlationcurvesCh2(:,get(handles.slider3,'Value')+1),'r+',...
 correlationcurvesChCC(:,1),correlationcurvesChCC(:,get(handles.slider3,'Value')+1),'b+')
hold (handles.axes1,'on')

semilogx(handles.axes1,correlationcurvesCh1(:,1),meancorrcurve(:,1),'g-',correlationcurvesCh2(:,1),meancorrcurve(:,2),'r-',correlationcurvesChCC(:,1),meancorrcurveCC,'b-')
xlim(handles.axes1,[min(correlationcurvesCh1(:,1)) max(correlationcurvesCh1(:,1))])
ylim(handles.axes1,[ylowlimit yuplimit]);
hold (handles.axes1,'off')


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
