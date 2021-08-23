clear all
close all
fclose all;

global movingaveragewindowsize blocksize spatialfilter depletioncorrection membranewidth

% Selection of GUIs and GUI parameters:
ACFselection=1; % 1 divide timetrace in ACFsegmentnumb segments and select or discard segment ACFs
intensityselection=1; % 1(default): GUI to adjust exponential fit for bleaching correction (by removing short intenisty segments)
if ACFselection==1
    global x0 y0 fixed curveincl_Chs correlationcurves_Chs correlationcurvesCC_Chs ...
    sigmascurves_Chs sigmascurvesCC_Chs ItraceII_Chs...
    colormatrix_Chs colormatrixCC_Chs goon
end
if intensityselection==1
    global segmentsincl1 segmentsincl2 segmentsincl3 segmentsincl4 Ifull1 Ifull2 Ifull3 Ifull4 timeline1binned timeline2binned timeline3binned timeline4binned numberofsegments segmentlength correctionfitseg1 correctionfitseg2 correctionfitseg3 correctionfitseg4 ff1 ff2 ff3 ff4
end
ACFsegmentnumb=20; % Number of segments for segment-wise analysis
numberofsegments=20; % Number of intensity segments for exponential fit adjustment
movingaveragewindowsize=1000; % Number of blocks for lateral alignment of lines in kymograph
blocksize=movingaveragewindowsize;

% Acquisition parameters:
specchannelnumb=23; % specify number of spectral bins, e.g. 23
spectralchannels=[495 504 513 522 531 540 548 557 566 575 584 593 602 611 620 629 637 646 655 664 673 682 691]; % Spectral bins
scantime_lines=403.20e-06; % Scan time (time to scan one line and move back the scanner). Minimum scan time on Zeiss LSM780 is 472.73e-06
pixeltime=1.23e-06; % Time to scan one pixel
timebreak=0; % Intervall in between subsequent scans (e.g. line_break_line_break...), default 0
scantime=scantime_lines+timebreak;
S=6.71; % Structure parameter as obtained from pFCS calibration of focal volume

% Analysis parameters
spatialfilter=2.5; % Factor of stdev above and below mean mean membrane position (to select membrane pixels after Gaussian fitting)
depletioncorrection=1; % Apply bleaching correction? Default 1
binningwindow=100; % Intensity binning factor for binned time series (for display purposes only)
binningwindow_GUI=25; % Intensity binning factor for binned time series segments in GUI (for display purposes only)
backgroundcorrection=0; %1: background will be subtracted (as selected by rectangular ROI), 0: no subtraction
binning=0; % Timebinning?
binningfactor=1; % Binning factor for time binning (binningfactor=n -> 2^n binning in time)

% Fit parameters:
fixed=[false false true]; % Fix fit parameters (N, tau, S)? No: false, Yes: true
x0=[500,0.01,S]; % Initial conditions for ACF fits (N0,taud0,S)
y0=[5000,1,S]; % Initial conditions for CCF fits (N0,taud0,S) -> if no binding is expected, high N0 and taud0 allow the fit to converge in all cases

% Import of .tif files:
path= uigetdir; 
files=dir([path '/*.tif']);
[namedata,remain]=strtok(files(1).name,'.');
name1file=files(1).name
inputfilename=name1file;
linearray1=double(imread([path '/' name1file]));

scrsz=   get(0,'ScreenSize');

% Select directory in which to save results:
path2=uigetdir('D:\Shared\Valentin\LSM Zeiss Data');

% Plot kymograph (first three spectral bins):
figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Line Fluorescence')
linearraypreview=sum(linearray1,3);
imagesc(linearraypreview)
xlabel('pixelposition')
ylabel('time (linenumber)')
rect=getrect;
xlb=round(rect(1));
xub=round(rect(1)+rect(3));
linearray1=linearray1(:,xlb:xub,:);
speclinearray=zeros(size(linearray1,1),size(linearray1,2),specchannelnumb);
speclinearray(:,:,1:3)=linearray1;
clear linearray1 linearraypreview

% Continue import:
for i=2:size(files,1)-1
    namefile=files(i).name
    linearray=double(imread([path '/' namefile]));
    speclinearray(:,:,(i-1)*3+1:i*3)=linearray(:,xlb:xub,:);
end
namelastfile=files(size(files,1)).name
linearraylast=double(imread([path '/' namelastfile]));
%speclinearray(:,:,end-1:end)=linearraylast(:,xlb:xub,1:2);
rest=mod(specchannelnumb,3);
if rest==0
    speclinearray(:,:,end-2:end)=linearraylast(:,xlb:xub,:);
else
    speclinearray(:,:,end-rest+1:end)=linearraylast(:,xlb:xub,1:rest);
end
clear linearraylast
fclose all;

% Check whether linearray can be divided in blocks: 
if mod(length(linearray(:,1)),movingaveragewindowsize)>0 || mod(movingaveragewindowsize,2)==1;
    fprintf('Window Size not valid!\n')
    return
end

% Figure: Fluorescence of all scanned lines:    
figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Line Fluorescence')
linearray=sum(speclinearray,3);
imagesc(linearray)
xlabel('pixelposition')
ylabel('time (linenumber)')
    
% Polygonal Selection:
mask=roipoly;
linearraymasked = mask.*linearray;
speclinearraymasked=speclinearray;
for i=1:size(speclinearray,3)
    i
    speclinearraymasked(:,:,i)=mask.*speclinearray(:,:,i);
end
    
% Background correction:
if backgroundcorrection==1
    specbackgroundmasked=speclinearraymasked;
    mask_bgr=roipoly;
    for i=1:size(speclinearray,3)
    specbackgroundmasked(:,:,i)=mask_bgr.*speclinearray(:,:,i);
    end
end

clear speclinearray linearray linedata
fclose all;
    
% Selected kympgraph pixels:
figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Selected Lines')
imagesc(linearraymasked)

% Lateral alignment:
if backgroundcorrection==1
    [speclinefluorescenceseries,speclinearrayalignedba]=speclinealignfuncbr(speclinearraymasked,specbackgroundmasked);
    clear specbackgroundmasked
else
    [speclinefluorescenceseries,speclinearrayalignedba]=speclinealignfunc(speclinearraymasked);
end
    
% Aligned kymograph:    
figure('OuterPosition',[scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Line Fluorescence (aligned)')
linearrayalignedba=sum(speclinearrayalignedba,3);
imagesc(linearrayalignedba);
xlabel('pixelposition (a.u.)')
ylabel('time (linenumber)') 

% Figure binned fluorescence time series (before decomposition):
timeline=1:1:size(speclinefluorescenceseries,1);
timeline=timeline'*scantime;
linefluorescenceseries=sum(speclinefluorescenceseries,3);
f=fit(timeline,linefluorescenceseries,'exp2');
lastbinlength=mod(length(linefluorescenceseries),binningwindow)+binningwindow;
linefluorescenceseriesbinned=mean(reshape(linefluorescenceseries(1:end-lastbinlength),binningwindow,length(linefluorescenceseries(1:end-lastbinlength))/binningwindow),1);
linefluorescenceseriesbinned=[linefluorescenceseriesbinned mean(linefluorescenceseries(end-lastbinlength+1:end))];
correctionfit=f.a.*exp(f.b*timeline)+f.c.*exp(f.d*timeline); 
correctionfitbinned=mean(reshape(correctionfit(1:end-lastbinlength),binningwindow,length(correctionfit(1:end-lastbinlength))/binningwindow),1);
correctionfitbinned=[correctionfitbinned mean(correctionfit(end-lastbinlength+1:end))];
linefluorescencebinnedresiduals=linefluorescenceseriesbinned-correctionfitbinned;

hI=figure('OuterPosition',[scrsz(3) scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2],'Name','Binned Membrane Fluorescence time series');
timeline_binned=1:1:length(linefluorescenceseriesbinned);
timeline_binned=timeline_binned*binningwindow*scantime;
subplot(2,1,1),plot(timeline_binned,linefluorescenceseriesbinned)
hold on
plot(f)
xlabel('time')
ylabel('line fluorescence')
subplot(2,1,2),plot(timeline_binned,linefluorescencebinnedresiduals)
hold on
plot(timeline_binned,0*linefluorescencebinnedresiduals)
xlabel('time')
ylabel('line fluorescence residuals')
    
% Calculation of spectrum:     
emissionspec=squeeze(sum(speclinefluorescenceseries,1)./sum(sum(speclinefluorescenceseries,1),3));
emspec=figure('Name','Spectrum');
plot(spectralchannels,emissionspec)
xlabel('Channel [nm]')
ylabel('Norm.emission')
savefig(emspec,[path2 '\' inputfilename(3:end-4) '_Flpectrum.fig'])
specoutput=zeros(specchannelnumb,2);
specoutput(:,1)=spectralchannels';
specoutput(:,2)=emissionspec';

fidspec=fopen([path2 '\' inputfilename(3:end-4) '_flspectrum.txt'],'a'); 
fprintf(fidspec,'%e\t %e\n',specoutput');

clear speclinearraymasked
clear speclinearrayalignedba

% Spectral decomposition based on spectral filters: 
[inputfilenamespecpatterns, path]=uigetfile('*.txt'); % Load spectral_patterns.txt
spectral_patterns=load([path '/' inputfilenamespecpatterns]);
spectral_patterns=spectral_patterns';
speclinefluorescenceseries=squeeze(speclinefluorescenceseries);
scantime_raw=scantime;
if binning==0
    linefluorescenceseries=linefluorescenceseries';
    linefluorescenceseries_GFP=zeros(size(linefluorescenceseries));
    linefluorescenceseries_YFP=zeros(size(linefluorescenceseries));
    linefluorescenceseries_Ch=zeros(size(linefluorescenceseries));
    linefluorescenceseries_A=zeros(size(linefluorescenceseries));
    Iinv=zeros(1,specchannelnumb);
    Ittt=mean(speclinefluorescenceseries,1);
    Iinv(Ittt>0)=1./Ittt(Ittt>0);
    D=diag(Iinv);
    w=(spectral_patterns*D*spectral_patterns')\spectral_patterns*D;
    for ttt=1:length(linefluorescenceseries)
        linefluorescenceseries_GFP(ttt)=sum(w(1,:).*speclinefluorescenceseries(ttt,:));
        linefluorescenceseries_YFP(ttt)=sum(w(2,:).*speclinefluorescenceseries(ttt,:));
        linefluorescenceseries_A(ttt)=sum(w(3,:).*speclinefluorescenceseries(ttt,:));
        linefluorescenceseries_Ch(ttt)=sum(w(4,:).*speclinefluorescenceseries(ttt,:));
    end
    lastbinlength=mod(length(linefluorescenceseries_GFP),binningwindow)+binningwindow;
    linefluorescenceseriesbinned_GFP=nanmean(reshape(linefluorescenceseries_GFP(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_GFP(1:end-lastbinlength))/binningwindow),1);
    linefluorescenceseriesbinned_GFP=[linefluorescenceseriesbinned_GFP nanmean(linefluorescenceseries_GFP(end-lastbinlength+1:end))];
    timeline_binned_GFP=1:1:length(linefluorescenceseriesbinned_GFP);
    timeline_binned_GFP=timeline_binned_GFP*binningwindow*scantime;
    timeline_binned_YFP=timeline_binned_GFP;
    timeline_binned_Ch=timeline_binned_GFP;
    linefluorescenceseriesbinned_YFP=nanmean(reshape(linefluorescenceseries_YFP(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_GFP(1:end-lastbinlength))/binningwindow),1);
    linefluorescenceseriesbinned_YFP=[linefluorescenceseriesbinned_YFP nanmean(linefluorescenceseries_YFP(end-lastbinlength+1:end))];
    linefluorescenceseriesbinned_A=nanmean(reshape(linefluorescenceseries_A(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_A(1:end-lastbinlength))/binningwindow),1);
    linefluorescenceseriesbinned_A=[linefluorescenceseriesbinned_A nanmean(linefluorescenceseries_A(end-lastbinlength+1:end))];
    linefluorescenceseriesbinned_Ch=nanmean(reshape(linefluorescenceseries_Ch(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_Ch(1:end-lastbinlength))/binningwindow),1);
    linefluorescenceseriesbinned_Ch=[linefluorescenceseriesbinned_Ch nanmean(linefluorescenceseries_Ch(end-lastbinlength+1:end))];

    % Plot of decomposed time series:
    hI=figure('OuterPosition',[scrsz(3) scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2],'Name','Spectrally decomposed line fluorescence');
    plot(timeline_binned_GFP,linefluorescenceseriesbinned_GFP,'Color',[26/255 150/255 65/255]);
    hold on
    plot(timeline_binned_YFP,linefluorescenceseriesbinned_YFP,'Color',[253/255 174/255 97/255]);
    hold on
    plot(timeline_binned_A,linefluorescenceseriesbinned_A,'Color',[253/255 141/255 60/255]);
    hold on
    plot(timeline_binned_Ch,linefluorescenceseriesbinned_Ch,'Color',[151/255 25/255 28/255]);

    linefluorescenceseries=linefluorescenceseries';
else
    % Binning in time to improve SNR and subsequent decomposition:    
    for binningstep=1:binningfactor
        scantime=scantime_raw*2^binningstep;
        binningwindowSNR=2;
        lastbinlengthSNR=mod(length(linefluorescenceseries),binningwindowSNR)+binningwindowSNR;
        linefluorescenceseriesbinnedSNR=mean(reshape(linefluorescenceseries(1:end-lastbinlengthSNR),binningwindowSNR,length(linefluorescenceseries(1:end-lastbinlengthSNR))/binningwindowSNR),1);
        linefluorescenceseriesbinnedSNR=[linefluorescenceseriesbinnedSNR mean(linefluorescenceseries(end-lastbinlengthSNR+1:end))];
        speclinefluorescenceseries_binned=zeros(length(linefluorescenceseriesbinnedSNR),specchannelnumb);
        for ch=1:specchannelnumb
            linefluorescenceseriesch=squeeze(speclinefluorescenceseries(:,ch));
            lastbinlengthSNRch=mod(length(linefluorescenceseriesch),binningwindowSNR)+binningwindowSNR;
            linefluorescenceseriesbinnedSNRch=mean(reshape(linefluorescenceseriesch(1:end-lastbinlengthSNRch),binningwindowSNR,length(linefluorescenceseriesch(1:end-lastbinlengthSNRch))/binningwindowSNR),1);
            linefluorescenceseriesbinnedSNRch=[linefluorescenceseriesbinnedSNRch mean(linefluorescenceseriesch(end-lastbinlengthSNRch+1:end))];
            speclinefluorescenceseries_binned(:,ch)=linefluorescenceseriesbinnedSNRch;
        end
        linefluorescenceseries_GFP=zeros(size(linefluorescenceseriesbinnedSNR));
        linefluorescenceseries_YFP=zeros(size(linefluorescenceseriesbinnedSNR));
        linefluorescenceseries_A=zeros(size(linefluorescenceseriesbinnedSNR));
        linefluorescenceseries_Ch=zeros(size(linefluorescenceseriesbinnedSNR));
        Iinv=zeros(1,specchannelnumb);
        Ittt=mean(speclinefluorescenceseries_binned,1);
        Iinv(Ittt>0)=1./Ittt(Ittt>0);
        D=diag(Iinv);
        w=(spectral_patterns*D*spectral_patterns')\spectral_patterns*D;
        for ttt=1:length(linefluorescenceseriesbinnedSNR)
            linefluorescenceseries_GFP(ttt)=sum(w(1,:).*speclinefluorescenceseries_binned(ttt,:));
            linefluorescenceseries_YFP(ttt)=sum(w(2,:).*speclinefluorescenceseries_binned(ttt,:));
            linefluorescenceseries_A(ttt)=sum(w(3,:).*speclinefluorescenceseries_binned(ttt,:));
            linefluorescenceseries_Ch(ttt)=sum(w(4,:).*speclinefluorescenceseries_binned(ttt,:));
        end
        lastbinlength=mod(length(linefluorescenceseries_GFP),binningwindow)+binningwindow;
        linefluorescenceseriesbinned_GFP=nanmean(reshape(linefluorescenceseries_GFP(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_GFP(1:end-lastbinlength))/binningwindow),1);
        linefluorescenceseriesbinned_GFP=[linefluorescenceseriesbinned_GFP nanmean(linefluorescenceseries_GFP(end-lastbinlength+1:end))];
        timeline_binned_GFP=1:1:length(linefluorescenceseriesbinned_GFP);
        timeline_binned_GFP=timeline_binned_GFP*binningwindow*scantime;
        timeline_binned_YFP=timeline_binned_GFP;
        timeline_binned_A=timeline_binned_GFP;
        timeline_binned_Ch=timeline_binned_GFP;
        linefluorescenceseriesbinned_YFP=nanmean(reshape(linefluorescenceseries_YFP(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_GFP(1:end-lastbinlength))/binningwindow),1);
        linefluorescenceseriesbinned_YFP=[linefluorescenceseriesbinned_YFP nanmean(linefluorescenceseries_YFP(end-lastbinlength+1:end))];
        linefluorescenceseriesbinned_A=nanmean(reshape(linefluorescenceseries_A(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_A(1:end-lastbinlength))/binningwindow),1);
        linefluorescenceseriesbinned_A=[linefluorescenceseriesbinned_A nanmean(linefluorescenceseries_A(end-lastbinlength+1:end))];
        linefluorescenceseriesbinned_Ch=nanmean(reshape(linefluorescenceseries_Ch(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_Ch(1:end-lastbinlength))/binningwindow),1);
        linefluorescenceseriesbinned_Ch=[linefluorescenceseriesbinned_Ch nanmean(linefluorescenceseries_Ch(end-lastbinlength+1:end))];     

        % Plot of decomposed time series:
        hI=figure('OuterPosition',[scrsz(3) scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2],'Name','Spectrally decomposed line fluorescence');
        plot(timeline_binned_GFP,linefluorescenceseriesbinned_GFP,'Color',[26/255 150/255 65/255]);
        hold on
        plot(timeline_binned_YFP,linefluorescenceseriesbinned_YFP,'Color',[253/255 174/255 97/255]);
        hold on
        plot(timeline_binned_A,linefluorescenceseriesbinned_A,'Color',[253/255 141/255 60/255]);
        hold on
        plot(timeline_binned_Ch,linefluorescenceseriesbinned_Ch,'Color',[151/255 25/255 28/255]);

        linefluorescenceseries=linefluorescenceseriesbinnedSNR';
        speclinefluorescenceseries=speclinefluorescenceseries_binned;
    end
end
% Plot of average spectra and calculation of spectral fractions:    
figure('Name','Fluorescent Spectra')
spectraldecomposfitfunc=@(x,t) x(1)*spectral_patterns(1,:)+x(2)*spectral_patterns(2,:)+x(3)*spectral_patterns(3,:)+x(4)*spectral_patterns(4,:);%+x(4)*avgspectrum_Card;
x0spec=[0.25 0.25 0.25 0.25];
fixedspec=[false false false false];
[fitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixedspec,spectralchannels,emissionspec',spectraldecomposfitfunc,x0spec);
fractionGFP=fitparameters(1)
fractionYFP=fitparameters(2)
fractionA=fitparameters(3)
fractionCh=fitparameters(4)
fractions=[fractionGFP fractionYFP fractionA fractionCh];
avgspectrum_all4_fit=spectraldecomposfitfunc(fitparameters,spectralchannels);  
fidfrac=fopen([path2  '\specfractions_GFP_YFP_A_Ch2.txt'],'a'); % adjust file name if needed, file to save spectral fractions
fprintf(fidfrac,inputfilename(1:end-4));
fprintf(fidfrac,'\t %e\t %e\t %e\t %e\n',fractions');
    
fidw=fopen([path2 '\' inputfilename(1:end-4) '_photonweights.txt'],'a'); % adjust file name if needed, photon weights are saved as .txt
weightsoutput=zeros(length(spectralchannels),5);
weightsoutput(:,1)=spectralchannels';
weightsoutput(:,2:5)=w';
fprintf(fidw,'%e\t %e\t %e\t %e\t %e\n',weightsoutput');
    
plot(spectralchannels,emissionspec);
hold on
plot(spectralchannels,avgspectrum_all4_fit,'--');
xlabel('Channel [nm]')
ylabel('Norm.emission')

% Dialogue: If only one species is present in the measurement, the analysis
% can be interrupted here
str_fractions='Spectral fractions are: %.2f (G), %.2f (Y), %.2f (A), %.2f (Ch2). Continue with FCS analysis?';
questionstring=sprintf(str_fractions,fractionGFP,fractionYFP,fractionA,fractionCh);
answer=questdlg(questionstring, ...
'Continue?', ...
'Yes','No','Yes');
switch answer
    case 'Yes'
        continueFCS=1;
    case 'No'
        continueFCS=0;
end
clear speclinefluorescenceseries speclinefluorescenceseries_binned
if continueFCS==0
    return
end
    
% Calculation of bleaching fractions and bleaching correction: 
line1fluorescenceseries=linefluorescenceseries_GFP;
line2fluorescenceseries=linefluorescenceseries_YFP;
line3fluorescenceseries=linefluorescenceseries_A;
line4fluorescenceseries=linefluorescenceseries_Ch;
bleachingfraction1=1-nanmean(linefluorescenceseries_GFP(end-250+1:end))/nanmean(linefluorescenceseries_GFP(1:250));
bleachingfraction2=1-nanmean(linefluorescenceseries_YFP(end-250+1:end))/nanmean(linefluorescenceseries_YFP(1:250));
bleachingfraction3=1-nanmean(linefluorescenceseries_A(end-250+1:end))/nanmean(linefluorescenceseries_A(1:250));
bleachingfraction4=1-nanmean(linefluorescenceseries_Ch(end-250+1:end))/nanmean(linefluorescenceseries_Ch(1:250));
timeline1=1:1:length(linefluorescenceseries_GFP);
timeline1=timeline1*scantime;
timeline2=timeline1;
timeline3=timeline1;
timeline4=timeline1;
line1fluorescenceseriesbinned=linefluorescenceseriesbinned_GFP;
line2fluorescenceseriesbinned=linefluorescenceseriesbinned_YFP;
line3fluorescenceseriesbinned=linefluorescenceseriesbinned_A;
line4fluorescenceseriesbinned=linefluorescenceseriesbinned_Ch;
timeline1_binned=timeline_binned_GFP;
timeline2_binned=timeline_binned_YFP;
timeline3_binned=timeline_binned_A;
timeline4_binned=timeline_binned_Ch;

% Depletion correction of time trace:
if intensityselection==1
    Ifull1=line1fluorescenceseriesbinned;
    Ifull2=line2fluorescenceseriesbinned;
    Ifull3=line3fluorescenceseriesbinned;
    Ifull4=line4fluorescenceseriesbinned;
    timeline1binned=timeline_binned_GFP;
    timeline2binned=timeline_binned_YFP;
    timeline3binned=timeline_binned_A;
    timeline4binned=timeline_binned_Ch;
    segmentlength=length(Ifull1)/numberofsegments;
    segmentsincl1=ones(1,numberofsegments);
    segmentsincl2=ones(1,numberofsegments);
    segmentsincl3=ones(1,numberofsegments);
    segmentsincl4=ones(1,numberofsegments);
        
    % GUI for selection of intensity segments:
    Intselection4colors
    goon=0;
    while goon==0
        pause(5)
    end
    line1fluorescenceseries_corrected=depletioncorrectionfunc(timeline1,line1fluorescenceseries,ff1);
    line2fluorescenceseries_corrected=depletioncorrectionfunc(timeline2,line2fluorescenceseries,ff2);
    line3fluorescenceseries_corrected=depletioncorrectionfunc(timeline3,line3fluorescenceseries,ff3);
    line4fluorescenceseries_corrected=depletioncorrectionfunc(timeline4,line4fluorescenceseries,ff4);
    lastbinlength=mod(length(line1fluorescenceseries),binningwindow)+binningwindow;
    line1fluorescenceseries_corrected_binned=mean(reshape(line1fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line1fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line1fluorescenceseries_corrected_binned=[line1fluorescenceseries_corrected_binned mean(line1fluorescenceseries_corrected(end-lastbinlength+1:end))];
    line2fluorescenceseries_corrected_binned=mean(reshape(line2fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line2fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line2fluorescenceseries_corrected_binned=[line2fluorescenceseries_corrected_binned mean(line2fluorescenceseries_corrected(end-lastbinlength+1:end))];
    line3fluorescenceseries_corrected_binned=mean(reshape(line3fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line3fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line3fluorescenceseries_corrected_binned=[line3fluorescenceseries_corrected_binned mean(line3fluorescenceseries_corrected(end-lastbinlength+1:end))];
    line4fluorescenceseries_corrected_binned=mean(reshape(line4fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line4fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line4fluorescenceseries_corrected_binned=[line4fluorescenceseries_corrected_binned mean(line4fluorescenceseries_corrected(end-lastbinlength+1:end))];
    figure('OuterPosition',[scrsz(3)/3 scrsz(4)/4 2*scrsz(3)/3 scrsz(4)/4],'Name','Corrected Membrane Fluorescence time series')
    subplot(4,1,1),plot(timeline1_binned,line1fluorescenceseries_corrected_binned); % Diesen Plot auch gebinnt darstellen!
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch1')
    subplot(4,1,2),plot(timeline2_binned,line2fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch2')
    subplot(4,1,3),plot(timeline3_binned,line3fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch3')
    subplot(4,1,4),plot(timeline4_binned,line4fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch4')
else
    if depletioncorrection==1
        line1fluorescenceseries_corrected=depletioncorrectionfunc(timeline1,line1fluorescenceseries,f1); 
        line2fluorescenceseries_corrected=depletioncorrectionfunc(timeline2,line2fluorescenceseries,f2);
        line3fluorescenceseries_corrected=depletioncorrectionfunc(timeline3,line3fluorescenceseries,f3);
        line4fluorescenceseries_corrected=depletioncorrectionfunc(timeline4,line4fluorescenceseries,f4);
        lastbinlength=mod(length(line1fluorescenceseries),binningwindow)+binningwindow;
        line1fluorescenceseries_corrected_binned=mean(reshape(line1fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line1fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
        line1fluorescenceseries_corrected_binned=[line1fluorescenceseries_corrected_binned mean(line1fluorescenceseries_corrected(end-lastbinlength+1:end))];
        line2fluorescenceseries_corrected_binned=mean(reshape(line2fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line2fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
        line2fluorescenceseries_corrected_binned=[line2fluorescenceseries_corrected_binned mean(line2fluorescenceseries_corrected(end-lastbinlength+1:end))];
        line3fluorescenceseries_corrected_binned=mean(reshape(line3fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line3fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
        line3fluorescenceseries_corrected_binned=[line3fluorescenceseries_corrected_binned mean(line3fluorescenceseries_corrected(end-lastbinlength+1:end))];
        line4fluorescenceseries_corrected_binned=mean(reshape(line4fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line4fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
        line4fluorescenceseries_corrected_binned=[line4fluorescenceseries_corrected_binned mean(line4fluorescenceseries_corrected(end-lastbinlength+1:end))];
        figure('Name','Corrected Membrane Fluorescence time series')
        subplot(1,4,1),plot(timeline1_binned,line1fluorescenceseries_corrected_binned);
        xlabel('time (linenumber)')
        ylabel('line fluorescence Ch1')
        hold on
        subplot(1,4,2),plot(timeline2_binned,line2fluorescenceseries_corrected_binned);
        xlabel('time (linenumber)')
        ylabel('line fluorescence Ch2')
        subplot(1,4,3),plot(timeline3_binned,line3fluorescenceseries_corrected_binned);
        xlabel('time (linenumber)')
        ylabel('line fluorescence Ch3')
        subplot(1,4,4),plot(timeline4_binned,line4fluorescenceseries_corrected_binned);
        xlabel('time (linenumber)')
        ylabel('line fluorescence Ch4')
    end
 end
 
 % Calculation of correlation functions:
 fprintf('Calculating correlation using multiple tau...\n');
 if depletioncorrection==1
     [tcorr1,fautocorr1,sigmas1]=autocorrFCSmultipletau(line1fluorescenceseries_corrected);
     [tcorr2,fautocorr2,sigmas2]=autocorrFCSmultipletau(line2fluorescenceseries_corrected);
     [tcorr3,fautocorr3,sigmas3]=autocorrFCSmultipletau(line3fluorescenceseries_corrected);
     [tcorr4,fautocorr4,sigmas4]=autocorrFCSmultipletau(line4fluorescenceseries_corrected);
     [tcorrcc12,fcrosscorr12,sigmascc12]=crosscorrFCSmultipletau(line1fluorescenceseries_corrected,line2fluorescenceseries_corrected);
     [tcorrcc13,fcrosscorr13,sigmascc13]=crosscorrFCSmultipletau(line1fluorescenceseries_corrected,line3fluorescenceseries_corrected);
     [tcorrcc14,fcrosscorr14,sigmascc14]=crosscorrFCSmultipletau(line1fluorescenceseries_corrected,line4fluorescenceseries_corrected);
     [tcorrcc23,fcrosscorr23,sigmascc23]=crosscorrFCSmultipletau(line2fluorescenceseries_corrected,line3fluorescenceseries_corrected);
     [tcorrcc24,fcrosscorr24,sigmascc24]=crosscorrFCSmultipletau(line2fluorescenceseries_corrected,line4fluorescenceseries_corrected);
     [tcorrcc34,fcrosscorr34,sigmascc34]=crosscorrFCSmultipletau(line3fluorescenceseries_corrected,line4fluorescenceseries_corrected);
 else
     [tcorr1,fautocorr1,sigmas1]=autocorrFCSmultipletau(line1fluorescenceseries);
     [tcorr2,fautocorr2,sigmas2]=autocorrFCSmultipletau(line2fluorescenceseries);
     [tcorr3,fautocorr3,sigmas3]=autocorrFCSmultipletau(line3fluorescenceseries);
     [tcorr4,fautocorr4,sigmas4]=autocorrFCSmultipletau(line4fluorescenceseries);
     [tcorrcc12,fcrosscorr12,sigmascc12]=crosscorrFCSmultipletau(line1fluorescenceseries,line2fluorescenceseries);
     [tcorrcc13,fcrosscorr13,sigmascc13]=crosscorrFCSmultipletau(line1fluorescenceseries,line3fluorescenceseries);
     [tcorrcc14,fcrosscorr14,sigmascc14]=crosscorrFCSmultipletau(line1fluorescenceseries,line4fluorescenceseries);
     [tcorrcc23,fcrosscorr23,sigmascc23]=crosscorrFCSmultipletau(line2fluorescenceseries,line3fluorescenceseries);
     [tcorrcc24,fcrosscorr24,sigmascc24]=crosscorrFCSmultipletau(line2fluorescenceseries,line4fluorescenceseries);
     [tcorrcc34,fcrosscorr34,sigmascc34]=crosscorrFCSmultipletau(line3fluorescenceseries,line4fluorescenceseries);
 end
 tcorr=tcorr1*scantime; %all tcorr's are equal anyways
 tcorrcc=tcorrcc12*scantime; %In this case, all tcorr's are equal anyways
    
weights1=fautocorr1./sigmas1;
weights2=fautocorr2./sigmas2;
weights3=fautocorr3./sigmas3;
weights4=fautocorr4./sigmas4;
weightscc12=fcrosscorr12./sigmascc12;
weightscc13=fcrosscorr13./sigmascc13;
weightscc14=fcrosscorr14./sigmascc14;
weightscc23=fcrosscorr23./sigmascc23;
weightscc24=fcrosscorr24./sigmascc24;
weightscc34=fcrosscorr34./sigmascc34;
    

% Plot of CFs:
figure('OuterPosition',[scrsz(3) scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Correlation functions')
semilogx(tcorr,fautocorr1,'.','Color',[26/255 150/255 65/255])
pause
hold on
semilogx(tcorr,fautocorr2,'.','Color',[253/255 174/255 97/255])
pause
hold on
semilogx(tcorr,fautocorr3,'.','Color',[253/255 141/255 60/255])
pause
hold on
semilogx(tcorr,fautocorr4,'.','Color',[151/255 25/255 28/255])
pause
hold on
semilogx(tcorrcc,fcrosscorr12,'-','Color', [37/255 52/255 148/255])
pause
hold on
semilogx(tcorrcc,fcrosscorr13,'-','Color', [65/255 182/255 196/255])
pause
hold on
semilogx(tcorrcc,fcrosscorr14,'-','Color', [252/255 141/255 89/255])
pause
hold on
semilogx(tcorrcc,fcrosscorr23,'-','Color', [254/255 153/255 41/255])
pause
hold on
semilogx(tcorrcc,fcrosscorr24,'-','Color', [227/255 74/255 51/255])
pause
hold on
semilogx(tcorrcc,fcrosscorr34,'-','Color', [179/255 0/255 0/255])

if depletioncorrection==1
    line1fluorescenceseries=line1fluorescenceseries_corrected;
    line2fluorescenceseries=line2fluorescenceseries_corrected;
    line3fluorescenceseries=line3fluorescenceseries_corrected;
    line4fluorescenceseries=line4fluorescenceseries_corrected;
end

lastbinlength=mod(length(line1fluorescenceseries),binningwindow_GUI)+binningwindow_GUI;
line1fluorescenceseriesbinned_GUI=mean(reshape(line1fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line1fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line1fluorescenceseriesbinned_GUI=[line1fluorescenceseriesbinned_GUI mean(line1fluorescenceseries(end-lastbinlength+1:end))];
line2fluorescenceseriesbinned_GUI=mean(reshape(line2fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line2fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line2fluorescenceseriesbinned_GUI=[line2fluorescenceseriesbinned_GUI mean(line2fluorescenceseries(end-lastbinlength+1:end))];
line3fluorescenceseriesbinned_GUI=mean(reshape(line3fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line3fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line3fluorescenceseriesbinned_GUI=[line3fluorescenceseriesbinned_GUI mean(line3fluorescenceseries(end-lastbinlength+1:end))];
line4fluorescenceseriesbinned_GUI=mean(reshape(line4fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line4fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line4fluorescenceseriesbinned_GUI=[line4fluorescenceseriesbinned_GUI mean(line4fluorescenceseries(end-lastbinlength+1:end))];

linefluorescenceseries_Chs=[line1fluorescenceseries' line2fluorescenceseries' line3fluorescenceseries' line4fluorescenceseries'];
linefluorescenceseriesbinned_Chs_GUI=[line1fluorescenceseriesbinned_GUI' line2fluorescenceseriesbinned_GUI' line3fluorescenceseriesbinned_GUI' line4fluorescenceseriesbinned_GUI'];

% Segment-wise ACF analysis: A plot of CFs is generated which shows the CFs
% of all segments (overlayed)
if ACFselection==1
    segmentlength=length(line1fluorescenceseries)/ACFsegmentnumb;
    Isegmentlength_GUI=segmentlength/binningwindow_GUI;
    Itrace_Chs=zeros(segmentlength,ACFsegmentnumb+1,4);
    Itrace_Chs_GUI=zeros(Isegmentlength_GUI,ACFsegmentnumb+1,4);
    segsACF=figure('OuterPosition',[4*scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','ACFs segments');
    segsCCF=figure('OuterPosition',[5*scrsz(3)/3 scrsz(4)/4 scrsz(3)/3 3*scrsz(4)/4],'Name','CCFs segments');
    [tcorr11,fcorr11,sigmas11]=autocorrFCSmultipletau(linefluorescenceseries_Chs(1:segmentlength,1));
    correlationcurves_Chs=zeros(size(fcorr11,2),ACFsegmentnumb+1,4);
    correlationcurvesCC_Chs=zeros(size(fcorr11,2),ACFsegmentnumb+1,6);
    sigmascurves_Chs=correlationcurves_Chs;
    sigmascurvesCC_Chs=correlationcurvesCC_Chs;
    colormatrix_Chs=[26/255 150/255 65/255;253/255 174/255 97/255;153/255 142/255 195/255;151/255 25/255 28/255];
    colormatrixCC_Chs=[37/255 52/255 148/255;65/255 182/255 196/255;252/255 141/255 89/255;254/255 153/255 41/255;227/255 74/255 51/255;179/255 0/255 0/255];
    for i=1:size(Itrace_Chs,2)-1
        ccindex=1;
        for j=1:4
            Itrace_Chs(:,i+1,j)=linefluorescenceseries_Chs(segmentlength*(i-1)+1:i*segmentlength,j);
            Itrace_Chs_GUI(:,i+1,j)=linefluorescenceseriesbinned_Chs_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI,j);
            [tcorrij,fcorrij,sigmasij]=autocorrFCSmultipletau(Itrace_Chs(:,i+1,j));
            correlationcurves_Chs(:,i+1,j)=fcorrij';
            sigmascurves_Chs(:,i,j)=sigmasij';
            figure(segsACF)
            hold on
            subplot(2,2,j),semilogx(tcorrij,fcorrij,'Color',colormatrix_Chs(j,:),'LineStyle','none','Marker','.','MarkerFaceColor',colormatrix_Chs(j,:));
            xlabel('\tau [s]')
            ylabel('Correlation')
            if j<4
                for k=j+1:4
                    Itrace_Chs(:,i+1,k)=linefluorescenceseries_Chs(segmentlength*(i-1)+1:i*segmentlength,k);
                    Itrace_Chs_GUI(:,i+1,k)=linefluorescenceseriesbinned_Chs_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI,k);
                    [tcorrccijk,fcrosscorrijk,sigmasccijk]=crosscorrFCSmultipletau(Itrace_Chs(:,i+1,j),Itrace_Chs(:,i+1,k));
                    correlationcurvesCC_Chs(:,i+1,ccindex)=fcrosscorrijk';
                    sigmascurvesCC_Chs(:,i,ccindex)=sigmasccijk';
                    figure(segsCCF)
                    hold on
                    subplot(3,2,ccindex),semilogx(tcorrccijk,fcrosscorrijk,'Color',colormatrixCC_Chs(ccindex,:),'LineStyle','none','Marker','.','MarkerFaceColor',colormatrixCC_Chs(ccindex,:));
                    xlabel('\tau [s]')
                    ylabel('Correlation')
                    ccindex=ccindex+1;
                end
            end
            tcorrij=tcorrij*scantime; %All tcorri's equal anyways
            tcorrccijk=tcorrccijk*scantime; %All tcorrcci's equal anyways
        end
        lbi=min(tcorrij);
        ubi=max(tcorrij);
        %pause
        hold on
    end
    ItraceII_Chs=Itrace_Chs_GUI;
    segmenttime_GUI=scantime:scantime*binningwindow_GUI:scantime*segmentlength;
    for i=1:size(ItraceII_Chs,3)
        ItraceII_Chs(:,1,i)=segmenttime_GUI';
    end
    for i=1:size(correlationcurves_Chs,3)
        correlationcurves_Chs(:,1,i)=tcorrij';
    end
    for i=1:size(correlationcurvesCC_Chs,3)
        correlationcurvesCC_Chs(:,1,i)=tcorrccijk';
    end
    curveincl_Chs=ones(ACFsegmentnumb,4);

    % GUI for selection of segment ACFs:
    Corrselection4Channels
    goon=0;
    while goon==0
        pause(5)
    end
    fprintf('Done!')
    for i=1:size(curveincl_Chs,1)
        for j=1:4
            if curveincl_Chs(i,j)==0
                linefluorescenceseries_Chs(segmentlength*(i-1)+1:i*segmentlength,j)=NaN;
            end
        end     
    end
    
    % Calculation of final CFs:
    [tcorr1final,fcorr1final,sigmas1final]=autocorrFCSmultipletau(linefluorescenceseries_Chs(:,1));
    fcorrfinal_Chs=zeros(length(fcorr1final),4);
    sigmasfinal_Chs=fcorrfinal_Chs;
    fccfinal_Chs=zeros(length(fcorr1final),6);
    sigmasccfinal_Chs=fccfinal_Chs;
    ccindex=1;
    for i=1:4
        [tcorrifinal,fcorrifinal,sigmasifinal]=autocorrFCSmultipletau(linefluorescenceseries_Chs(:,i));
        fcorrfinal_Chs(:,i)=fcorrifinal;
        sigmasfinal_Chs(:,i)=sigmasifinal;
        if i<4
            for j=i+1:4
                [tcorrccijfinal,fccijfinal,sigmasccijfinal]=crosscorrFCSmultipletau(linefluorescenceseries_Chs(:,i),linefluorescenceseries_Chs(:,j));
                fccfinal_Chs(:,ccindex)=fccijfinal;
                sigmasccfinal_Chs(:,ccindex)=sigmasccijfinal;
                ccindex=ccindex+1;
            end   
        end
    end
    tcorrfinal=tcorrifinal*scantime;
    tcorrccfinal=tcorrccijfinal*scantime;
    weightsfinal_Chs=abs(fcorrfinal_Chs./sigmasfinal_Chs);
    weightsccfinal_Chs=abs(fccfinal_Chs./sigmasccfinal_Chs);
    
    % Average of segments:
    fcorrsegfinal_Chs=zeros(length(tcorrij),4);
    sigmassegfinal_Chs=fcorrsegfinal_Chs;
    fccsegfinal_Chs=zeros(length(tcorrccijk),6);
    sigmasccsegfinal_Chs=fccsegfinal_Chs;
    for i=1:size(curveincl_Chs,1);
        ccindex=1;
        for j=1:4
            if curveincl_Chs(i,j)==1
                fcorrsegfinal_Chs(:,j)=fcorrsegfinal_Chs(:,j)+correlationcurves_Chs(:,1+i,j);
                sigmassegfinal_Chs(:,j)=sigmassegfinal_Chs(:,j)+sigmascurves_Chs(:,i,j);
            end
            if j<4
                for k=j+1:4
                    if (curveincl_Chs(i,j)+curveincl_Chs(i,k))==2
                        fccsegfinal_Chs(:,ccindex)=fccsegfinal_Chs(:,ccindex)+correlationcurvesCC_Chs(:,i+1,ccindex);
                        sigmasccsegfinal_Chs(:,ccindex)=sigmasccsegfinal_Chs(:,ccindex)+sigmascurvesCC_Chs(:,i,ccindex);
                    end
                    ccindex=ccindex+1;
                end  
            end
        end
    end
    ccindex=1;
    for i=1:4
        fcorrsegfinal_Chs(:,i)=fcorrsegfinal_Chs(:,i)/sum(curveincl_Chs(:,i));
        sigmassegfinal_Chs(:,i)=sigmassegfinal_Chs(:,i)/sum(curveincl_Chs(:,i));
        if i<4
            for j=i+1:4
                fccsegfinal_Chs(:,ccindex)=fccsegfinal_Chs(:,ccindex)/sum((curveincl_Chs(:,i)+curveincl_Chs(:,j))==2);
                sigmasccsegfinal_Chs(:,ccindex)=sigmasccsegfinal_Chs(:,ccindex)/sum((curveincl_Chs(:,i)+curveincl_Chs(:,j))==2);
                ccindex=ccindex+1;
            end
        end
    end
  
    weightssegfinal_Chs=abs(fcorrsegfinal_Chs./sigmassegfinal_Chs);
    tcorrsegfinal=tcorrij;
    weightsccsegfinal_Chs=abs(fccsegfinal_Chs./sigmasccsegfinal_Chs);
    tcorrccsegfinal=tcorrij; %for now
    
    % Fit of final CFs with diffusion model:
    fprintf('Final Fitting...\n');
    lb=scantime;
    lbcc=lb;
    lbseg=lb;
    lbccseg=lb;
    ub=max(tcorrfinal);
    ubcc=ub;
    ubseg=max(tcorrsegfinal);
    ubccseg=ubseg;
    
    weightsfinal_Chs(weightsfinal_Chs==0)=10^-4;
    weightsccfinal_Chs(weightsccfinal_Chs==0)=10^-4;
    
    lsautofitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
    flscrossfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);
    
    Nfinal_Chs=zeros(4,1);
    taudfinal_Chs=Nfinal_Chs;
    Sfitfinal_Chs=Nfinal_Chs;
    CIfinal_Chs=zeros(12,2);
    Nccfinal_Chs=zeros(6,1);
    taudccfinal_Chs=Nccfinal_Chs;
    Sfitccfinal_Chs=Nccfinal_Chs;
    CIccfinal_Chs=zeros(18,2);
    
    % Fit of whole time series (try) and plot:
    try
        h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curves');
        positionvector1=[0.1 0.35 0.8 0.55];
        positionvector2=[0.1 0.1 0.8 0.15];
        for i=1:4
            [Nfinal,taudfinal,Sfitfinal,CIfinal,fitcurvefinal,residualsfinal] = autocorrfit2Ddiff(tcorrfinal,fcorrfinal_Chs(:,i)',lb,ub,weightsfinal_Chs(:,i)',lsautofitfunc,x0,fixed);
            Nfinal_Chs(i)=Nfinal;
            taudfinal_Chs(i)=taudfinal;
            Sfitfinal_Chs(i)=Sfitfinal;
            CIfinal_Chs((i-1)*3+1:i*3,:)=CIfinal;
            subplot('Position',positionvector1),semilogx(tcorrfinal,fitcurvefinal,'Color',colormatrix_Chs(i,:),'LineWidth',2)
            hold on
            subplot('Position',positionvector1),semilogx(tcorrfinal,fcorrfinal_Chs(:,i),'d','Color',colormatrix_Chs(i,:),'LineWidth',1)
            hold on
            xlabel('Lag time [s]')
            ylabel('Correlation')
            subplot('Position',positionvector2),semilogx(tcorrfinal,residualsfinal,'d','Color',colormatrix_Chs(i,:),'LineWidth',1)
            hold on
        end
        for i=1:6 
            [Nccfinal,taudccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiff(tcorrccfinal,fccfinal_Chs(:,i)',lbcc,ubcc,weightsccfinal_Chs(:,i)',lsautofitfunc,y0,fixed);
            Nccfinal_Chs(i)=Nccfinal;
            taudccfinal_Chs(i)=taudccfinal;
            Sfitccfinal_Chs(i)=Sfitccfinal;
            CIccfinal_Chs((i-1)*3+1:i*3,:)=CIccfinal;
            subplot('Position',positionvector1),semilogx(tcorrccfinal,fitcurveccfinal,'Color',colormatrixCC_Chs(i,:),'LineWidth',2)
            hold on
            subplot('Position',positionvector1),semilogx(tcorrccfinal,fccfinal_Chs(:,i),'s','Color',colormatrixCC_Chs(i,:),'LineWidth',1)
            hold on
            subplot('Position',positionvector2),semilogx(tcorrccfinal,residualsccfinal,'s','Color',colormatrixCC_Chs(i,:),'LineWidth',1)
            hold on
            xlabel('Lag time [s]')
            ylabel('Residuals')
        end
        subplot('Position',positionvector2),semilogx(tcorrfinal,zeros(size(tcorrfinal)),'-k','LineWidth',1.5);
        
        % Relative cross-correlation amplitude:
        relCC12=max([Nfinal_Chs(1)/Nccfinal_Chs(1) Nfinal_Chs(2)/Nccfinal_Chs(1)]);
        relCC13=max([Nfinal_Chs(1)/Nccfinal_Chs(2) Nfinal_Chs(3)/Nccfinal_Chs(2)]);
        relCC14=max([Nfinal_Chs(1)/Nccfinal_Chs(3) Nfinal_Chs(4)/Nccfinal_Chs(3)]);
        relCC23=max([Nfinal_Chs(2)/Nccfinal_Chs(4) Nfinal_Chs(3)/Nccfinal_Chs(4)]);
        relCC24=max([Nfinal_Chs(2)/Nccfinal_Chs(5) Nfinal_Chs(4)/Nccfinal_Chs(5)]);
        relCC34=max([Nfinal_Chs(3)/Nccfinal_Chs(6) Nfinal_Chs(4)/Nccfinal_Chs(6)]);

        errorfull=0;
    catch
        errorfull=1;
    end
    
    Nsegfinal_Chs=zeros(4,1);
    taudsegfinal_Chs=Nsegfinal_Chs;
    Sfitsegfinal_Chs=Nsegfinal_Chs;
    CIsegfinal_Chs=zeros(12,2);
    Nccsegfinal_Chs=zeros(6,1);
    taudccsegfinal_Chs=Nccsegfinal_Chs;
    Sfitccsegfinal_Chs=Nccsegfinal_Chs;
    CIccsegfinal_Chs=zeros(18,2);
    
    % Fit of segment-averaged time series and plot:
    hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curves');
    positionvector1=[0.1 0.35 0.8 0.55];
    positionvector2=[0.1 0.1 0.8 0.15];
    for i=1:4
        [Nsegfinal,taudsegfinal,Sfitsegfinal,CIsegfinal,fitcurvesegfinal,residualssegfinal] = autocorrfit2Ddiff(tcorrsegfinal,fcorrsegfinal_Chs(:,i)',lbseg,ubseg,weightssegfinal_Chs(:,i)',lsautofitfunc,x0,fixed);
        Nsegfinal_Chs(i)=Nsegfinal;
        taudsegfinal_Chs(i)=taudsegfinal;
        Sfitsegfinal_Chs(i)=Sfitsegfinal;
        CIsegfinal_Chs((i-1)*3+1:i*3,:)=CIsegfinal;
        subplot('Position',positionvector1),semilogx(tcorrsegfinal,fitcurvesegfinal,'Color',colormatrix_Chs(i,:),'LineWidth',2)
        hold on
        subplot('Position',positionvector1),semilogx(tcorrsegfinal,fcorrsegfinal_Chs(:,i),'d','Color',colormatrix_Chs(i,:),'LineWidth',1)
        hold on
        xlabel('Lag time [s]')
        ylabel('Correlation')
        subplot('Position',positionvector2),semilogx(tcorrsegfinal,residualssegfinal,'d','Color',colormatrix_Chs(i,:),'LineWidth',1)
        hold on
    end
    for i=1:6 
        [Nccsegfinal,taudccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiff(tcorrccsegfinal,fccsegfinal_Chs(:,i)',lbccseg,ubccseg,weightsccsegfinal_Chs(:,i)',lsautofitfunc,y0,fixed);
        Nccsegfinal_Chs(i)=Nccsegfinal;
        taudccsegfinal_Chs(i)=taudccsegfinal;
        Sfitccsegfinal_Chs(i)=Sfitccsegfinal;
        CIccsegfinal_Chs((i-1)*3+1:i*3,:)=CIccsegfinal;
        subplot('Position',positionvector1),semilogx(tcorrccsegfinal,fitcurveccsegfinal,'Color',colormatrixCC_Chs(i,:),'LineWidth',2)
        hold on
        subplot('Position',positionvector1),semilogx(tcorrccsegfinal,fccsegfinal_Chs(:,i),'s','Color',colormatrixCC_Chs(i,:),'LineWidth',1)
        hold on
        subplot('Position',positionvector2),semilogx(tcorrccsegfinal,residualsccsegfinal,'s','Color',colormatrixCC_Chs(i,:),'LineWidth',1)
        hold on
        xlabel('Lag time [s]')
        ylabel('Residuals')
    end
    subplot('Position',positionvector2),semilogx(tcorrsegfinal,zeros(size(tcorrsegfinal)),'-k','LineWidth',1.5);

    % Relative cross-correlation amplitude:
    relCC12seg=max([Nsegfinal_Chs(1)/Nccsegfinal_Chs(1) Nsegfinal_Chs(2)/Nccsegfinal_Chs(1)]);
    relCC13seg=max([Nsegfinal_Chs(1)/Nccsegfinal_Chs(2) Nsegfinal_Chs(3)/Nccsegfinal_Chs(2)]);
    relCC14seg=max([Nsegfinal_Chs(1)/Nccsegfinal_Chs(3) Nsegfinal_Chs(4)/Nccsegfinal_Chs(3)]);
    relCC23seg=max([Nsegfinal_Chs(2)/Nccsegfinal_Chs(4) Nsegfinal_Chs(3)/Nccsegfinal_Chs(4)]);
    relCC24seg=max([Nsegfinal_Chs(2)/Nccsegfinal_Chs(5) Nsegfinal_Chs(4)/Nccsegfinal_Chs(5)]);
    relCC34seg=max([Nsegfinal_Chs(3)/Nccsegfinal_Chs(6) Nsegfinal_Chs(4)/Nccsegfinal_Chs(6)]);
    
    % Save output (fit parameters, ACFs etc):
    fid3=fopen([path2 '\' inputfilename(1:end-4) '_final_ACFseg.txt'],'a'); % adjust path if necessary!
    fid4=fopen([path2 '\' inputfilename(1:end-4) '_final_fitparameters_seg.txt'],'a'); % adjust path if necessary!
    fid5=fopen([path2 '\' inputfilename(1:end-4) 'channel_Is.mat'],'a'); % adjust path if necessary!
    fid6=fopen([path2 '\' inputfilename(1:end-4) '_final_curveincls.txt'],'a'); % adjust path if necessary!
    
    linechannelfluorescenceseries=[line1fluorescenceseries' line2fluorescenceseries' line3fluorescenceseries'];
    save([path2 '\' inputfilename(1:end-4) 'channel_Is.mat'],'linechannelfluorescenceseries')
    
    
    if errorfull==0
        % Full analysis:
        fid1=fopen([path2 '\' inputfilename(1:end-4) '_final_ACF.txt'],'a'); % adjust path if necessary!
        fid2=fopen([path2 '\' inputfilename(1:end-4) '_final_fitparameters.txt'],'a'); % adjust path if necessary!
        saveas(h,[path2 '\' inputfilename(1:end-4) ' CCFs_full.fig'])
        output=zeros(length(tcorrfinal),22);    

        output(:,1)=tcorrfinal;
        for i=1:4
            output(:,2*i)=fcorrfinal_Chs(:,i);
            output(:,2*i+1)=weightsfinal_Chs(:,i);
        end
        output(:,10)=tcorrccfinal;
        for i=1:6
            output(:,9+2*i)=fccfinal_Chs(:,i);
            output(:,10+2*i)=weightsccfinal_Chs(:,i);
        end
        
        outputparameters=zeros(38,3);
        for i=1:4
            outputparameters((i-1)*3+1:i*3,1)=[Nfinal_Chs(i);taudfinal_Chs(i);Sfitfinal_Chs(i)];
            outputparameters((i-1)*3+1:i*3,2:end)=CIfinal_Chs((i-1)*3+1:i*3,:);
        end
        for i=1:6
            outputparameters((i-1)*3+13:i*3+12,1)=[Nccfinal_Chs(i);taudccfinal_Chs(i);Sfitccfinal_Chs(i)];
            outputparameters((i-1)*3+13:i*3+12,2:end)=CIccfinal_Chs((i-1)*3+1:i*3,:);
        end

        outputparameters(31:34,1)=[bleachingfraction1;bleachingfraction2;bleachingfraction3;bleachingfraction4];
        outputparameters(35:38,1)=nanmean(linefluorescenceseries_Chs,1)';
        
        fprintf(fid1,'%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output');
        fprintf(fid2,'%e\t %e\t %e\n',outputparameters');
    end
    
    % Segment analysis:
    outputseg=zeros(length(tcorrsegfinal),22);    
    outputseg(:,1)=tcorrsegfinal;
    for i=1:4
        outputseg(:,2*i)=fcorrsegfinal_Chs(:,i);
        outputseg(:,2*i+1)=weightssegfinal_Chs(:,i);
    end
    outputseg(:,10)=tcorrccsegfinal;
    for i=1:6
        outputseg(:,9+2*i)=fccsegfinal_Chs(:,i);
        outputseg(:,10+2*i)=weightsccsegfinal_Chs(:,i);
    end

    outputparametersseg=zeros(38,3);
    for i=1:4
        outputparametersseg((i-1)*3+1:i*3,1)=[Nsegfinal_Chs(i);taudsegfinal_Chs(i);Sfitsegfinal_Chs(i)];
        outputparametersseg((i-1)*3+1:i*3,2:end)=CIsegfinal_Chs((i-1)*3+1:i*3,:);
    end
    for i=1:6
        outputparametersseg((i-1)*3+13:i*3+12,1)=[Nccsegfinal_Chs(i);taudccsegfinal_Chs(i);Sfitccsegfinal_Chs(i)];
        outputparametersseg((i-1)*3+13:i*3+12,2:end)=CIccsegfinal_Chs((i-1)*3+1:i*3,:);
    end

    outputparametersseg(31:34,1)=[bleachingfraction1;bleachingfraction2;bleachingfraction3;bleachingfraction4];
    outputparametersseg(35:38,1)=nanmean(linefluorescenceseries_Chs,1)';
    
    fprintf(fid3,'%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',outputseg');
    fprintf(fid4,'%e\t %e\t %e\n',outputparametersseg');
    fprintf(fid6,'%e\t %e\t %e\t %e\n',curveincl_Chs');
    
    saveas(hh,[path2 '\' inputfilename(1:end-4) ' CCFs_seg.fig'])
    saveas(hI,[path2 '\' inputfilename(1:end-4) ' I_decomposed.fig'])
end

