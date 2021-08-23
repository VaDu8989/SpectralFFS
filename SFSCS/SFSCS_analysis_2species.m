clear all
close all
fclose all;

global movingaveragewindowsize blocksize spatialfilter depletioncorrection membranewidth

% Selection of GUIs and GUI parameters:
ACFselection=1; % 1 divide timetrace in ACFsegmentnumb segments and select or discard segment ACFs, default 1. If not desired, keep on 1 and set ACFsegmentnumb=2.
intensityselection=1; % 1(default): GUI to adjust exponential fit for bleaching correction (by removing short intenisty segments)
if ACFselection==1
    global x0 y0 fixed curveincl1 curveincl2 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesChCC  sigmas1curves sigmas2curves sigmascccurves corfit ItraceIICh1 ItraceIICh2 goon
end
if intensityselection==1
    global segmentsincl1 segmentsincl2 Ifull1 Ifull2 timeline1binned timeline2binned numberofsegments segmentlength correctionfitseg1 correctionfitseg2 ff1 ff2
end

ACFsegmentnumb=10; % Number of segments for segment-wise analysis
numberofsegments=20; % Number of intensity segments for exponential fit adjustment
movingaveragewindowsize=1000; % Number of blocks for lateral alignment of lines in kymograph
blocksize=movingaveragewindowsize;

% Acquisition parameters:
specchannelnumb=15; % specify number of spectral bins, e.g. 14 for A+Ch2, 15 for G+Y
spectralchannels=[495 504 513 522 531 540 548 557 566 575 584 593 602 611 620]; % For red FPs e.g. [575 584 593 602 611 620 629 637 646 655 664 673 682 691]; 
scantime_lines=403.20e-06; % Scan time (time to scan one line and move back the scanner). Minimum scan time on Zeiss LSM780 is 472.73e-06
pixeltime=1.23e-06; % Time to scan one pixel
timebreak=0; % Intervall in between subsequent scans (e.g. line_break_line_break...), default 0
scantime=scantime_lines+timebreak;
S=6.71; % Structure parameter as obtained from pFCS calibration of focal volume

% Analysis parameters:
spatialfilter=2.5; % Factor of stdev above and below mean mean membrane position (to select membrane pixels after Gaussian fitting)
depletioncorrection=1; % Apply bleaching correction? Default 1
binningwindow=100; % Intensity binning factor for binned time series (for display purposes only)
binningwindow_GUI=25; % Intensity binning factor for binned time series segments in GUI (for display purposes only)
backgroundcorrection=0; %1: background will be subtracted (as selected by ROI), 0: no subtraction
binning=0; % Timebinning?
binningfactor=1; % Binning factor for time binning (binningfactor=n -> 2^n binning in time)

% Fit parameters:
fixed=[false false true]; % Fix fit parameters (N, tau, S)? No: false, Yes: true
x0=[500,0.01,S]; % Initial conditions for ACF fits (N0,taud0,S)
y0=[5000,1,S]; % Initial conditions for CCF fits (N0,taud0,S) -> if no binding is expected, high N0 and taud0 allow the fit to converge in all cases
   
% Import of .tif files:
path= uigetdir; % Select directory containing spectral .tif files
files=dir([path '/*.tif']);
[namedata,remain]=strtok(files(1).name,'.');
name1file=files(1).name % loading first file
%         [namedata,remain]=strtok(files(end).name,'.');
%         name1file=files(end).name
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
rect=getrect; % Make rectangular selection
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
    fprintf('Window size not valid!\n')
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
savefig(emspec,[path2 '\' inputfilename(1:end-4) '_Flpectrum.fig'])
specoutput=zeros(specchannelnumb,2);
specoutput(:,1)=spectralchannels';
specoutput(:,2)=emissionspec';
fidspec=fopen([path2 '\' inputfilename(1:end-4) '_flspectrum.txt'],'a');
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
    Iinv=zeros(1,specchannelnumb);
    Ittt=mean(speclinefluorescenceseries,1);
    Iinv(Ittt>0)=1./Ittt(Ittt>0);
    D=diag(Iinv);
    w=(spectral_patterns*D*spectral_patterns')\spectral_patterns*D;
    for ttt=1:length(linefluorescenceseries)
        linefluorescenceseries_GFP(ttt)=sum(w(1,:).*speclinefluorescenceseries(ttt,:));
        linefluorescenceseries_YFP(ttt)=sum(w(2,:).*speclinefluorescenceseries(ttt,:));
    end
    lastbinlength=mod(length(linefluorescenceseries_GFP),binningwindow)+binningwindow;
    linefluorescenceseriesbinned_GFP=nanmean(reshape(linefluorescenceseries_GFP(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_GFP(1:end-lastbinlength))/binningwindow),1);
    linefluorescenceseriesbinned_GFP=[linefluorescenceseriesbinned_GFP nanmean(linefluorescenceseries_GFP(end-lastbinlength+1:end))];
    timeline_binned_GFP=1:1:length(linefluorescenceseriesbinned_GFP);
    timeline_binned_GFP=timeline_binned_GFP*binningwindow*scantime;
    timeline_binned_YFP=timeline_binned_GFP;
    lastbinlength=mod(length(linefluorescenceseries_GFP),binningwindow)+binningwindow;
    linefluorescenceseriesbinned_YFP=nanmean(reshape(linefluorescenceseries_YFP(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_GFP(1:end-lastbinlength))/binningwindow),1);
    linefluorescenceseriesbinned_YFP=[linefluorescenceseriesbinned_YFP nanmean(linefluorescenceseries_YFP(end-lastbinlength+1:end))];

    % Plot of decomposed time series:
    hI=figure('OuterPosition',[scrsz(3) scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2],'Name','Spectrally decomposed line fluorescence');
    plot(timeline_binned_GFP,linefluorescenceseriesbinned_GFP,'Color',[253/255 141/255 60/255]);
    hold on
    plot(timeline_binned_YFP,linefluorescenceseriesbinned_YFP,'Color',[151/255 25/255 28/255]);

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

        Iinv=zeros(1,specchannelnumb);
        Ittt=mean(speclinefluorescenceseries_binned,1);
        Iinv(Ittt>0)=1./Ittt(Ittt>0);
        D=diag(Iinv);
        w=(spectral_patterns*D*spectral_patterns')\spectral_patterns*D;
        for ttt=1:length(linefluorescenceseriesbinnedSNR)
            linefluorescenceseries_GFP(ttt)=sum(w(1,:).*speclinefluorescenceseries_binned(ttt,:));
            linefluorescenceseries_YFP(ttt)=sum(w(2,:).*speclinefluorescenceseries_binned(ttt,:));
        end
        lastbinlength=mod(length(linefluorescenceseries_GFP),binningwindow)+binningwindow;
        linefluorescenceseriesbinned_GFP=nanmean(reshape(linefluorescenceseries_GFP(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_GFP(1:end-lastbinlength))/binningwindow),1);
        linefluorescenceseriesbinned_GFP=[linefluorescenceseriesbinned_GFP nanmean(linefluorescenceseries_GFP(end-lastbinlength+1:end))];
        timeline_binned_GFP=1:1:length(linefluorescenceseriesbinned_GFP);
        timeline_binned_GFP=timeline_binned_GFP*binningwindow*scantime;
        timeline_binned_YFP=timeline_binned_GFP;
        lastbinlength=mod(length(linefluorescenceseries_GFP),binningwindow)+binningwindow;
        linefluorescenceseriesbinned_YFP=nanmean(reshape(linefluorescenceseries_YFP(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_GFP(1:end-lastbinlength))/binningwindow),1);
        linefluorescenceseriesbinned_YFP=[linefluorescenceseriesbinned_YFP nanmean(linefluorescenceseries_YFP(end-lastbinlength+1:end))];

        % Plot of decomposed time series:
        hI=figure('OuterPosition',[scrsz(3) scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2],'Name','Spectrally decomposed line fluorescence');
        plot(timeline_binned_GFP,linefluorescenceseriesbinned_GFP,'Color',[253/255 141/255 60/255]);
        hold on
        plot(timeline_binned_YFP,linefluorescenceseriesbinned_YFP,'Color',[151/255 25/255 28/255]);

        linefluorescenceseries=linefluorescenceseriesbinnedSNR';
        speclinefluorescenceseries=speclinefluorescenceseries_binned;
    end
end
    
% Plot of average spectra and calculation of spectral fractions:   
figure('Name','Fluorescent Spectra')
spectraldecomposfitfunc=@(x,t) x(1)*spectral_patterns(1,:)+x(2)*spectral_patterns(2,:);
x0spec=[0.5 0.5];
fixedspec=[false false];
[fitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixedspec,spectralchannels,emissionspec',spectraldecomposfitfunc,x0spec);
fractionGFP=fitparameters(1)
fractionYFP=fitparameters(2)
fractions=[fractionGFP fractionYFP];
avgspectrum_all2_fit=spectraldecomposfitfunc(fitparameters,spectralchannels);
fidfrac=fopen([path2  '\specfractions_G_Y.txt'],'a'); % adjust file name if needed, file to save spectral fractions
fprintf(fidfrac,inputfilename(1:end-4));
fprintf(fidfrac,'\t %e\t %e\n',fractions');
fidw=fopen([path2 '\' inputfilename(1:end-4) '_photonweights.txt'],'a'); % adjust path if needed, photon weights are saved as .txt
weightsoutput=zeros(length(spectralchannels),3);
weightsoutput(:,1)=spectralchannels';
weightsoutput(:,2:3)=w';
fprintf(fidw,'%e\t %e\t %e\n',weightsoutput');
    
plot(spectralchannels,emissionspec);
hold on
plot(spectralchannels,avgspectrum_all2_fit,'--');
xlabel('Channel [nm]')
ylabel('Norm.emission')
    
% Dialogue: If only one species is present in the measurement, the analysis
% can be interrupted here
str_fractions='Spectral fractions are: %.2f (G), %.2f (Y). Continue with FCS analysis?';
questionstring=sprintf(str_fractions,fractionGFP,fractionYFP);
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
bleachingfraction1=1-nanmean(linefluorescenceseries_GFP(end-250+1:end))/nanmean(linefluorescenceseries_GFP(1:250));
bleachingfraction2=1-nanmean(linefluorescenceseries_YFP(end-250+1:end))/nanmean(linefluorescenceseries_YFP(1:250));
timeline1=1:1:length(linefluorescenceseries_GFP);
timeline1=timeline1*scantime;
timeline2=timeline1;
line1fluorescenceseriesbinned=linefluorescenceseriesbinned_GFP;
line2fluorescenceseriesbinned=linefluorescenceseriesbinned_YFP;
timeline1_binned=timeline_binned_GFP;
timeline2_binned=timeline_binned_YFP;

% Depletion correction of time trace:
if intensityselection==1
    Ifull1=line1fluorescenceseriesbinned;
    Ifull2=line2fluorescenceseriesbinned;
    timeline1binned=timeline_binned_GFP;
    timeline2binned=timeline_binned_YFP;
    segmentlength=length(Ifull1)/numberofsegments;
    segmentsincl1=ones(1,numberofsegments);
    segmentsincl2=ones(1,numberofsegments);

    % GUI for selection of intensity segments:
    Intselection2colors
    goon=0;
    while goon==0
        pause(5)
    end
    line1fluorescenceseries_corrected=depletioncorrectionfunc(timeline1,line1fluorescenceseries,ff1);
    line2fluorescenceseries_corrected=depletioncorrectionfunc(timeline2,line2fluorescenceseries,ff2);
    lastbinlength=mod(length(line1fluorescenceseries),binningwindow)+binningwindow;
    line1fluorescenceseries_corrected_binned=mean(reshape(line1fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line1fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line1fluorescenceseries_corrected_binned=[line1fluorescenceseries_corrected_binned mean(line1fluorescenceseries_corrected(end-lastbinlength+1:end))];
    line2fluorescenceseries_corrected_binned=mean(reshape(line2fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line2fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line2fluorescenceseries_corrected_binned=[line2fluorescenceseries_corrected_binned mean(line2fluorescenceseries_corrected(end-lastbinlength+1:end))];
    figure('OuterPosition',[scrsz(3)/3 scrsz(4)/4 2*scrsz(3)/3 scrsz(4)/4],'Name','Corrected Membrane Fluorescence time series')
    subplot(2,1,1),plot(timeline1_binned,line1fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch1')
    subplot(2,1,2),plot(timeline2_binned,line2fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch2')      
else
    if depletioncorrection==1
    line1fluorescenceseries_corrected=depletioncorrectionfunc(timeline1,line1fluorescenceseries,f1); 
    line2fluorescenceseries_corrected=depletioncorrectionfunc(timeline2,line2fluorescenceseries,f2);
    lastbinlength=mod(length(line1fluorescenceseries),binningwindow)+binningwindow;
    line1fluorescenceseries_corrected_binned=mean(reshape(line1fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line1fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line1fluorescenceseries_corrected_binned=[line1fluorescenceseries_corrected_binned mean(line1fluorescenceseries_corrected(end-lastbinlength+1:end))];
    line2fluorescenceseries_corrected_binned=mean(reshape(line2fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line2fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line2fluorescenceseries_corrected_binned=[line2fluorescenceseries_corrected_binned mean(line2fluorescenceseries_corrected(end-lastbinlength+1:end))];
    figure('Name','Corrected Membrane Fluorescence time series')
    subplot(1,2,1),plot(timeline1_binned,line1fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch1')
    hold on
    subplot(1,2,2),plot(timeline2_binned,line2fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch2')
    end
 end

% Calculation of correlation functions:
fprintf('Calculating correlation using multiple tau...\n');
if depletioncorrection==1
    [tcorr1,fautocorr1,sigmas1]=autocorrFCSmultipletau(line1fluorescenceseries_corrected);
    [tcorr2,fautocorr2,sigmas2]=autocorrFCSmultipletau(line2fluorescenceseries_corrected);
        [tcorrcc,fcrosscorr,sigmascc]=crosscorrFCSmultipletau(line1fluorescenceseries_corrected,line2fluorescenceseries_corrected);
else
    [tcorr1,fautocorr1,sigmas1]=autocorrFCSmultipletau(line1fluorescenceseries);
    [tcorr2,fautocorr2,sigmas2]=autocorrFCSmultipletau(line2fluorescenceseries);
        [tcorrcc,fcrosscorr,sigmascc]=crosscorrFCSmultipletau(line1fluorescenceseries,line2fluorescenceseries);
end
tcorr1=tcorr1*scantime;
tcorr2=tcorr2*scantime;
tcorrcc=tcorrcc*scantime;
weights1=fautocorr1./sigmas1;
weights2=fautocorr2./sigmas2;
weightscc=fcrosscorr./sigmascc;
    
% Plot of CFs:
figure('OuterPosition',[scrsz(3) scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Correlation functions')
semilogx(tcorr1,fautocorr1,'g.')
hold on
semilogx(tcorr2,fautocorr2,'r.')
hold on
semilogx(tcorrcc,fcrosscorr,'b.')

if depletioncorrection==1
    line1fluorescenceseries=line1fluorescenceseries_corrected;
    line2fluorescenceseries=line2fluorescenceseries_corrected;
end

lastbinlength=mod(length(line1fluorescenceseries),binningwindow_GUI)+binningwindow_GUI;
line1fluorescenceseriesbinned_GUI=mean(reshape(line1fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line1fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line1fluorescenceseriesbinned_GUI=[line1fluorescenceseriesbinned_GUI mean(line1fluorescenceseries(end-lastbinlength+1:end))];
line2fluorescenceseriesbinned_GUI=mean(reshape(line2fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line2fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line2fluorescenceseriesbinned_GUI=[line2fluorescenceseriesbinned_GUI mean(line2fluorescenceseries(end-lastbinlength+1:end))];

% Segment-wise ACF analysis: A plot of CFs is generated which shows the CFs
% of all segments (overlayed)
if ACFselection==1
    segmentlength=length(line1fluorescenceseries)/ACFsegmentnumb;
    Isegmentlength_GUI=segmentlength/binningwindow_GUI;
    Itrace1=zeros(segmentlength,ACFsegmentnumb);
    Itrace2=Itrace1;
    Itrace1_GUI=zeros(Isegmentlength_GUI,ACFsegmentnumb);
    Itrace2_GUI=Itrace1_GUI;
    figure('OuterPosition',[5*scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','ACFs segments')
    for i=1:size(Itrace1,2)
        Itrace1(:,i)=line1fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)';
        Itrace2(:,i)=line2fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)';
        Itrace1_GUI(:,i)=line1fluorescenceseriesbinned_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI)';
        Itrace2_GUI(:,i)=line2fluorescenceseriesbinned_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI)';
        [tcorr1i,fcorr1i,sigmas1i]=autocorrFCSmultipletau(Itrace1(:,i));
        [tcorr2i,fcorr2i,sigmas2i]=autocorrFCSmultipletau(Itrace2(:,i));
        [tcorrcci,fcrosscorri,sigmascci]=crosscorrFCSmultipletau(Itrace1(:,i),Itrace2(:,i));
        tcorr1i=tcorr1i*scantime;
        tcorr2i=tcorr2i*scantime;
        tcorrcci=tcorrcci*scantime;
        correlationcurvesCh1(:,i)=fcorr1i';
        sigmas1curves(:,i)=sigmas1i';
        correlationcurvesCh2(:,i)=fcorr2i';
        sigmas2curves(:,i)=sigmas2i';
        correlationcurvesChCC(:,i)=fcrosscorri';
        sigmascccurves(:,i)=sigmascci';
        lbi=min(tcorr1i);
        ubi=max(tcorr1i);

        p1=subplot(2,2,1);semilogx(tcorr1i,fcorr1i,'Color',[26/255 150/255 65/255],'LineStyle','none','Marker','.','MarkerFaceColor',[26/255 150/255 65/255]);
        title('ACF 1')
        xlabel('\tau [s]')
        ylabel('Correlation')
        hold on
        p2=subplot(2,2,2);semilogx(tcorr2i,fcorr2i,'Color',[253/255 174/255 97/255],'LineStyle','none','Marker','.','MarkerFaceColor',[253/255 174/255 97/255]);
        title('ACF 2')
        xlabel('\tau [s]')
        ylabel('Correlation')
        hold on
        p3=subplot(2,2,3);semilogx(tcorrcci,fcrosscorri,'Color', [84/255 39/255 143/255],'LineStyle','none','Marker','.','MarkerFaceColor',[84/255 39/255 143/255]);
        title('CCF 12')
        xlabel('\tau [s]')
        ylabel('Correlation')
        %pause
        hold on
    end
    corfit1=zeros(length(tcorr1i),ACFsegmentnumb);
    corfit2=corfit1;
    corfitcc=corfit1;
    ItraceIICh1=Itrace1_GUI;
    ItraceIICh2=Itrace2_GUI;
    segmenttime1_GUI=scantime:scantime*binningwindow_GUI:scantime*segmentlength;
    segmenttime2_GUI=segmenttime1_GUI;
    ItraceIICh1=[segmenttime1_GUI' ItraceIICh1];
    ItraceIICh2=[segmenttime2_GUI' ItraceIICh2];

    correlationcurvesCh1=[tcorr1i' correlationcurvesCh1];
    size(correlationcurvesCh1)
    correlationcurvesCh2=[tcorr2i' correlationcurvesCh2];
    correlationcurvesChCC=[tcorrcci' correlationcurvesChCC];
    corfit1=[tcorr1i' corfit1];
    corfit2=[tcorr2i' corfit2];
    corfitcc=[tcorrcci' corfitcc];
    curveincl1=ones(size(corfit1,2)-1,1);
    curveincl2=curveincl1;

    % GUI for selection of segment ACFs:

    %Corrselection2Channelsindividuell
    Corrselection2Channelsindividuellpreview
    goon=0;
    while goon==0
        pause(5)
    end
    fprintf('Done!')
    for i=1:length(curveincl1)
        if curveincl1(i)==0
            line1fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=NaN;
        end
        if curveincl2(i)==0
            line2fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=NaN;
        end
    end

    % Calculation of final CFs:
    [tcorr1final,fcorr1final,sigmas1final]=autocorrFCSmultipletau(line1fluorescenceseries);
    [tcorr2final,fcorr2final,sigmas2final]=autocorrFCSmultipletau(line2fluorescenceseries);
    [tcorrccfinal,fcrosscorrfinal,sigmasccfinal]=crosscorrFCSmultipletau(line1fluorescenceseries,line2fluorescenceseries);
    tcorr1final=tcorr1final*scantime;
    tcorr2final=tcorr2final*scantime;
    tcorrccfinal=tcorrccfinal*scantime;
    weights1final=abs(fcorr1final./sigmas1final);
    weights2final=abs(fcorr2final./sigmas2final);
    weightsccfinal=abs(fcrosscorrfinal./sigmasccfinal);

    % Average of segments:
    fcorr1segfinal=zeros(length(tcorr1i),1);
    sigmas1segfinal=fcorr1segfinal;
    weight1seg=fcorr1segfinal;
    fcorr2segfinal=zeros(length(tcorr2i),1);
    sigmas2segfinal=fcorr2segfinal;
    weight2seg=fcorr2segfinal;

    fcrosscorrsegfinal=zeros(length(tcorrcci),1);
    sigmasccsegfinal=fcrosscorrsegfinal;
    weightsccseg=fcrosscorrsegfinal;
    for i=1:length(curveincl1);
        if curveincl1(i)==1
             fcorr1segfinal=fcorr1segfinal+correlationcurvesCh1(:,1+i);
             sigmas1segfinal=sigmas1segfinal+sigmas1curves(:,i);
        end
        if curveincl2(i)==1
            fcorr2segfinal=fcorr2segfinal+correlationcurvesCh2(:,1+i);
            sigmas2segfinal=sigmas2segfinal+sigmas2curves(:,i);
        end
        if (curveincl1(i)+curveincl2(i))==2
        fcrosscorrsegfinal=fcrosscorrsegfinal+correlationcurvesChCC(:,1+i);
        sigmasccsegfinal=sigmasccsegfinal+sigmascccurves(:,i);
        end
    end
    fcorr1segfinal=fcorr1segfinal/sum(curveincl1);
    fcorr2segfinal=fcorr2segfinal/sum(curveincl2);
    fcrosscorrsegfinal=fcrosscorrsegfinal/sum((curveincl1+curveincl2)==2);
    sigmas1segfinal=sigmas1segfinal/sum(curveincl1);
    sigmas2segfinal=sigmas2segfinal/sum(curveincl2);
    sigmasccsegfinal=sigmasccsegfinal/sum((curveincl1+curveincl2)==2);
    weights1segfinal=abs(fcorr1segfinal./sigmas1segfinal);
    tcorr1segfinal=tcorr1i;
    weights2segfinal=abs(fcorr2segfinal./sigmas2segfinal);
    tcorr2segfinal=tcorr2i;
    weightsccsegfinal=abs(fcrosscorrsegfinal./sigmasccsegfinal);
    tcorrccsegfinal=tcorrcci;

    % Fit of final CFs with diffusion model:
    fprintf('Final Fitting...\n');
    lb=scantime;
    lbcc=lb;
    lbseg=lb;
    lbccseg=lb;
    ub=max(tcorr1final);
    ubcc=ub;
    tfit=tcorr1final;
    ubseg=max(tcorr1segfinal);
    ubccseg=ubseg;
    tcorr1fit=tcorr1final(1:length(tfit));
    fcorr1fit=fcorr1final(1:length(tfit));
    weights1fit=weights1final(1:length(tfit));
    weights1fit(weights1fit==0)=10^-4;

    tcorr2fit=tcorr2final(1:length(tfit));
    fcorr2fit=fcorr2final(1:length(tfit));
    weights2fit=weights2final(1:length(tfit));
    weights2fit(weights2fit==0)=10^-4;

    tcorrccfit=tcorrccfinal(1:length(tfit));
    fcrosscorrfit=fcrosscorrfinal(1:length(tfit));
    weightsccfit=weightsccfinal(1:length(tfit));
    weightsccfit(weightsccfit==0)=10^-4;
    tcorrccsegfinalfit=tcorrccsegfinal(1:length(tcorr1segfinal));
    fcrosscorrsegfinalfit=fcrosscorrsegfinal(1:length(tcorr1segfinal));
    weightsccsegfinalfit=weightsccsegfinal(1:length(tcorr1segfinal));

    lsautofitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
    flscrossfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);

    % Fit of whole time series:
    [N1final,taud1final,Sfit1final,CI1final,fitcurve1final,residuals1final] = autocorrfit2Ddiff(tcorr1fit,fcorr1fit,lb,ub,weights1fit,lsautofitfunc,x0,fixed);
    [N2final,taud2final,Sfit2final,CI2final,fitcurve2final,residuals2final] = autocorrfit2Ddiff(tcorr2fit,fcorr2fit,lb,ub,weights2fit,lsautofitfunc,x0,fixed);
    [Nccfinal,tauccfinal,Sfitccfinal,CIccfinal,fitcurveccfinal,residualsccfinal] =autocorrfit2Ddiff(tcorrccfit,fcrosscorrfit,lbcc,ubcc,weightsccfit,lsautofitfunc,y0,fixed);
    
    % Fit of segment-averaged time series:
    [N1segfinal,taud1segfinal,Sfit1segfinal,CI1segfinal,fitcurve1segfinal,residuals1segfinal] = autocorrfit2Ddiff(tcorr1segfinal,fcorr1segfinal',lbseg,ubseg,weights1segfinal',lsautofitfunc,x0,fixed);
    [N2segfinal,taud2segfinal,Sfit2segfinal,CI2segfinal,fitcurve2segfinal,residuals2segfinal] = autocorrfit2Ddiff(tcorr2segfinal,fcorr2segfinal',lbseg,ubseg,weights2segfinal',lsautofitfunc,x0,fixed);
    [Nccsegfinal,tauccsegfinal,Sfitccsegfinal,CIccsegfinal,fitcurveccsegfinal,residualsccsegfinal] =autocorrfit2Ddiff(tcorrccsegfinalfit,fcrosscorrsegfinalfit',lbccseg,ubccseg,weightsccsegfinalfit',lsautofitfunc,y0,fixed);

    % Relative cross-correlation amplitude:
    relCC=max([N1final/Nccfinal N2final/Nccfinal]);
    relCCsegnew=max([N1segfinal/Nccsegfinal N2segfinal/Nccsegfinal]);
   
    % Confidencen intervals:
    us1=0.5*(CI1final(1:3,2)-CI1final(1:3,1));
    us2=0.5*(CI2final(1:3,2)-CI2final(1:3,1));
    uscc=0.5*(CIccfinal(1:3,2)-CIccfinal(1:3,1));
    
    us1seg=0.5*(CI1segfinal(1:3,2)-CI1segfinal(1:3,1));
    us2seg=0.5*(CI2segfinal(1:3,2)-CI2segfinal(1:3,1));
    usccseg=0.5*(CIccsegfinal(1:3,2)-CIccsegfinal(1:3,1));
    
    % Figure: Plot of data and fit
    % Full:
    h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curve');
    positionvector1=[0.1 0.35 0.8 0.55];
    positionvector2=[0.1 0.1 0.8 0.15];
    subplot('Position',positionvector1),semilogx(tcorr1fit,fitcurve1final,'-g',tcorr2fit,fitcurve2final,'-r',tcorrccfit,fitcurveccfinal,'-b')
    hold on
    subplot('Position',positionvector1),semilogx(tcorr1fit,fcorr1fit,'gs',tcorr2fit,fcorr2fit,'rd',tcorrccfit,fcrosscorrfit,'bx')
    xlabel('time')
    ylabel('autocorrelation')
    subplot('Position',positionvector2),semilogx(tcorr1fit,residuals1final,'gs',tcorr2fit,residuals2final,'rd',tcorrccfit,residualsccfinal,'bx')
    xlabel('time')
    ylabel('residuals')
    hold on
    subplot('Position',positionvector2),semilogx(tcorr1fit,zeros(size(tcorr1fit)),'-k');
    textN1=['N_1 = %.1f ' char(177) ' %.1f'];
    strtextN1=sprintf(textN1,N1final,us1(1));
    textN2=['N_2 = %.1f ' char(177) ' %.1f'];
    strtextN2=sprintf(textN2,N2final,us2(1));
    textrelcc='rel.cc. = %.2f ';
    strtextrelcc=sprintf(textrelcc,relCC);
    texttaud1=[' = %.2f ' char(177) ' %.2f'];
    strtexttaud1=['\tau_1' sprintf(texttaud1,taud1final,us1(2)) ' s'] ;
    texttaud2=[' = %.2f ' char(177) ' %.2f'];
    strtexttaud2=['\tau_2' sprintf(texttaud2,taud2final,us2(2)) ' s'] ;
    texttaudcc=[' = %.2f ' char(177) ' %.2f'];
    strtexttaudcc=['\tau_{cc}' sprintf(texttaudcc,tauccfinal,uscc(2)) ' s'] ;
    dim = [0.68 0.58 0.3 0.3];
    str = {strtextN1,strtextN2,strtextrelcc,strtexttaud1,strtexttaud2,strtexttaudcc};
    t=annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(t,'FontSize',18)
    
    % Segments:
    hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curve');
    positionvector1=[0.1 0.35 0.8 0.55];
    positionvector2=[0.1 0.1 0.8 0.15];
    subplot('Position',positionvector1),semilogx(tcorr1segfinal,fitcurve1segfinal,'-g',tcorr2segfinal,fitcurve2segfinal,'-r',tcorrccsegfinalfit,fitcurveccsegfinal,'-b')
    hold on
    subplot('Position',positionvector1),semilogx(tcorr1segfinal,fcorr1segfinal,'gs',tcorr2segfinal,fcorr2segfinal,'rd',tcorrccsegfinalfit,fcrosscorrsegfinalfit,'bx')
    xlabel('time')
    ylabel('autocorrelation')
    subplot('Position',positionvector2),semilogx(tcorr1segfinal,residuals1segfinal,'gs',tcorr2segfinal,residuals2segfinal,'rd',tcorrccsegfinalfit,residualsccsegfinal,'bx')
    xlabel('time')
    ylabel('residuals')
    hold on
    subplot('Position',positionvector2),semilogx(tcorr1segfinal,zeros(size(tcorr1segfinal)),'-k');
    textN1seg=['N_1 = %.1f ' char(177) ' %.1f'];
    strtextN1seg=sprintf(textN1seg,N1segfinal,us1seg(1));
    textN2seg=['N_2 = %.1f ' char(177) ' %.1f'];
    strtextN2seg=sprintf(textN2seg,N2segfinal,us2seg(1));
    textrelccseg='rel.cc. = %.2f ';
    strtextrelccseg=sprintf(textrelccseg,relCCsegnew);
    texttaud1seg=[' = %.2f ' char(177) ' %.2f'];
    strtexttaud1seg=['\tau_1' sprintf(texttaud1seg,taud1segfinal,us1seg(2)) ' s'] ;
    texttaud2seg=[' = %.2f ' char(177) ' %.2f'];
    strtexttaud2seg=['\tau_2' sprintf(texttaud2seg,taud2segfinal,us2seg(2)) ' s'] ;
    texttaudccseg=[' = %.2f ' char(177) ' %.2f'];
    strtexttaudccseg=['\tau_{cc}' sprintf(texttaudccseg,tauccsegfinal,usccseg(2)) ' s'] ;
    dim = [0.68 0.58 0.3 0.3];
    strseg = {strtextN1seg,strtextN2seg,strtextrelccseg,strtexttaud1seg,strtexttaud2seg,strtexttaudccseg};
    t=annotation('textbox',dim,'String',strseg,'FitBoxToText','on');
    set(t,'FontSize',18)
   
    % Save output (fit parameters, ACFs etc)
    fid1=fopen([path2 '\' inputfilename(1:end-4) '_final_ACF.txt'],'a');
    fid2=fopen([path2 '\' inputfilename(1:end-4) '_final_fitparameters.txt'],'a'); 
    fid3=fopen([path2 '\' inputfilename(1:end-4) '_final_ACFseg.txt'],'a');
    fid4=fopen([path2 '\' inputfilename(1:end-4) '_final_fitparameters_seg.txt'],'a');
    fid5=fopen([path2 '\' inputfilename(1:end-4) 'channel_Is.mat'],'a'); 
    fid6=fopen([path2 '\' inputfilename(1:end-4) '_final_curveincls.txt'],'a');
    
    linechannelfluorescenceseries=[line1fluorescenceseries' line2fluorescenceseries'];
    save([path2 '\' inputfilename(1:end-4) 'channel_Is.mat'],'linechannelfluorescenceseries')
    
    % Segment analysis:
    outputseg=zeros(length(tcorr1segfinal),9);    
    outputseg(:,1)=tcorr1segfinal;
    outputseg(:,2)=fcorr1segfinal;
    outputseg(:,3)=weights1segfinal;
    outputseg(:,4)=tcorr2segfinal;
    outputseg(:,5)=fcorr2segfinal;
    outputseg(:,6)=weights2segfinal;
    outputseg(:,7)=tcorrccsegfinal;
    outputseg(:,8)=fcrosscorrsegfinal;
    outputseg(:,9)=weightsccsegfinal;
    
    outputparametersseg=zeros(13,3);
    outputparametersseg(1:3,1)=[N1segfinal;taud1segfinal;Sfit1segfinal];
    outputparametersseg(1:3,2:end)=CI1segfinal;
    outputparametersseg(4:6,1)=[N2segfinal;taud2segfinal;Sfit2segfinal];
    outputparametersseg(4:6,2:end)=CI2segfinal;
    outputparametersseg(7:9,1)=[Nccsegfinal;tauccsegfinal;Sfitccsegfinal];
    outputparametersseg(7:9,2:end)=CIccsegfinal(1:3,:);
    outputparametersseg(10,1)=bleachingfraction1;
    outputparametersseg(11,1)=bleachingfraction2;
    outputparametersseg(12,1)=nanmean(line1fluorescenceseries);
    outputparametersseg(13,1)=nanmean(line2fluorescenceseries);
    
    fprintf(fid3,'%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',outputseg');
    fprintf(fid4,'%e\t %e\t %e\n',outputparametersseg');
    
    % Full analysis:
    output=zeros(length(tcorr2),9);    
    output(:,1)=tcorr1final;
    output(:,2)=fcorr1final;
    output(:,3)=weights1final;
    output(:,4)=tcorr2final;
    output(:,5)=fcorr2final;
    output(:,6)=weights2final;
    output(:,7)=tcorrccfinal;
    output(:,8)=fcrosscorrfinal;
    output(:,9)=weightsccfinal;
    
    outputparameters=zeros(13,3);
    outputparameters(1:3,1)=[N1final;taud1final;Sfit1final];
    outputparameters(1:3,2:end)=CI1final;
    outputparameters(4:6,1)=[N2final;taud2final;Sfit2final];
    outputparameters(4:6,2:end)=CI2final;
    outputparameters(7:9,1)=[Nccfinal;tauccfinal;Sfitccfinal];
    outputparameters(7:9,2:end)=CIccfinal(1:3,:);
    outputparameters(10,1)=bleachingfraction1;
    outputparameters(11,1)=bleachingfraction2;
    outputparameters(12,1)=nanmean(line1fluorescenceseries);
    outputparameters(13,1)=nanmean(line2fluorescenceseries);
    
    fprintf(fid1,'%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output');
    fprintf(fid2,'%e\t %e\t %e\n',outputparameters');
    
    curveincls=[curveincl1 curveincl2];
    fprintf(fid6,'%e\t %e\n',curveincls');
    
    % Save figures:
    saveas(hh,[path2 '\' inputfilename(1:end-4) ' CCFs_seg.fig'])
    saveas(h,[path2 '\' inputfilename(1:end-4) ' CCFs_full.fig'])
    saveas(hI,[path2 '\' inputfilename(1:end-4) ' I_decomposed.fig'])
end

