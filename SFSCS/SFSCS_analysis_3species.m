clear all
close all
fclose all;

global movingaveragewindowsize blocksize spatialfilter depletioncorrection membranewidth

% Selection of GUIs and GUI parameters:
ACFselection=1; % 1 (default): GUI to divide timetrace in ACFsegmentnumb segments and select or discard segment ACFs/CCFs. Final CFs are the average of all kept segments
intensityselection=1; % 1(default): GUI to adjust exponential fit for bleaching correction (by removing short intenisty segments)
if ACFselection==1
    global x0 y0 fixed curveincl1 curveincl2 curveincl3 correlationcurvesCh1 correlationcurvesCh2 correlationcurvesCh3 correlationcurvesChCC12 correlationcurvesChCC13 correlationcurvesChCC23  sigmas1curves sigmas2curves sigmas3curves sigmascc12curves sigmascc13curves sigmascc23curves corfit ItraceIICh1 ItraceIICh2 ItraceIICh3 goon
end
if intensityselection==1
    global segmentsincl1 segmentsincl2 segmentsincl3 Ifull1 Ifull2 Ifull3 timeline1binned timeline2binned timeline3binned numberofsegments segmentlength correctionfitseg1 correctionfitseg2 correctionfitseg3 ff1 ff2 ff3
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
path= uigetdir; % Select directory containing spectral .tif files
files=dir([path '/*.tif']);
[namedata,remain]=strtok(files(1).name,'.');
name1file=files(1).name % loading first file
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
    Iinv=zeros(1,specchannelnumb);
    Ittt=mean(speclinefluorescenceseries,1);
    Iinv(Ittt>0)=1./Ittt(Ittt>0);
    D=diag(Iinv);
    w=(spectral_patterns*D*spectral_patterns')\spectral_patterns*D;
    for ttt=1:length(linefluorescenceseries)
        linefluorescenceseries_GFP(ttt)=sum(w(1,:).*speclinefluorescenceseries(ttt,:));
        linefluorescenceseries_YFP(ttt)=sum(w(2,:).*speclinefluorescenceseries(ttt,:));
        linefluorescenceseries_Ch(ttt)=sum(w(3,:).*speclinefluorescenceseries(ttt,:));
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
    linefluorescenceseriesbinned_Ch=nanmean(reshape(linefluorescenceseries_Ch(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_Ch(1:end-lastbinlength))/binningwindow),1);
    linefluorescenceseriesbinned_Ch=[linefluorescenceseriesbinned_Ch nanmean(linefluorescenceseries_Ch(end-lastbinlength+1:end))];

    % Plot of decomposed time series:
    hI=figure('OuterPosition',[scrsz(3) scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2],'Name','Spectrally decomposed line fluorescence');
    plot(timeline_binned_GFP,linefluorescenceseriesbinned_GFP,'Color',[26/255 150/255 65/255]);
    hold on
    plot(timeline_binned_YFP,linefluorescenceseriesbinned_YFP,'Color',[253/255 174/255 97/255]);
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
        linefluorescenceseries_Ch=zeros(size(linefluorescenceseriesbinnedSNR));

        Iinv=zeros(1,specchannelnumb);
        Ittt=mean(speclinefluorescenceseries_binned,1);
        Iinv(Ittt>0)=1./Ittt(Ittt>0);
        D=diag(Iinv);
        w=(spectral_patterns*D*spectral_patterns')\spectral_patterns*D;
        
        for ttt=1:length(linefluorescenceseriesbinnedSNR)
            linefluorescenceseries_GFP(ttt)=sum(w(1,:).*speclinefluorescenceseries_binned(ttt,:));
            linefluorescenceseries_YFP(ttt)=sum(w(2,:).*speclinefluorescenceseries_binned(ttt,:));
            linefluorescenceseries_Ch(ttt)=sum(w(3,:).*speclinefluorescenceseries_binned(ttt,:));
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
        linefluorescenceseriesbinned_Ch=nanmean(reshape(linefluorescenceseries_Ch(1:end-lastbinlength),binningwindow,length(linefluorescenceseries_Ch(1:end-lastbinlength))/binningwindow),1);
        linefluorescenceseriesbinned_Ch=[linefluorescenceseriesbinned_Ch nanmean(linefluorescenceseries_Ch(end-lastbinlength+1:end))];

        % Plot of decomposed time series:
        hI=figure('OuterPosition',[scrsz(3) scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2],'Name','Spectrally decomposed line fluorescence');
        plot(timeline_binned_GFP,linefluorescenceseriesbinned_GFP,'Color',[26/255 150/255 65/255]);
        hold on
        plot(timeline_binned_YFP,linefluorescenceseriesbinned_YFP,'Color',[253/255 174/255 97/255]);
        hold on
        plot(timeline_binned_Ch,linefluorescenceseriesbinned_Ch,'Color',[151/255 25/255 28/255]);

        linefluorescenceseries=linefluorescenceseriesbinnedSNR';
        speclinefluorescenceseries=speclinefluorescenceseries_binned;
    end
end

% Plot of average spectra and calculation of spectral fractions: 
figure('Name','Fluorescent Spectra')
spectraldecomposfitfunc=@(x,t) x(1)*spectral_patterns(1,:)+x(2)*spectral_patterns(2,:)+x(3)*spectral_patterns(3,:);
x0spec=[0.5 0.5 0.5];
fixedspec=[false false false];
[fitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixedspec,spectralchannels,emissionspec',spectraldecomposfitfunc,x0spec);
fractionGFP=fitparameters(1)
fractionYFP=fitparameters(2)
fractionCh=fitparameters(3)
fractions=[fractionGFP fractionYFP fractionCh];
avgspectrum_all3_fit=spectraldecomposfitfunc(fitparameters,spectralchannels);
fidfrac=fopen([path2  '\specfractions_GFP_YFP_Ch2.txt'],'a'); % adjust file name if needed, file to save spectral fractions
fprintf(fidfrac,inputfilename(1:end-4));
fprintf(fidfrac,'\t %e\t %e\t %e\n',fractions');

fidw=fopen([path2 '\' inputfilename(1:end-4) '_photonweights.txt'],'a'); % adjust file name if needed, photon weights are saved as .txt
weightsoutput=zeros(length(spectralchannels),4);
weightsoutput(:,1)=spectralchannels';
weightsoutput(:,2:4)=w';
fprintf(fidw,'%e\t %e\t %e\t %e\n',weightsoutput');
    
plot(spectralchannels,emissionspec);
hold on
plot(spectralchannels,avgspectrum_all3_fit,'--');
xlabel('Channel [nm]')
ylabel('Norm.emission')
    
% Dialogue: If only one species is present in the measurement, the analysis
% can be interrupted here
str_fractions='Spectral fractions are: %.2f (G), %.2f (Y), %.2f (Ch2). Continue with FCS analysis?';
questionstring=sprintf(str_fractions,fractionGFP,fractionYFP,fractionCh);
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
line3fluorescenceseries=linefluorescenceseries_Ch;
bleachingfraction1=1-nanmean(linefluorescenceseries_GFP(end-250+1:end))/nanmean(linefluorescenceseries_GFP(1:250));
bleachingfraction2=1-nanmean(linefluorescenceseries_YFP(end-250+1:end))/nanmean(linefluorescenceseries_YFP(1:250));
bleachingfraction3=1-nanmean(linefluorescenceseries_Ch(end-250+1:end))/nanmean(linefluorescenceseries_Ch(1:250));
timeline1=1:1:length(linefluorescenceseries_GFP);
timeline1=timeline1*scantime;
timeline2=timeline1;
timeline3=timeline1;
line1fluorescenceseriesbinned=linefluorescenceseriesbinned_GFP;
line2fluorescenceseriesbinned=linefluorescenceseriesbinned_YFP;
line3fluorescenceseriesbinned=linefluorescenceseriesbinned_Ch;

timeline1_binned=timeline_binned_GFP;
timeline2_binned=timeline_binned_YFP;
timeline3_binned=timeline_binned_Ch;

% Depletion correction of time trace:
if intensityselection==1
    Ifull1=line1fluorescenceseriesbinned;
    Ifull2=line2fluorescenceseriesbinned;
    Ifull3=line3fluorescenceseriesbinned;
    timeline1binned=timeline_binned_GFP;
    timeline2binned=timeline_binned_YFP;
    timeline3binned=timeline_binned_Ch;
    segmentlength=length(Ifull1)/numberofsegments;
    segmentsincl1=ones(1,numberofsegments);
    segmentsincl2=ones(1,numberofsegments);
    segmentsincl3=ones(1,numberofsegments);
        
    % GUI for selection of intensity segments:
    Intselection3colors
    goon=0;
    while goon==0
        pause(5)
    end
    line1fluorescenceseries_corrected=depletioncorrectionfunc(timeline1,line1fluorescenceseries,ff1);
    line2fluorescenceseries_corrected=depletioncorrectionfunc(timeline2,line2fluorescenceseries,ff2);
    line3fluorescenceseries_corrected=depletioncorrectionfunc(timeline3,line3fluorescenceseries,ff3);
    lastbinlength=mod(length(line1fluorescenceseries),binningwindow)+binningwindow;
    line1fluorescenceseries_corrected_binned=mean(reshape(line1fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line1fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line1fluorescenceseries_corrected_binned=[line1fluorescenceseries_corrected_binned mean(line1fluorescenceseries_corrected(end-lastbinlength+1:end))];
    line2fluorescenceseries_corrected_binned=mean(reshape(line2fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line2fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line2fluorescenceseries_corrected_binned=[line2fluorescenceseries_corrected_binned mean(line2fluorescenceseries_corrected(end-lastbinlength+1:end))];
    line3fluorescenceseries_corrected_binned=mean(reshape(line3fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line3fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line3fluorescenceseries_corrected_binned=[line3fluorescenceseries_corrected_binned mean(line3fluorescenceseries_corrected(end-lastbinlength+1:end))];
    figure('OuterPosition',[scrsz(3)/3 scrsz(4)/4 2*scrsz(3)/3 scrsz(4)/4],'Name','Corrected Membrane Fluorescence time series')
    subplot(3,1,1),plot(timeline1_binned,line1fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch1')
    subplot(3,1,2),plot(timeline2_binned,line2fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch2')
    subplot(3,1,3),plot(timeline3_binned,line3fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch3')
else
    if depletioncorrection==1
    line1fluorescenceseries_corrected=depletioncorrectionfunc(timeline1,line1fluorescenceseries,f1);
    line2fluorescenceseries_corrected=depletioncorrectionfunc(timeline2,line2fluorescenceseries,f2);
    line3fluorescenceseries_corrected=depletioncorrectionfunc(timeline3,line3fluorescenceseries,f3);
    lastbinlength=mod(length(line1fluorescenceseries),binningwindow)+binningwindow;
    line1fluorescenceseries_corrected_binned=mean(reshape(line1fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line1fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line1fluorescenceseries_corrected_binned=[line1fluorescenceseries_corrected_binned mean(line1fluorescenceseries_corrected(end-lastbinlength+1:end))];
    line2fluorescenceseries_corrected_binned=mean(reshape(line2fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line2fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line2fluorescenceseries_corrected_binned=[line2fluorescenceseries_corrected_binned mean(line2fluorescenceseries_corrected(end-lastbinlength+1:end))];
    line3fluorescenceseries_corrected_binned=mean(reshape(line3fluorescenceseries_corrected(1:end-lastbinlength),binningwindow,length(line3fluorescenceseries_corrected(1:end-lastbinlength))/binningwindow),1);
    line3fluorescenceseries_corrected_binned=[line3fluorescenceseries_corrected_binned mean(line3fluorescenceseries_corrected(end-lastbinlength+1:end))];
    figure('Name','Corrected Membrane Fluorescence time series')
    subplot(1,3,1),plot(timeline1_binned,line1fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch1')
    hold on
    subplot(1,3,2),plot(timeline2_binned,line2fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch2')
    subplot(1,3,3),plot(timeline3_binned,line3fluorescenceseries_corrected_binned);
    xlabel('time (linenumber)')
    ylabel('line fluorescence Ch3')
    end
end

% Calculation of correlation functions:
fprintf('Calculating correlation using multiple tau...\n');
if depletioncorrection==1
    [tcorr1,fautocorr1,sigmas1]=autocorrFCSmultipletau(line1fluorescenceseries_corrected);
    [tcorr2,fautocorr2,sigmas2]=autocorrFCSmultipletau(line2fluorescenceseries_corrected);
    [tcorr3,fautocorr3,sigmas3]=autocorrFCSmultipletau(line3fluorescenceseries_corrected);
    [tcorrcc12,fcrosscorr12,sigmascc12]=crosscorrFCSmultipletau(line1fluorescenceseries_corrected,line2fluorescenceseries_corrected);
    [tcorrcc13,fcrosscorr13,sigmascc13]=crosscorrFCSmultipletau(line1fluorescenceseries_corrected,line3fluorescenceseries_corrected);
    [tcorrcc23,fcrosscorr23,sigmascc23]=crosscorrFCSmultipletau(line2fluorescenceseries_corrected,line3fluorescenceseries_corrected);
else
    [tcorr1,fautocorr1,sigmas1]=autocorrFCSmultipletau(line1fluorescenceseries);
    [tcorr2,fautocorr2,sigmas2]=autocorrFCSmultipletau(line2fluorescenceseries);
    [tcorr3,fautocorr3,sigmas3]=autocorrFCSmultipletau(line3fluorescenceseries);
    [tcorrcc12,fcrosscorr12,sigmascc12]=crosscorrFCSmultipletau(line1fluorescenceseries,line2fluorescenceseries);
    [tcorrcc13,fcrosscorr13,sigmascc13]=crosscorrFCSmultipletau(line1fluorescenceseries,line3fluorescenceseries);
    [tcorrcc23,fcrosscorr23,sigmascc23]=crosscorrFCSmultipletau(line2fluorescenceseries,line3fluorescenceseries);
end
tcorr=tcorr1*scantime; %all tcorr's are equal anyways
tcorrcc=tcorrcc12*scantime; %In this case, all tcorr's are equal anyways 
weights1=fautocorr1./sigmas1;
weights2=fautocorr2./sigmas2;
weights3=fautocorr3./sigmas3;
weightscc12=fcrosscorr12./sigmascc12;
weightscc13=fcrosscorr13./sigmascc13;
weightscc23=fcrosscorr23./sigmascc23;

% Plot of CFs:
figure('OuterPosition',[scrsz(3) scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Correlation functions')
semilogx(tcorr,fautocorr1,'.','Color',[26/255 150/255 65/255])
hold on
semilogx(tcorr,fautocorr2,'.','Color',[253/255 174/255 97/255])
hold on
semilogx(tcorr,fautocorr3,'.','Color',[151/255 25/255 28/255])
hold on
semilogx(tcorrcc,fcrosscorr12,'.','Color', [84/255 39/255 143/255])
hold on
semilogx(tcorrcc,fcrosscorr13,'d','Color', [43/255 140/255 190/255])
hold on
semilogx(tcorrcc,fcrosscorr23,'x','Color', [150/255 150/255 150/255])

if depletioncorrection==1
    line1fluorescenceseries=line1fluorescenceseries_corrected;
    line2fluorescenceseries=line2fluorescenceseries_corrected;
    line3fluorescenceseries=line3fluorescenceseries_corrected;
end

lastbinlength=mod(length(line1fluorescenceseries),binningwindow_GUI)+binningwindow_GUI;
line1fluorescenceseriesbinned_GUI=mean(reshape(line1fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line1fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line1fluorescenceseriesbinned_GUI=[line1fluorescenceseriesbinned_GUI mean(line1fluorescenceseries(end-lastbinlength+1:end))];
line2fluorescenceseriesbinned_GUI=mean(reshape(line2fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line2fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line2fluorescenceseriesbinned_GUI=[line2fluorescenceseriesbinned_GUI mean(line2fluorescenceseries(end-lastbinlength+1:end))];
line3fluorescenceseriesbinned_GUI=mean(reshape(line3fluorescenceseries(1:end-lastbinlength),binningwindow_GUI,length(line3fluorescenceseries(1:end-lastbinlength))/binningwindow_GUI),1);
line3fluorescenceseriesbinned_GUI=[line3fluorescenceseriesbinned_GUI mean(line3fluorescenceseries(end-lastbinlength+1:end))];


% Segment-wise ACF analysis: A plot of CFs is generated which shows the CFs
% of all segments (overlayed)
if ACFselection==1
    segmentlength=length(line1fluorescenceseries)/ACFsegmentnumb;
    Isegmentlength_GUI=segmentlength/binningwindow_GUI;
    Itrace1=zeros(segmentlength,ACFsegmentnumb);
    Itrace2=Itrace1;
    Itrace3=Itrace1;
    Itrace1_GUI=zeros(Isegmentlength_GUI,ACFsegmentnumb);
    Itrace2_GUI=Itrace1_GUI;
    Itrace3_GUI=Itrace1_GUI;
    figure('OuterPosition',[5*scrsz(3)/3 scrsz(4)/4 scrsz(3)/3 3*scrsz(4)/4],'Name','ACFs segments')
    for i=1:size(Itrace1,2)
        Itrace1(:,i)=line1fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)';
        Itrace2(:,i)=line2fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)';
        Itrace3(:,i)=line3fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)';
        Itrace1_GUI(:,i)=line1fluorescenceseriesbinned_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI)';
        Itrace2_GUI(:,i)=line2fluorescenceseriesbinned_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI)';
        Itrace3_GUI(:,i)=line3fluorescenceseriesbinned_GUI(Isegmentlength_GUI*(i-1)+1:i*Isegmentlength_GUI)';
        [tcorr1i,fcorr1i,sigmas1i]=autocorrFCSmultipletau(Itrace1(:,i));
        [tcorr2i,fcorr2i,sigmas2i]=autocorrFCSmultipletau(Itrace2(:,i));
        [tcorr3i,fcorr3i,sigmas3i]=autocorrFCSmultipletau(Itrace3(:,i));
        [tcorrcci12,fcrosscorri12,sigmascci12]=crosscorrFCSmultipletau(Itrace1(:,i),Itrace2(:,i));
        [tcorrcci13,fcrosscorri13,sigmascci13]=crosscorrFCSmultipletau(Itrace1(:,i),Itrace3(:,i));
        [tcorrcci23,fcrosscorri23,sigmascci23]=crosscorrFCSmultipletau(Itrace2(:,i),Itrace3(:,i));
        tcorri=tcorr1i*scantime; %All tcorri's equal anyways
        tcorrcci=tcorrcci12*scantime; %All tcorrcci's equal anyways
        correlationcurvesCh1(:,i)=fcorr1i';
        sigmas1curves(:,i)=sigmas1i';
        correlationcurvesCh2(:,i)=fcorr2i';
        sigmas2curves(:,i)=sigmas2i';
        correlationcurvesCh3(:,i)=fcorr3i';
        sigmas3curves(:,i)=sigmas3i';
        correlationcurvesChCC12(:,i)=fcrosscorri12';
        sigmascc12curves(:,i)=sigmascci12';
        correlationcurvesChCC13(:,i)=fcrosscorri13';
        sigmascc13curves(:,i)=sigmascci13';
        correlationcurvesChCC23(:,i)=fcrosscorri23';
        sigmascc23curves(:,i)=sigmascci23';
        
        lbi=min(tcorri);
        ubi=max(tcorri);
        p1=subplot(3,2,1);semilogx(tcorri,fcorr1i,'Color',[26/255 150/255 65/255],'LineStyle','none','Marker','.','MarkerFaceColor',[26/255 150/255 65/255]);
        title('ACF 1')
        xlabel('\tau [s]')
        ylabel('Correlation')
        hold on
        p2=subplot(3,2,2);semilogx(tcorri,fcorr2i,'Color',[253/255 174/255 97/255],'LineStyle','none','Marker','.','MarkerFaceColor',[253/255 174/255 97/255]);
        title('ACF 2')
        xlabel('\tau [s]')
        ylabel('Correlation')
        hold on
        p3=subplot(3,2,3);semilogx(tcorri,fcorr3i,'Color',[151/255 25/255 28/255],'LineStyle','none','Marker','.','MarkerFaceColor',[151/255 25/255 28/255]);
        title('ACF 3')
        xlabel('\tau [s]')
        ylabel('Correlation')
        hold on
        p4=subplot(3,2,4);semilogx(tcorrcci,fcrosscorri12,'Color', [84/255 39/255 143/255],'LineStyle','none','Marker','.','MarkerFaceColor',[84/255 39/255 143/255]);
        title('CCF 12')
        xlabel('\tau [s]')
        ylabel('Correlation')
        hold on
        p5=subplot(3,2,5);semilogx(tcorrcci,fcrosscorri13,'Color', [43/255 140/255 190/255],'LineStyle','none','Marker','.','MarkerFaceColor',[43/255 140/255 190/255]);
        title('CCF 13')
        xlabel('\tau [s]')
        ylabel('Correlation')
        hold on
        p6=subplot(3,2,6);semilogx(tcorrcci,fcrosscorri23,'Color', [150/255 150/255 150/255],'LineStyle','none','Marker','.','MarkerFaceColor',[150/255 150/255 150/255]);
        title('CCF 23')
        xlabel('\tau [s]')
        ylabel('Correlation')
        %pause
        hold on
    end
    corfit1=zeros(length(tcorr1i),ACFsegmentnumb); %Are the corfits needed? Not modified yet!
    corfit2=corfit1;
    corfit3=corfit1;
    corfitcc=corfit1;
    ItraceIICh1=Itrace1_GUI;
    ItraceIICh2=Itrace2_GUI;
    ItraceIICh3=Itrace3_GUI;
    segmenttime_GUI=scantime:scantime*binningwindow_GUI:scantime*segmentlength;
    ItraceIICh1=[segmenttime_GUI' ItraceIICh1];
    ItraceIICh2=[segmenttime_GUI' ItraceIICh2];
    ItraceIICh3=[segmenttime_GUI' ItraceIICh3];
    correlationcurvesCh1=[tcorri' correlationcurvesCh1];
    correlationcurvesCh2=[tcorri' correlationcurvesCh2];
    correlationcurvesCh3=[tcorri' correlationcurvesCh3];
    correlationcurvesChCC12=[tcorrcci' correlationcurvesChCC12];
    correlationcurvesChCC13=[tcorrcci' correlationcurvesChCC13];
    correlationcurvesChCC23=[tcorrcci' correlationcurvesChCC23];
    
    corfit1=[tcorr1i' corfit1]; %Are these used anywhere?
    corfit2=[tcorr2i' corfit2];
    corfitcc=[tcorrcci' corfitcc];
    
    curveincl1=ones(size(corfit1,2)-1,1);
    curveincl2=curveincl1;
    curveincl3=curveincl1;

    % GUI for selection of segment ACFs:
    
    %Corrselection2Channelsindividuell
    Corrselection3Channelsindividuellpreview
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
        if curveincl3(i)==0
            line3fluorescenceseries(segmentlength*(i-1)+1:i*segmentlength)=NaN;
        end
    end
    
    % Calculation of final CFs:
    [tcorr1final,fcorr1final,sigmas1final]=autocorrFCSmultipletau(line1fluorescenceseries);
    [tcorr2final,fcorr2final,sigmas2final]=autocorrFCSmultipletau(line2fluorescenceseries);
    [tcorr3final,fcorr3final,sigmas3final]=autocorrFCSmultipletau(line3fluorescenceseries);
    [tcorrcc12final,fcrosscorr12final,sigmascc12final]=crosscorrFCSmultipletau(line1fluorescenceseries,line2fluorescenceseries);
    [tcorrcc13final,fcrosscorr13final,sigmascc13final]=crosscorrFCSmultipletau(line1fluorescenceseries,line3fluorescenceseries);
    [tcorrcc23final,fcrosscorr23final,sigmascc23final]=crosscorrFCSmultipletau(line2fluorescenceseries,line3fluorescenceseries);
    tcorrfinal=tcorr1final*scantime;
    tcorrccfinal=tcorrcc12final*scantime;
    weights1final=abs(fcorr1final./sigmas1final);
    weights2final=abs(fcorr2final./sigmas2final);
    weights3final=abs(fcorr3final./sigmas3final);
    weightscc12final=abs(fcrosscorr12final./sigmascc12final);
    weightscc13final=abs(fcrosscorr13final./sigmascc13final);
    weightscc23final=abs(fcrosscorr23final./sigmascc23final);
    
    % Average of segments:
    fcorr1segfinal=zeros(length(tcorri),1);
    sigmas1segfinal=fcorr1segfinal;
    weight1seg=fcorr1segfinal;
    fcorr2segfinal=zeros(length(tcorri),1);
    sigmas2segfinal=fcorr2segfinal;
    weight2seg=fcorr2segfinal;
    fcorr3segfinal=zeros(length(tcorri),1);
    sigmas3segfinal=fcorr3segfinal;
    weight3seg=fcorr3segfinal;
    
    fcrosscorr12segfinal=zeros(length(tcorrcci),1);
    sigmascc12segfinal=fcrosscorr12segfinal;
    weightscc12seg=fcrosscorr12segfinal;
    fcrosscorr13segfinal=zeros(length(tcorrcci),1);
    sigmascc13segfinal=fcrosscorr13segfinal;
    weightscc13seg=fcrosscorr13segfinal;
    fcrosscorr23segfinal=zeros(length(tcorrcci),1);
    sigmascc23segfinal=fcrosscorr23segfinal;
    weightscc23seg=fcrosscorr23segfinal;
    for i=1:length(curveincl1);
        if curveincl1(i)==1
            fcorr1segfinal=fcorr1segfinal+correlationcurvesCh1(:,1+i);
            sigmas1segfinal=sigmas1segfinal+sigmas1curves(:,i);
        end
        if curveincl2(i)==1
            fcorr2segfinal=fcorr2segfinal+correlationcurvesCh2(:,1+i);
            sigmas2segfinal=sigmas2segfinal+sigmas2curves(:,i);
        end
        if curveincl3(i)==1
            fcorr3segfinal=fcorr3segfinal+correlationcurvesCh3(:,1+i);
            sigmas3segfinal=sigmas3segfinal+sigmas3curves(:,i);
        end
        if (curveincl1(i)+curveincl2(i))==2
        fcrosscorr12segfinal=fcrosscorr12segfinal+correlationcurvesChCC12(:,1+i);
        sigmascc12segfinal=sigmascc12segfinal+sigmascc12curves(:,i);
        end
        if (curveincl1(i)+curveincl3(i))==2
        fcrosscorr13segfinal=fcrosscorr13segfinal+correlationcurvesChCC13(:,1+i);
        sigmascc13segfinal=sigmascc13segfinal+sigmascc13curves(:,i);
        end
        if (curveincl2(i)+curveincl3(i))==2
        fcrosscorr23segfinal=fcrosscorr23segfinal+correlationcurvesChCC23(:,1+i);
        sigmascc23segfinal=sigmascc23segfinal+sigmascc23curves(:,i);
        end
    end
    fcorr1segfinal=fcorr1segfinal/sum(curveincl1);
    fcorr2segfinal=fcorr2segfinal/sum(curveincl2);
    fcorr3segfinal=fcorr3segfinal/sum(curveincl3);
    fcrosscorr12segfinal=fcrosscorr12segfinal/sum((curveincl1+curveincl2)==2);
    fcrosscorr13segfinal=fcrosscorr13segfinal/sum((curveincl1+curveincl3)==2);
    fcrosscorr23segfinal=fcrosscorr23segfinal/sum((curveincl2+curveincl3)==2);
    sigmas1segfinal=sigmas1segfinal/sum(curveincl1);
    sigmas2segfinal=sigmas2segfinal/sum(curveincl2);
    sigmas3segfinal=sigmas3segfinal/sum(curveincl3);
    sigmascc12segfinal=sigmascc12segfinal/sum((curveincl1+curveincl2)==2);
    sigmascc13segfinal=sigmascc13segfinal/sum((curveincl1+curveincl3)==2);
    sigmascc23segfinal=sigmascc23segfinal/sum((curveincl2+curveincl3)==2);
    weights1segfinal=abs(fcorr1segfinal./sigmas1segfinal);
    tcorrsegfinal=tcorri;
    weights2segfinal=abs(fcorr2segfinal./sigmas2segfinal);
    weights3segfinal=abs(fcorr3segfinal./sigmas3segfinal);
    weightscc12segfinal=abs(fcrosscorr12segfinal./sigmascc12segfinal);
    tcorrccsegfinal=tcorrcci;
    weightscc13segfinal=abs(fcrosscorr13segfinal./sigmascc13segfinal);
    weightscc23segfinal=abs(fcrosscorr23segfinal./sigmascc23segfinal);
   
    % Fit of final CFs with diffusion model:
    fprintf('Final Fitting...\n');
    lb=scantime;
    lbcc=lb;
    lbseg=lb;
    lbccseg=lb;
    ub=max(tcorrfinal);
    ubcc=ub;
    tfit=tcorrfinal;
    ubseg=max(tcorrsegfinal);
    ubccseg=ubseg;
 
    tcorrfit=tcorrfinal(1:length(tfit));
    fcorr1fit=fcorr1final(1:length(tfit));
    fcorr2fit=fcorr2final(1:length(tfit));
    fcorr3fit=fcorr3final(1:length(tfit));
    weights1fit=weights1final(1:length(tfit));
    weights1fit(weights1fit==0)=10^-4;
    weights2fit=weights2final(1:length(tfit));
    weights2fit(weights2fit==0)=10^-4;
    weights3fit=weights3final(1:length(tfit));
    weights3fit(weights3fit==0)=10^-4;
    
    tcorrccfit=tcorrccfinal(1:length(tfit));
    fcrosscorr12fit=fcrosscorr12final(1:length(tfit));
    weightscc12fit=weightscc12final(1:length(tfit));
    weightscc12fit(weightscc12fit==0)=10^-4;
    fcrosscorr13fit=fcrosscorr13final(1:length(tfit));
    weightscc13fit=weightscc13final(1:length(tfit));
    weightscc13fit(weightscc13fit==0)=10^-4;
    fcrosscorr23fit=fcrosscorr23final(1:length(tfit));
    weightscc23fit=weightscc23final(1:length(tfit));
    weightscc23fit(weightscc23fit==0)=10^-4;
    
    tcorrccsegfinalfit=tcorrccsegfinal(1:length(tcorrsegfinal));
    fcrosscorr12segfinalfit=fcrosscorr12segfinal(1:length(tcorrsegfinal));
    weightscc12segfinalfit=weightscc12segfinal(1:length(tcorrsegfinal));
    fcrosscorr13segfinalfit=fcrosscorr13segfinal(1:length(tcorrsegfinal));
    weightscc13segfinalfit=weightscc13segfinal(1:length(tcorrsegfinal));
    fcrosscorr23segfinalfit=fcrosscorr23segfinal(1:length(tcorrsegfinal));
    weightscc23segfinalfit=weightscc23segfinal(1:length(tcorrsegfinal));
    
    lsautofitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5;
    flscrossfitfunc=@(x,t)1/x(1)*((1+t./x(2)).^-0.5).*(1+t./(x(2)*x(3)^2)).^-0.5+x(4);
    
    % Fit of whole time series (try):
    try
        [N1final,taud1final,Sfit1final,CI1final,fitcurve1final,residuals1final] = autocorrfit2Ddiff(tcorrfit,fcorr1fit,lb,ub,weights1fit,lsautofitfunc,x0,fixed);
        [N2final,taud2final,Sfit2final,CI2final,fitcurve2final,residuals2final] = autocorrfit2Ddiff(tcorrfit,fcorr2fit,lb,ub,weights2fit,lsautofitfunc,x0,fixed);
        [N3final,taud3final,Sfit3final,CI3final,fitcurve3final,residuals3final] = autocorrfit2Ddiff(tcorrfit,fcorr3fit,lb,ub,weights3fit,lsautofitfunc,x0,fixed);
        [Ncc12final,taucc12final,Sfitcc12final,CIcc12final,fitcurvecc12final,residualscc12final] =autocorrfit2Ddiff(tcorrccfit,fcrosscorr12fit,lbcc,ubcc,weightscc12fit,lsautofitfunc,y0,fixed);
        [Ncc13final,taucc13final,Sfitcc13final,CIcc13final,fitcurvecc13final,residualscc13final] =autocorrfit2Ddiff(tcorrccfit,fcrosscorr13fit,lbcc,ubcc,weightscc13fit,lsautofitfunc,y0,fixed);
        [Ncc23final,taucc23final,Sfitcc23final,CIcc23final,fitcurvecc23final,residualscc23final] =autocorrfit2Ddiff(tcorrccfit,fcrosscorr23fit,lbcc,ubcc,weightscc23fit,lsautofitfunc,y0,fixed);

        % Relative cross-correlation amplitude
        relCC12=max([N1final/Ncc12final N2final/Ncc12final]);
        relCC13=max([N1final/Ncc13final N3final/Ncc13final]);
        relCC23=max([N2final/Ncc23final N3final/Ncc23final]);

        % Confidencen intervals:
        us1=0.5*(CI1final(1:3,2)-CI1final(1:3,1));
        us2=0.5*(CI2final(1:3,2)-CI2final(1:3,1));
        us3=0.5*(CI3final(1:3,2)-CI3final(1:3,1));
        uscc12=0.5*(CIcc12final(1:3,2)-CIcc12final(1:3,1));
        uscc13=0.5*(CIcc13final(1:3,2)-CIcc13final(1:3,1));
        uscc23=0.5*(CIcc23final(1:3,2)-CIcc23final(1:3,1));

        % Figure: Plot of data and fit
        h=figure('OuterPosition',[scrsz(3) 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Full Curves');
        positionvector1=[0.1 0.35 0.8 0.55];
        positionvector2=[0.1 0.1 0.8 0.15];
        subplot('Position',positionvector1),semilogx(tcorrfit,fitcurve1final,'Color',[26/255 150/255 65/255],'LineWidth',2)
        hold on 
        semilogx(tcorrfit,fitcurve2final,'Color',[253/255 174/255 97/255],'LineWidth',2)
        hold on
        semilogx(tcorrfit,fitcurve3final,'Color',[151/255 25/255 28/255],'LineWidth',2)
        hold on
        subplot('Position',positionvector1),semilogx(tcorrccfit,fitcurvecc12final,'Color', [84/255 39/255 143/255],'LineWidth',2)
        hold on
        semilogx(tcorrccfit,fitcurvecc13final,'Color', [43/255 140/255 190/255],'LineWidth',2)
        hold on
        semilogx(tcorrccfit,fitcurvecc23final,'Color', [150/255 150/255 150/255],'LineWidth',2)
        hold on
        subplot('Position',positionvector1),semilogx(tcorrfit,fcorr1fit,'s','Color',[26/255 150/255 65/255],'LineWidth',1)
        hold on
        semilogx(tcorrfit,fcorr2fit,'d','Color',[253/255 174/255 97/255],'LineWidth',1)
        hold on
        semilogx(tcorrfit,fcorr3fit,'x','Color',[151/255 25/255 28/255],'LineWidth',1)
        hold on
        subplot('Position',positionvector1),semilogx(tcorrccfit,fcrosscorr12fit,'s','Color', [84/255 39/255 143/255],'LineWidth',1)
        hold on
        semilogx(tcorrccfit,fcrosscorr13fit,'d','Color', [43/255 140/255 190/255],'LineWidth',1)
        hold on
        semilogx(tcorrccfit,fcrosscorr23fit,'x','Color', [150/255 150/255 150/255],'LineWidth',1)
        xlabel('time')
        ylabel('autocorrelation')
        subplot('Position',positionvector2),semilogx(tcorrfit,residuals1final,'s','Color',[26/255 150/255 65/255],'LineWidth',1)
        hold on
        semilogx(tcorrfit,residuals2final,'d','Color',[253/255 174/255 97/255],'LineWidth',1)
        hold on
        semilogx(tcorrfit,residuals3final,'x','Color',[151/255 25/255 28/255],'LineWidth',1)
        hold on
        subplot('Position',positionvector2),semilogx(tcorrccfit,residualscc12final,'s','Color', [84/255 39/255 143/255],'LineWidth',1)
        hold on
        semilogx(tcorrccfit,residualscc13final,'d','Color', [43/255 140/255 190/255],'LineWidth',1)
        hold on
        semilogx(tcorrccfit,residualscc23final,'x','Color', [150/255 150/255 150/255],'LineWidth',1)
        xlabel('time')
        ylabel('residuals')
        hold on
        subplot('Position',positionvector2),semilogx(tcorrfit,zeros(size(tcorrfit)),'-k','LineWidth',1.5);
        errorfull=0;
    catch
        errorfull=1;
    end
    
    % Fit of segment-averaged time series:
    [N1segfinal,taud1segfinal,Sfit1segfinal,CI1segfinal,fitcurve1segfinal,residuals1segfinal] = autocorrfit2Ddiff(tcorrsegfinal,fcorr1segfinal',lbseg,ubseg,weights1segfinal',lsautofitfunc,x0,fixed);
    [N2segfinal,taud2segfinal,Sfit2segfinal,CI2segfinal,fitcurve2segfinal,residuals2segfinal] = autocorrfit2Ddiff(tcorrsegfinal,fcorr2segfinal',lbseg,ubseg,weights2segfinal',lsautofitfunc,x0,fixed);
    [N3segfinal,taud3segfinal,Sfit3segfinal,CI3segfinal,fitcurve3segfinal,residuals3segfinal] = autocorrfit2Ddiff(tcorrsegfinal,fcorr3segfinal',lbseg,ubseg,weights3segfinal',lsautofitfunc,x0,fixed);
    [Ncc12segfinal,taucc12segfinal,Sfitcc12segfinal,CIcc12segfinal,fitcurvecc12segfinal,residualscc12segfinal] =autocorrfit2Ddiff(tcorrccsegfinalfit,fcrosscorr12segfinalfit',lbccseg,ubccseg,weightscc12segfinalfit',lsautofitfunc,y0,fixed);
    [Ncc13segfinal,taucc13segfinal,Sfitcc13segfinal,CIcc13segfinal,fitcurvecc13segfinal,residualscc13segfinal] =autocorrfit2Ddiff(tcorrccsegfinalfit,fcrosscorr13segfinalfit',lbccseg,ubccseg,weightscc13segfinalfit',lsautofitfunc,y0,fixed);
    [Ncc23segfinal,taucc23segfinal,Sfitcc23segfinal,CIcc23segfinal,fitcurvecc23segfinal,residualscc23segfinal] =autocorrfit2Ddiff(tcorrccsegfinalfit,fcrosscorr23segfinalfit',lbccseg,ubccseg,weightscc23segfinalfit',lsautofitfunc,y0,fixed);
    
    % Relative cross-correlation amplitude:
    relCC12segnew=max([N1segfinal/Ncc12segfinal N2segfinal/Ncc12segfinal]);
    relCC13segnew=max([N1segfinal/Ncc13segfinal N3segfinal/Ncc13segfinal]);
    relCC23segnew=max([N2segfinal/Ncc23segfinal N3segfinal/Ncc23segfinal]);
    
    % Confidencen intervals:
    us1seg=0.5*(CI1segfinal(1:3,2)-CI1segfinal(1:3,1));
    us2seg=0.5*(CI2segfinal(1:3,2)-CI2segfinal(1:3,1));
    us3seg=0.5*(CI3segfinal(1:3,2)-CI3segfinal(1:3,1));
    uscc12seg=0.5*(CIcc12segfinal(1:3,2)-CIcc12segfinal(1:3,1));
    uscc13seg=0.5*(CIcc13segfinal(1:3,2)-CIcc13segfinal(1:3,1));
    uscc23seg=0.5*(CIcc23segfinal(1:3,2)-CIcc23segfinal(1:3,1));
    
     % Figure: Plot of data and fit
    hh=figure('OuterPosition',[4*scrsz(3)/3 50 scrsz(3)/3 scrsz(4)/2],'Name','Final Fit Segment Averaged Curves');
    positionvector1=[0.1 0.35 0.8 0.55];
    positionvector2=[0.1 0.1 0.8 0.15];
    subplot('Position',positionvector1),semilogx(tcorrsegfinal,fitcurve1segfinal,'Color',[26/255 150/255 65/255],'LineWidth',2)
    hold on 
    semilogx(tcorrsegfinal,fitcurve2segfinal,'Color',[253/255 174/255 97/255],'LineWidth',2)
    hold on
    semilogx(tcorrsegfinal,fitcurve3segfinal,'Color',[151/255 25/255 28/255],'LineWidth',2)
    hold on
    subplot('Position',positionvector1),semilogx(tcorrccsegfinalfit,fitcurvecc12segfinal,'Color', [84/255 39/255 143/255],'LineWidth',2)
    hold on
    semilogx(tcorrccsegfinalfit,fitcurvecc13segfinal,'Color', [43/255 140/255 190/255],'LineWidth',2)
    hold on
    semilogx(tcorrccsegfinalfit,fitcurvecc23segfinal,'Color', [150/255 150/255 150/255],'LineWidth',2)
    hold on
    subplot('Position',positionvector1),semilogx(tcorrsegfinal,fcorr1segfinal,'s','Color',[26/255 150/255 65/255],'LineWidth',1)
    hold on
    semilogx(tcorrsegfinal,fcorr2segfinal,'d','Color',[253/255 174/255 97/255],'LineWidth',1)
    hold on
    semilogx(tcorrsegfinal,fcorr3segfinal,'x','Color',[151/255 25/255 28/255],'LineWidth',1)
    hold on
    subplot('Position',positionvector1),semilogx(tcorrccsegfinalfit,fcrosscorr12segfinalfit,'s','Color', [84/255 39/255 143/255],'LineWidth',1)
    hold on
    semilogx(tcorrccsegfinalfit,fcrosscorr13segfinalfit,'d','Color', [43/255 140/255 190/255],'LineWidth',1)
    hold on
    semilogx(tcorrccsegfinalfit,fcrosscorr23segfinalfit,'x','Color', [150/255 150/255 150/255],'LineWidth',1)
    xlabel('time')
    ylabel('autocorrelation')
    subplot('Position',positionvector2),semilogx(tcorrsegfinal,residuals1segfinal,'s','Color',[26/255 150/255 65/255],'LineWidth',1)
    hold on
    semilogx(tcorrsegfinal,residuals2segfinal,'d','Color',[253/255 174/255 97/255],'LineWidth',1)
    hold on
    semilogx(tcorrsegfinal,residuals3segfinal,'x','Color',[151/255 25/255 28/255],'LineWidth',1)
    hold on
    subplot('Position',positionvector2),semilogx(tcorrccsegfinalfit,residualscc12segfinal,'s','Color', [84/255 39/255 143/255],'LineWidth',1)
    hold on
    semilogx(tcorrccsegfinalfit,residualscc13segfinal,'d','Color', [43/255 140/255 190/255],'LineWidth',1)
    hold on
    semilogx(tcorrccsegfinalfit,residualscc23segfinal,'x','Color', [150/255 150/255 150/255],'LineWidth',1)
    xlabel('time')
    ylabel('residuals')
    hold on
    subplot('Position',positionvector2),semilogx(tcorrsegfinal,zeros(size(tcorrsegfinal)),'-k','LineWidth',1.5);
    
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
        output=zeros(length(tcorr2),14);    
        output(:,1)=tcorr1final;
        output(:,2)=fcorr1final;
        output(:,3)=weights1final;
        output(:,4)=fcorr2final;
        output(:,5)=weights2final;
        output(:,6)=fcorr3final;
        output(:,7)=weights3final;
        output(:,8)=tcorrccfinal;
        output(:,9)=fcrosscorr12final;
        output(:,10)=weightscc12final;
        output(:,11)=fcrosscorr13final;
        output(:,12)=weightscc13final;
        output(:,13)=fcrosscorr23final;
        output(:,14)=weightscc23final;

        outputparameters=zeros(24,3);
        outputparameters(1:3,1)=[N1final;taud1final;Sfit1final];
        outputparameters(1:3,2:end)=CI1final;
        outputparameters(4:6,1)=[N2final;taud2final;Sfit2final];
        outputparameters(4:6,2:end)=CI2final;
        outputparameters(7:9,1)=[N3final;taud3final;Sfit3final];
        outputparameters(7:9,2:end)=CI3final;
        outputparameters(10:12,1)=[Ncc12final;taucc12final;Sfitcc12final];
        outputparameters(10:12,2:end)=CIcc12final(1:3,:);
        outputparameters(13:15,1)=[Ncc13final;taucc13final;Sfitcc13final];
        outputparameters(13:15,2:end)=CIcc13final(1:3,:);
        outputparameters(16:18,1)=[Ncc23final;taucc23final;Sfitcc23final];
        outputparameters(16:18,2:end)=CIcc23final(1:3,:);
        outputparameters(19,1)=bleachingfraction1;
        outputparameters(20,1)=bleachingfraction2;
        outputparameters(21,1)=bleachingfraction3;
        outputparameters(22,1)=nanmean(line1fluorescenceseries);
        outputparameters(23,1)=nanmean(line2fluorescenceseries);
        outputparameters(24,1)=nanmean(line3fluorescenceseries);
        
        fprintf(fid1,'%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output');
        fprintf(fid2,'%e\t %e\t %e\n',outputparameters');
    end
    
    % Segment analysis:
    outputseg=zeros(length(tcorrsegfinal),14);    
    outputseg(:,1)=tcorrsegfinal;
    outputseg(:,2)=fcorr1segfinal;
    outputseg(:,3)=weights1segfinal;
    outputseg(:,4)=fcorr2segfinal;
    outputseg(:,5)=weights2segfinal;
    outputseg(:,6)=fcorr3segfinal;
    outputseg(:,7)=weights3segfinal;
    outputseg(:,8)=tcorrccsegfinal;
    outputseg(:,9)=fcrosscorr12segfinal;
    outputseg(:,10)=weightscc12segfinal;
    outputseg(:,11)=fcrosscorr13segfinal;
    outputseg(:,12)=weightscc13segfinal;
    outputseg(:,13)=fcrosscorr23segfinal;
    outputseg(:,14)=weightscc23segfinal;
    
    outputparametersseg=zeros(24,3);
    outputparametersseg(1:3,1)=[N1segfinal;taud1segfinal;Sfit1segfinal];
    outputparametersseg(1:3,2:end)=CI1segfinal;
    outputparametersseg(4:6,1)=[N2segfinal;taud2segfinal;Sfit2segfinal];
    outputparametersseg(4:6,2:end)=CI2segfinal;
    outputparametersseg(7:9,1)=[N3segfinal;taud3segfinal;Sfit3segfinal];
    outputparametersseg(7:9,2:end)=CI3segfinal;
    outputparametersseg(10:12,1)=[Ncc12segfinal;taucc12segfinal;Sfitcc12segfinal];
    outputparametersseg(10:12,2:end)=CIcc12segfinal(1:3,:);
    outputparametersseg(13:15,1)=[Ncc13segfinal;taucc13segfinal;Sfitcc13segfinal];
    outputparametersseg(13:15,2:end)=CIcc13segfinal(1:3,:);
    outputparametersseg(16:18,1)=[Ncc23segfinal;taucc23segfinal;Sfitcc23segfinal];
    outputparametersseg(16:18,2:end)=CIcc23segfinal(1:3,:);
    
    outputparametersseg(19,1)=bleachingfraction1;
    outputparametersseg(20,1)=bleachingfraction2;
    outputparametersseg(21,1)=bleachingfraction3;
    outputparametersseg(22,1)=nanmean(line1fluorescenceseries);
    outputparametersseg(23,1)=nanmean(line2fluorescenceseries);
    outputparametersseg(24,1)=nanmean(line3fluorescenceseries);
    
    fprintf(fid3,'%e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',outputseg');
    fprintf(fid4,'%e\t %e\t %e\n',outputparametersseg');
    
    curveincls=[curveincl1 curveincl2 curveincl3];
    fprintf(fid6,'%e\t %e\t %e\n',curveincls');
    
    % Save figures:
    saveas(hh,[path2 '\' inputfilename(1:end-4) ' CCFs_seg.fig'])
    saveas(hI,[path2 '\' inputfilename(1:end-4) ' I_decomposed.fig'])
end
