% Pool extracted single species spectra and generate spectral patterns
% (.txt file is generated as ouput)
close all

specchannelnumb=23; % specify here the number of spectral bins
pathGFP= uigetdir;  % select directory in which extracted spectra from species 1 were saved
pathYFP= uigetdir; % select directory in which extracted spectra from species 2 were saved
pathCh=uigetdir; % select directory in which extracted spectra from species 3 were saved
path2=uigetdir; % select directory in which spectral_patterns.txt file should be saved

filesGFP=dir([pathGFP '/*flspectrum.txt']);
filesYFP=dir([pathYFP '/*flspectrum.txt']);
filesCh=dir([pathCh '/*flspectrum.txt']);

spectralchannels=zeros(1,specchannelnumb);
spectra_GFP=zeros(specchannelnumb,size(filesGFP,1));
spectra_YFP=zeros(specchannelnumb,size(filesYFP,1));
spectra_Ch=zeros(specchannelnumb,size(filesCh,1));

% load all spectra from species 1:
for i=1:size(filesGFP,1)
    [namedata,remain]=strtok(filesGFP(i).name,'.');
    namefile=filesGFP(i).name;
    specdataGFP=load([pathGFP '/' namefile]);
    spectralchannels=specdataGFP(:,1);
    spectra_GFP(:,i)=specdataGFP(:,2);
end
% load all spectra from species 2:
for i=1:size(filesYFP,1)
    [namedata,remain]=strtok(filesYFP(i).name,'.');
    namefile=filesYFP(i).name;
    specdataYFP=load([pathYFP '/' namefile]);
    spectra_YFP(:,i)=specdataYFP(:,2);
end
% load all spectra from species 3:
for i=1:size(filesCh,1)
    [namedata,remain]=strtok(filesCh(i).name,'.');
    namefile=filesCh(i).name;
    specdataCh=load([pathCh '/' namefile]);
    spectra_Ch(:,i)=specdataCh(:,2);
end

avgspectrum_GFP=mean(spectra_GFP,2);
avgspectrum_YFP=mean(spectra_YFP,2);
avgspectrum_Ch=mean(spectra_Ch,2);

stdspectraGFP=std(spectra_GFP,1,2);
stdspectraYFP=std(spectra_YFP,1,2);
stdspectraCh=std(spectra_Ch,1,2);

% Plot of average spectra and their variation:
figure('Name','Fluorescent Spectra')
errorbar(spectralchannels,avgspectrum_GFP,stdspectraGFP/2,stdspectraGFP/2);
hold on
errorbar(spectralchannels,avgspectrum_YFP,stdspectraYFP/2,stdspectraYFP/2);
hold on
errorbar(spectralchannels,avgspectrum_Ch,stdspectraCh/2,stdspectraCh/2);
xlabel('Channel [nm]')
ylabel('Norm.emission')
legend('Spectrum mEGFP','Spectrum mEYFP','Spectrum mCherry')

% Generate output .txt file:
p_ij=[avgspectrum_GFP';avgspectrum_YFP';avgspectrum_Ch'];
fidpij=fopen([path2 '\spectral_patterns_G_Y_Ch.txt'],'a'); % adjust filename if necessary
fprintf(fidpij,'%e\t %e\t %e\n',p_ij);

    
