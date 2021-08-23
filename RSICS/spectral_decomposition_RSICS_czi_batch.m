clear all
close all

% Parameters for decomposition
channels=3; %Number of differently tagged species
timeseries=1; %Time series acquisition (if yes, set to 1) or single image (set to 0)?


% Load directory containing spectral images (.czi):
path= uigetdir; 
% Load directory where decomposed image stacks shall be saved:
path2=uigetdir;
% Load .txt file containing spectral patterns (for the specified number of
% channels):
[inputfilenamespecpatterns, path3]=uigetfile('*.txt');


% Loading files and performing spectral decomposition:
files=dir([path '/*.czi']);
[namedata,remain]=strtok(files(1).name,'.');
for i=1:size(files,1)
    inputfilename=files(i).name;
        
    inputimage=bfopen_b(path, inputfilename);
    
    close all % Close all graphs from previous file
    metadata=inputimage{1,2};


    lastchannel=str2double(metadata.get('Global Information|Image|SizeC #1'));
    startchannel=1;
    
    if timeseries==0
       numberframes=1;
       numberofchannels=23; % Specify here the number of channels for single image acquisition    
    else
       numberframes=str2double(metadata.get('Global Information|Image|SizeT #1'));
       numberofchannels=size(inputimage{1,1},1)/numberframes;
    end
    imagestack4D=zeros(size(inputimage{1,1}{1,1},1),size(inputimage{1,1}{1,1},2),numberframes,numberofchannels);
    channeli=zeros(size(inputimage{1,1}{1,1},1),size(inputimage{1,1}{1,1},2),numberframes);
    for ii=startchannel:lastchannel
        for iii=1:size(inputimage{1,1},1)
            number=(iii-ii+numberofchannels)/numberofchannels;
            if number==round(number)
                channeli(:,:,number)=inputimage{1,1}{iii,1};
            end
        end
        imagestack4D(:,:,:,ii)=channeli;
    end

    spectralbins=zeros(1,numberofchannels);
    spectralbins(1)=str2double(metadata.get('Global Information|Image|Channel|EmissionWavelength #01')); % first wavelength from metadata
    for iii=2:numberofchannels % read out all wavelengths from metadata
        if iii<10
           spectralbins(iii)=str2double(metadata.get(['Global Information|Image|Channel|EmissionWavelength #0'  num2str(iii)   ] ));
        end
        if iii>=10
           spectralbins(iii)=str2double(metadata.get(['Global Information|Image|Channel|EmissionWavelength #'  num2str(iii)   ] ));
        end
    end

    % Figure of time- and channel-averaged image:
    meanframe=squeeze(mean(mean(imagestack4D,4),3));
    figure('Name','Average image')
    imagesc(meanframe)

    % Calculation of pixel- and time-averaged spectrum:
    spectralcounts=size(spectralbins);
    for iii=1:numberofchannels
        spectralcounts(iii)=nansum(nansum(nansum(imagestack4D(:,:,:,iii))));
    end
    spectralweights=spectralcounts/sum(spectralcounts);

    % Plot of average spectrum
    emspec=figure('Name','Spectrum');
    plot(spectralbins,spectralweights)
    xlabel('Channel [nm]')
    ylabel('Norm.emission')
    
    % Calculation of spectral filters w (according to Schrimpf et al.,
    % Methods 2018):
    spectral_patterns=load([path3 '/' inputfilenamespecpatterns]);
    spectral_patterns=spectral_patterns';
    imagestack_GFP=zeros(size(imagestack4D,1),size(imagestack4D,2),numberframes);
    imagestack_YFP=imagestack_GFP;
    imagestack_Ch=imagestack_GFP;
    if channels==4
        imagestack_A=imagestack_GFP;
    end
    imagestack_all=zeros(size(imagestack4D,1),size(imagestack4D,2),numberframes,channels);
    Ittt=squeeze(nanmean(nanmean(nanmean(imagestack4D,3),2),1));
    Iinv(Ittt>0)=1./Ittt(Ittt>0);
    D=diag(Iinv);
    w=(spectral_patterns*D*spectral_patterns')\spectral_patterns*D;

    
    % Figure spectral weights:
    specfilt=figure('Name','Spectral filters');
    if channels==3
        plot(spectralbins,w(1,:),'-g',spectralbins,w(2,:),'-y',spectralbins,w(3,:),'-r');
    end
    if channels==4
        plot(spectralbins,w(1,:),'-g',spectralbins,w(2,:),'-y',spectralbins,w(3,:),'-m',spectralbins,w(4,:),'-r');
    end
    savefig(specfilt,[path2 '\' inputfilename(1:end-4) '_specfilters.fig'])

    % Spectral decomposition of each frame:
    for ttt=1:size(imagestack4D,3)
                imageframe=zeros(size(imagestack4D,1),size(imagestack4D,2));
                for iii=1:numberofchannels
                    imageframe_channel=w(1,iii).*squeeze(imagestack4D(:,:,ttt,iii));
                    imageframe=imageframe+imageframe_channel;
                end
                imagestack_GFP(:,:,ttt)=imageframe;
                imageframe=zeros(size(imagestack4D,1),size(imagestack4D,2));
                for iii=1:numberofchannels
                    imageframe_channel=w(2,iii).*squeeze(imagestack4D(:,:,ttt,iii));
                    imageframe=imageframe+imageframe_channel;
                end
                imagestack_YFP(:,:,ttt)=imageframe;
                
                if channels==3
                    imageframe=zeros(size(imagestack4D,1),size(imagestack4D,2));
                    for iii=1:numberofchannels
                        imageframe_channel=w(3,iii).*squeeze(imagestack4D(:,:,ttt,iii));
                        imageframe=imageframe+imageframe_channel;
                    end
                    imagestack_Ch(:,:,ttt)=imageframe;
                end
                if channels==4
                    imageframe=zeros(size(imagestack4D,1),size(imagestack4D,2));
                    for iii=1:numberofchannels
                        imageframe_channel=w(3,iii).*squeeze(imagestack4D(:,:,ttt,iii));
                        imageframe=imageframe+imageframe_channel;
                    end
                    imagestack_A(:,:,ttt)=imageframe;
                    imageframe=zeros(size(imagestack4D,1),size(imagestack4D,2));
                    for iii=1:numberofchannels
                        imageframe_channel=w(4,iii).*squeeze(imagestack4D(:,:,ttt,iii));
                        imageframe=imageframe+imageframe_channel;
                    end
                    imagestack_Ch(:,:,ttt)=imageframe;
                end   
    end

    % Save decomposed image stacks: 
    if channels==3
        meanframe_GFP=squeeze(mean(imagestack_GFP,3));
        meanframe_YFP=squeeze(mean(imagestack_YFP,3));
        meanframe_Ch=squeeze(mean(imagestack_Ch,3));

        imagestack_all(:,:,:,1)=imagestack_GFP;
        imagestack_all(:,:,:,2)=imagestack_YFP;
        imagestack_all(:,:,:,3)=imagestack_Ch;


        save([path2 '\' inputfilename(1:end-4) '_imagestack_GFP.mat'],'imagestack_GFP')
        save([path2 '\' inputfilename(1:end-4) '_imagestack_YFP.mat'],'imagestack_YFP')
        save([path2 '\' inputfilename(1:end-4) '_imagestack_Ch2.mat'],'imagestack_Ch')
    end
    if channels==4
        meanframe_GFP=squeeze(mean(imagestack_GFP,3));
        meanframe_YFP=squeeze(mean(imagestack_YFP,3));
        meanframe_A=squeeze(mean(imagestack_A,3));
        meanframe_Ch=squeeze(mean(imagestack_Ch,3));

        imagestack_all(:,:,:,1)=imagestack_GFP;
        imagestack_all(:,:,:,2)=imagestack_YFP;
        imagestack_all(:,:,:,3)=imagestack_A;
        imagestack_all(:,:,:,4)=imagestack_Ch;


        save([path2 '\' inputfilename(1:end-4) '_imagestack_GFP.mat'],'imagestack_GFP')
        save([path2 '\' inputfilename(1:end-4) '_imagestack_YFP.mat'],'imagestack_YFP')
        save([path2 '\' inputfilename(1:end-4) '_imagestack_A.mat'],'imagestack_A')
        save([path2 '\' inputfilename(1:end-4) '_imagestack_Ch2.mat'],'imagestack_Ch')
    end

    % Plot decomposed average images:
    load('gfpmap2.mat')
    load('yfpmap2.mat')
    load('mchmap2.mat')
    figstacks=figure('Name','Decomposed Image stacks');
    h1=subplot(2,2,1);imagesc(meanframe_GFP);
    colormap(h1,gfpmap2);
    title('Average Image G')
    h2=subplot(2,2,2);imagesc(meanframe_YFP)
    colormap(h2,yfpmap2);
    title('Average Image Y')
    if channels==3
        h3=subplot(2,2,3);imagesc(meanframe_Ch);
        colormap(h3,mchmap2);
        title('Average Image Ch2')
        h4=subplot(2,2,4);imagesc(meanframe_GFP+meanframe_YFP+meanframe_Ch);
        colormap(h4,gray)
        title('Average Image')
    end
    if channels==4
        subplot(2,2,3),imagesc(meanframe_A)
        title('Average Image A')
        subplot(2,2,4),imagesc(meanframe_Ch)
        title('Average Image Ch2')
    end
    saveas(figstacks,[path2 '\' inputfilename(1:end-4) '_decomposed_meanimages.fig']);

    % Plot decomposed fluorescence series (pixel averaged image intensity):
    figure('Name','Decomposed fluorescence series')
    plot(squeeze(mean(mean(imagestack_GFP))),'-g')
    hold on
    plot(squeeze(mean(mean(imagestack_YFP))),'-y')
    hold on
    if channels==3
        plot(squeeze(mean(mean(imagestack_Ch))),'-r')
        hold on
    end
    if channels==4
        plot(squeeze(mean(mean(imagestack_A))),'-m')
        hold on
        plot(squeeze(mean(mean(imagestack_Ch))),'-r')
        hold on
    end
    plot(squeeze(mean(mean(sum(imagestack_all,4)))),'-b')

    
    % Calculation of spectral fractions (by fit of linear combination of spectral patterns to average spectrum):
    if channels==3
        
        spectraldecomposfitfunc=@(x,t) x(1)*spectral_patterns(1,:)+x(2)*spectral_patterns(2,:)+x(3)*spectral_patterns(3,:);%+x(4)*avgspectrum_Card;
        x0=[0.5 0.5 0.5];
        fixed=[false false false];
        [fitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixed,spectralbins,spectralweights,spectraldecomposfitfunc,x0);
        fractionGFP=fitparameters(1)
        fractionYFP=fitparameters(2)
        fractionCh=fitparameters(3)
        %fractionCard=fitparameters(4)

        fractions=[fractionGFP fractionYFP fractionCh];

        fidfrac=fopen([path2  '\specfractions_GFP_YFP_Ch2.txt'],'a');
        %fprintf(fidfrac,'sample\t fracGFP\t fracYFP\t fracmCh\n');
        fprintf(fidfrac,inputfilename(1:end-4));
        fprintf(fidfrac,'\t %e\t %e\t %e\n',fractions');
    else
        if channels==4
            spectraldecomposfitfunc=@(x,t) x(1)*spectral_patterns(1,:)+x(2)*spectral_patterns(2,:)+x(3)*spectral_patterns(3,:)+x(4)*spectral_patterns(4,:);%+x(4)*avgspectrum_Card;
            x0=[0.5 0.5 0.5 0.5];
            fixed=[false false false false];
            [fitparameters,residuals,J,COVB,MSE] = nlinfitsome(fixed,spectralbins,spectralweights,spectraldecomposfitfunc,x0);
            fractionGFP=fitparameters(1)
            fractionYFP=fitparameters(2)
            fractionA=fitparameters(3)
            fractionCh=fitparameters(4)

            fractions=[fractionGFP fractionYFP fractionA fractionCh];

            fidfrac=fopen([path2  '\specfractions_GFP_YFP_A_Ch2.txt'],'a'); 
            %fprintf(fidfrac,'sample\t fracGFP\t fracYFP\t fracmCh\n');
            fprintf(fidfrac,inputfilename(1:end-4));
            fprintf(fidfrac,'\t %e\t %e\t %e\t %e\n',fractions');
        end
        
    end
    
    % Plot of average spectrum and fit
    emspec=figure('Name','Fluorescence Spectra');
    avgspectrum_all_fit=spectraldecomposfitfunc(fitparameters,spectralbins);
    plot(spectralbins,spectralweights);
    hold on
    plot(spectralbins,avgspectrum_all_fit,'--');
    xlabel('Channel [nm]')
    ylabel('Norm.emission')
    
    savefig(emspec,[path2 '\' inputfilename(1:end-4) '_Flpectrum.fig'])
    specoutput=zeros(numberofchannels,2);
    specoutput(:,1)=spectralbins';
    specoutput(:,2)=spectralweights';
    
    fidspec=fopen([path2 '\' inputfilename(1:end-4) '_flspectrum.txt'],'a'); % adjust path if necessary!
    fprintf(fidspec,'%e\t %e\n',specoutput');
end
