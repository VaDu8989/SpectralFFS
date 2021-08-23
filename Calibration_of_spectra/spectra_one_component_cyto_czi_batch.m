% Calculate single-species spectra from imagestack of membrane-anchored
% fluorescent proteins 

% Load directory of image files:
path_images= uigetdir; 
 
% Load path where spectra are saved:
path2=uigetdir;

% Loading files:
files=dir([path_images '/*.czi']); 
for i=1:size(files,1) % go through all acquired measurements/images
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    inputbackground=bfopen_b(path_images, namefile); 
    metadata=inputbackground{1,2}; % read out metadata

    lastchannel=str2double(metadata.get('Global Information|Image|SizeC #1'));
    startchannel=1;
    numberframes=str2double(metadata.get('Global Information|Image|SizeT #1'));
    numberofchannels=size(inputbackground{1,1},1)/numberframes;

    imagestack4D=zeros(size(inputbackground{1,1}{1,1},1),size(inputbackground{1,1}{1,1},2),numberframes,numberofchannels);
    channeli=zeros(size(inputbackground{1,1}{1,1},1),size(inputbackground{1,1}{1,1},2),numberframes);
    for ii=startchannel:lastchannel
        for iii=1:size(inputbackground{1,1},1)
            number=(iii-ii+numberofchannels)/numberofchannels;
            if number==round(number)
                channeli(:,:,number)=inputbackground{1,1}{iii,1};
            end
        end
        imagestack4D(:,:,:,ii)=channeli;
    end

    spectralbins=zeros(1,numberofchannels);
    spectralbins(1)=str2double(metadata.get('Global Information|Image|Channel|EmissionWavelength #01')); % wavelength of first spectral bin
    for iii=2:numberofchannels % Read out all wavelengths
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
    
    % Draw region-of-interest:
    mask=roipoly;
    meanframeroifinal=mask.*meanframe;
    
    % Plot selected region:
    figure('Name','Selection')
    imagesc(meanframeroifinal)

    % Applying mask to 4D stack:
    imagestack4Dmasked=imagestack4D;
    for ii=1:size(imagestack4D,3)
        for j=1:size(imagestack4D,4)
            imagestack4Dmasked(:,:,ii,j)=squeeze(imagestack4D(:,:,ii,j)).*mask;
        end
    end

    % Calculate spectrum:
    spectralcounts=size(spectralbins);
    for iii=1:numberofchannels
        spectralcounts(iii)=sum(sum(sum(imagestack4Dmasked(:,:,:,iii))));
    end
    spectralweights=spectralcounts/sum(spectralcounts);

    % Plot spectrum:
    emspec=figure('Name','Spectrum');
    plot(spectralbins,spectralweights)
    xlabel('Channel [nm]')
    ylabel('Norm.emission')

    % Save spectra (one for each file):
    savefig(emspec,[path2 '\' namefile(1:end-4) '_Flpectrum.fig'])
    specoutput=zeros(numberofchannels,2);
    specoutput(:,1)=spectralbins';
    specoutput(:,2)=spectralweights';
    fidspec=fopen([path2 '\' namefile(1:end-4) '_flspectrum.txt'],'a');
    fprintf(fidspec,'%e\t %e\n',specoutput');
    close all
    clear mask
end


