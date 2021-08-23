
clear 
close all
[filename, pathname]=uigetfile ( 'Z:\*.czi');
mapfirst=0; %not needed in this version
pathbatch=pathname;

shotnoisecorrection=1;   %how many maximum pixels have to be removed from the center of the correlation surface
filter=1; % 1=turning moving average filter on. Otherwise set to 0 (e.g. for slow membrane dynamics)
filterwindow=4; % frames for moving average filter
TICS_anal=0; %not needed in this version

%%% not needed in this version
automaticmapping=0; % 1= making maps of diffusion coefficents
mappingfirstframe=mapfirst; %0 means calculating D from average of all frames; 1 means calculating only from 1st frame
ssize=64; % size in pixels of the mapping subunit
mappingstep=8; %mapping density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



startframe=1; %initial frame


                                                                                                 
                                                                                                                 
                                                                                                                 
                                                                                                                 

                                                                                                         
                                                                                                                 
data=bfopen ([pathname  filename]);


%%%%% this section (lines 35-79) reads imaging settings from the czi file. It must be
%%%%% adapted to the specific microscopy setup used
imagemetadata=data{1,2};
allentries = arrayfun(@char, imagemetadata.keySet.toArray, 'UniformOutput', false);
for i=1:size(allentries,1)
allentries{i,2}=imagemetadata.get(allentries{i,1});
end                                                                                                                 
        
pixelnumb=str2double(imagemetadata.get('Global Information|Image|SizeX #1'));

SizeX=str2double(imagemetadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'))*pixelnumb*1e6;  %in µm
pixeltime=str2double(imagemetadata.get('Global Information|Image|Channel|LaserScanInfo|PixelTime #01'));
numberframes=str2double(imagemetadata.get('Global Information|Image|SizeT #1'));
pixelsize=str2double(imagemetadata.get('Global Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'))*1e+6;

if isnan(numberframes)
   numberframes=1;
end

 %%%% We noticed a problem with the way linetime and frametime are
 %%%% (sometime) written in czi files so what follows is an approximation
 %%%% that can be used in those cases
  
     if pixeltime>8.1e-06 && pixeltime<8.27e-06
    linetime=0.00492;
    frametime=1.26; %or more! check!
     else
         if pixeltime>4.0e-06 && pixeltime<4.3e-06
     linetime=0.00246;
     frametime=0.62914; %or more! check!
         else
             if pixeltime>2.0e-06 && pixeltime<2.1e-06
     linetime=0.00123;
    frametime=0.31457; %or more! check!
             else
                 if pixeltime>1.0e-06 && pixeltime<1.1e-06
     linetime=0.0006144;
    frametime=0.15729; %or more! check!
    else 
         linetime=pixeltime*2*pixelnumb;  %APPROXIMATED
        frametime=linetime*pixelnumb;
                 end
             end
          end
      end
    
endframe=numberframes;                                                                                                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                                                 
 
%The three channels are referred to (also in the other scripts) as GFP, YFP
%and Ch2 channels

 Ch2filename=[filename(1:end-4) '_imagestack_Ch2.mat'];
 GFPfilename=[filename(1:end-4) '_imagestack_GFP.mat']; 
 YFPfilename=[filename(1:end-4) '_imagestack_YFP.mat'];
 
 
data=load ([pathname  Ch2filename]);
nameoffield=fieldnames(data);
namefield=nameoffield{1,1};  
 

preframe=data.(namefield);
       
% checking that the parameters for filtering are OK
if endframe-startframe+1>numberframes || isnan(numberframes)
    fprintf(1,'There are not enough frames\n')
    
    
    if isnan(numberframes)
        fprintf(1,'There is only one frame\n')
    else
    
    end
    error('')
    
    
end

frame=double(preframe(:,:,startframe:endframe));
timetrace=squeeze(mean(mean(frame)));

% plotting the time trace
if automaticmapping==0
figure(50)
plot(timetrace)
pause
close 50
end

 figure('units','normalized','outerposition',[0 0 1 1])
 
% plotting the sum of all frames
imagesc(sum(frame,3))
colormap gray
k=1;

%starting the selection of image portions for later analysis
while true 
Use=roipoly;
Usenan=double(Use);
Usenan(Usenan==0)=NaN;
figure(50)
subplot(1,2,1)
try
imagesc(mean(frame,3).*Usenan)
catch 
    close all
    error('Done'); %ends the while cycle if the image was closed
end
subplot(1,2,2)
plot(squeeze(nanmean(nanmean(frame.*Usenan))))

button = questdlg( 'Do you want to keep this selection?', ' ','Yes','No','Yes');
switch button
    case 'Yes'
filename=Ch2filename;
note='';
prompt={'Add note'}; %any note from user that will be written later in the result file (e.g. "apoptotic cell")
dims=[1 35];
definput={note};
note=inputdlg(prompt,' ', dims, definput);
save([pathbatch '\ARICSbatch_' filename '_' num2str(k) '.mat'], '-regexp', '^(?!(data|frame|imagemetadata|simulationresult)$).');
filename=GFPfilename;
save([pathbatch '\ARICSbatch_' filename '_' num2str(k) '.mat'], '-regexp', '^(?!(data|frame|frame1|imagemetadata|simulationresult)$).'); 
filename=YFPfilename;
save([pathbatch '\ARICSbatch_' filename '_' num2str(k) '.mat'], '-regexp', '^(?!(data|frame|frame1|imagemetadata|simulationresult)$).');
k=k+1;
clear framea Use Usenan
close 50

    case 'No'
     clear framea Use Usenan
     close 50   
end
end
 