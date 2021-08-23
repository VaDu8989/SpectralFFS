clear  
close all
batchesdir=uigetdir;
listofbatches=dir([batchesdir '/*RICSbatch*']);
kk=1;
listofbatches2=zeros(size(listofbatches,1));
for k=1:size(listofbatches,1)
if listofbatches(k,1).name(end-12:end-10)=='YFP' %#ok<BDSCA>
listofbatches2(kk)=listofbatches(k);
kk=kk+1;
end
end

%cycles through the different batch files that are planned to be analyzed
for kkk=1:size(listofbatches2,2)
fid = fopen([batchesdir '\CCanalisi ARICS new fitmodel.txt'],'a');%where it writes results
pathpdf = [ batchesdir '\CCanalisi ARICS new fitmodel.pdf' ];
load([batchesdir '\' listofbatches2(1,kkk).name]);
lowlim= [0.19 0 0 -0.0001 3 0];  %%lower limit for fit:  w0   N   D   offset S
up = [0.3 100000 100 0.0001 7 13];  %upper limit
a1 = [0.25 500 50 0 5 11]; %starting parameters
                                                                                                 
data=load ([pathname  Ch2filename]);
nameoffield=fieldnames(data);
namefield=nameoffield{1,1};  
 
Ch2acq=data.(namefield);
    
                                                                                                         
data=load ([pathname  GFPfilename]);
nameoffield=fieldnames(data);
namefield=nameoffield{1,1};  
 
GFPacq=data.(namefield);
    
                                                                                                         
data=load ([pathname  YFPfilename]);
nameoffield=fieldnames(data);
namefield=nameoffield{1,1};  
 
YFPacq=data.(namefield);
    
data=load ([pathname  Afilename]);
nameoffield=fieldnames(data);
namefield=nameoffield{1,1};  
 
Aacq=data.(namefield);
    
    
Ch2frame=double(Ch2acq(:,:,startframe:endframe));
GFPframe=double(GFPacq(:,:,startframe:endframe));
YFPframe=double(YFPacq(:,:,startframe:endframe));
Aframe=double(Aacq(:,:,startframe:endframe));
clear frame



%%%%%%%%%%%%%%%%%%filtering slow dynamics, moving average (until line 90)
if endframe-startframe+1<filterwindow+1 && filter==1
   error('Error! You are trying to apply a running average with a window larger than the total amount of frames! \n Remove the filterning by setting filter=0 \n') 
end

Ch2nonfilteredframe=Ch2frame;
GFPnonfilteredframe=GFPframe;
YFPnonfilteredframe=YFPframe;
Anonfilteredframe=Aframe;

Ch2totalmean=mean(mean(mean(Ch2frame(:,:,:))));
GFPtotalmean=mean(mean(mean(GFPframe(:,:,:))));
YFPtotalmean=mean(mean(mean(YFPframe(:,:,:)))); 
Atotalmean=mean(mean(mean(Aframe(:,:,:))));

if filter==1 && endframe-startframe+1>filterwindow+1
    
Ch2frame2=RICSfilter(Ch2frame, filterwindow, startframe, endframe);
GFPframe2=RICSfilter(GFPframe, filterwindow, startframe, endframe);
YFPframe2=RICSfilter(YFPframe, filterwindow, startframe, endframe);
Aframe2=RICSfilter  (Aframe, filterwindow, startframe, endframe);
      
clear frame Ch2frame GFPframe YFPframe Aframe
Ch2frame=Ch2frame2+Ch2totalmean;
GFPframe=GFPframe2+GFPtotalmean;
YFPframe=YFPframe2+YFPtotalmean;
Aframe=Aframe2+Atotalmean;
clear frame2 Ch2frame2 GFPframe2 YFPframe2 Aframe2
end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Calculation of Cross-correlation between different channels, based on
%ARICS algorithm (DOI: 10.1016/j.bpj.2016.09.012)

ICS2DCCorr_Ch2_GFP = arbitrary_V0CC(Ch2frame, GFPframe,Use);
ICS2DCCorr_Ch2_YFP = arbitrary_V0CC(Ch2frame, YFPframe, Use);
ICS2DCCorr_GFP_YFP = arbitrary_V0CC(GFPframe, YFPframe, Use);
ICS2DCCorr_GFP_A = arbitrary_V0CC(GFPframe, Aframe, Use);  %A is the name of the fourth channel
ICS2DCCorr_A_YFP = arbitrary_V0CC(Aframe, YFPframe, Use);
ICS2DCCorr_A_Ch2 = arbitrary_V0CC(Aframe, Ch2frame, Use);
if rem(size(ICS2DCCorr_Ch2_GFP,1),2)~=0
ICS2DCCorr_Ch2_GFP=ICS2DCCorr_Ch2_GFP(2:end,2:end,:); 
end
if rem(size(ICS2DCCorr_Ch2_YFP,1),2)~=0
ICS2DCCorr_Ch2_YFP=ICS2DCCorr_Ch2_YFP(2:end,2:end,:); 
end
if rem(size(ICS2DCCorr_GFP_YFP,1),2)~=0
ICS2DCCorr_GFP_YFP=ICS2DCCorr_GFP_YFP(2:end,2:end,:); 
end
if rem(size(ICS2DCCorr_GFP_A,1),2)~=0
ICS2DCCorr_GFP_A=ICS2DCCorr_GFP_A(2:end,2:end,:); 
end
if rem(size(ICS2DCCorr_A_YFP,1),2)~=0
ICS2DCCorr_A_YFP=ICS2DCCorr_A_YFP(2:end,2:end,:); 
end
if rem(size(ICS2DCCorr_A_Ch2,1),2)~=0
ICS2DCCorr_A_Ch2=ICS2DCCorr_A_Ch2(2:end,2:end,:); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fitting the Ch2 vs. GFP CC curve and plotting figures

[fitergs, befrem, corr, ICS2DCCorrCrop]=CCfit_shift(ICS2DCCorr_Ch2_GFP, pixelsize, linetime, pixeltime, shotnoisecorrection, Ch2nonfilteredframe, a1, lowlim, up);
scrsz = get(0,'ScreenSize');
fiterg2=fitergs(1,:);Ch2GFPfiterg3=fitergs(2,:); %non-weighted fit output parameters
Ch2GFPfiterg2=fitergs(1,:); %fit output parameters
fiterg3=fitergs(2,:);
fig2=figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]);
meann=(mean(ICS2DCCorrCrop,3));
center=(size(ICS2DCCorrCrop,1)+1)/2;
meann(center,center)=0;


subplot(2,6,3)
s=surf(meann);
axis tight
colormap(jet)

zlabel('CC C-G (\xi,\eta)','FontSize',14)
set(s,'LineStyle','none')
title('CC Ch2-GFP')
                                        

subplot(2,6,9)


befremavg=mean(befrem,3);
befremavg(befremavg==max(befremavg))=0;
befremavg(befremavg==min(befremavg))=0;

plot(befremavg,'ko')
avgcorr=mean(corr,3);

hold on

plot(globalsalvo3drics_hor_shift(fiterg3(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'g-')   %not weighted
plot(globalsalvo3drics_hor_shift(fiterg2(1,:),[avgcorr(:,1:3) avgcorr(:,5:6) ]),'b-')


axis([size(avgcorr,1)/2-3*size(ICS2DCCorrCrop,1) size(avgcorr,1)/2+2*size(ICS2DCCorrCrop,1) min(avgcorr(:,4)) max(avgcorr(:,4))])
t=title( '(blue weighted fit, green not weighted fit)');
                                        
set(t, 'Position', get(t,'Position')+ [0 0.001 0 ])                                  
                                        
if filter==1                                       
xlabel(['AC=' num2str(filterwindow/(filterwindow-1)/fiterg2(2)) '   AC_n_w='  num2str(filterwindow/(filterwindow-1)/fiterg3(2))] )                                
else
xlabel(['AC=' num2str(1/fiterg2(2)) '   AC_n_w='  num2str(1/fiterg3(2))] )                                
end         
        
                                        
%fitting the Ch2 vs. YFP CC curve and plotting figures                           
                                        
[fitergs, befrem, corr, ICS2DCCorrCrop]=CCfit_shift(ICS2DCCorr_Ch2_YFP, pixelsize, linetime, pixeltime, shotnoisecorrection, Ch2nonfilteredframe, a1, lowlim, up);
fiterg2=fitergs(1,:);Ch2YFPfiterg3=fitergs(2,:);Ch2YFPfiterg2=fitergs(1,:);

fiterg3=fitergs(2,:);

meann=(mean(ICS2DCCorrCrop,3));
center=(size(ICS2DCCorrCrop,1)+1)/2;
meann(center,center)=0;

subplot(2,6,5)
s=surf(meann);
axis tight
colormap(jet)
zlabel('CC C-Y (\xi,\eta)','FontSize',14)
set(s,'LineStyle','none')
title('CC Ch2-YFP')
                                        

subplot(2,6,11)


befremavg=mean(befrem,3);
befremavg(befremavg==max(befremavg))=0;
befremavg(befremavg==min(befremavg))=0;


plot(befremavg,'ko')
avgcorr=mean(corr,3);

hold on

plot(globalsalvo3drics_hor_shift(fiterg3(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'g-')   %not weighted
plot(globalsalvo3drics_hor_shift(fiterg2(1,:),[avgcorr(:,1:3) avgcorr(:,5:6) ]),'b-')


axis([size(avgcorr,1)/2-3*size(ICS2DCCorrCrop,1) size(avgcorr,1)/2+2*size(ICS2DCCorrCrop,1) min(avgcorr(:,4)) max(avgcorr(:,4))])
                   
                                        
if filter==1                                       
xlabel(['AC=' num2str(filterwindow/(filterwindow-1)/fiterg2(2)) '   AC_n_w='  num2str(filterwindow/(filterwindow-1)/fiterg3(2))] )                                
else
xlabel(['AC=' num2str(1/fiterg2(2)) '   AC_n_w='  num2str(1/fiterg3(2))] )                                
end                                 
                                        
                                        
                                        
                                        
                                        
                                        
 %fitting the GFP vs. YFP CC curve and plotting figures            
[fitergs, befrem, corr, ICS2DCCorrCrop]=CCfit_shift(ICS2DCCorr_GFP_YFP, pixelsize, linetime, pixeltime, shotnoisecorrection, Ch2nonfilteredframe, a1, lowlim, up);
fiterg2=fitergs(1,:);GFPYFPfiterg3=fitergs(2,:);GFPYFPfiterg2=fitergs(1,:);

fiterg3=fitergs(2,:);
meann=(mean(ICS2DCCorrCrop,3));
center=(size(ICS2DCCorrCrop,1)+1)/2;
meann(center,center)=0;

subplot(2,6,1)
s=surf(meann);
axis tight
colormap(jet)
zlabel('CC G-Y (\xi,\eta)','FontSize',14)
set(s,'LineStyle','none')
title('CC GFP-YFP')
                                        

subplot(2,6,7)


befremavg=mean(befrem,3);
befremavg(befremavg==max(befremavg))=0;
befremavg(befremavg==min(befremavg))=0;

plot(befremavg,'ko')
avgcorr=mean(corr,3);

hold on

plot(globalsalvo3drics_hor_shift(fiterg3(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'g-')   %not weighted
plot(globalsalvo3drics_hor_shift(fiterg2(1,:),[avgcorr(:,1:3) avgcorr(:,5:6) ]),'b-')

axis([size(avgcorr,1)/2-3*size(ICS2DCCorrCrop,1) size(avgcorr,1)/2+2*size(ICS2DCCorrCrop,1) min(avgcorr(:,4)) max(avgcorr(:,4))])

                                        
if filter==1                                       
xlabel(['AC=' num2str(filterwindow/(filterwindow-1)/fiterg2(2)) '   AC_n_w='  num2str(filterwindow/(filterwindow-1)/fiterg3(2))] )                                
else
xlabel(['AC=' num2str(1/fiterg2(2)) '   AC_n_w='  num2str(1/fiterg3(2))] )                                
end                                     
                                        
                                         
                                        
 %fitting the GFP vs. A CC curve and plotting figures    
                                        
[fitergs, befrem, corr, ICS2DCCorrCrop]=CCfit_shift(ICS2DCCorr_GFP_A, pixelsize, linetime, pixeltime, shotnoisecorrection, Ch2nonfilteredframe, a1, lowlim, up);
fiterg2=fitergs(1,:);GFPAfiterg3=fitergs(2,:);GFPAfiterg2=fitergs(1,:);

fiterg3=fitergs(2,:);
meann=(mean(ICS2DCCorrCrop,3));
center=(size(ICS2DCCorrCrop,1)+1)/2;
meann(center,center)=0;
subplot(2,6,2)
s=surf(meann);
axis tight
colormap(jet)
zlabel('CC G-A (\xi,\eta)','FontSize',14)
set(s,'LineStyle','none')
title('CC GFP-A')
                                        

subplot(2,6,8)

befremavg=mean(befrem,3);
befremavg(befremavg==max(befremavg))=0;
befremavg(befremavg==min(befremavg))=0;

plot(befremavg,'ko')
avgcorr=mean(corr,3);

hold on

plot(globalsalvo3drics_hor_shift(fiterg3(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'g-')   %not weighted
plot(globalsalvo3drics_hor_shift(fiterg2(1,:),[avgcorr(:,1:3) avgcorr(:,5:6) ]),'b-')


axis([size(avgcorr,1)/2-3*size(ICS2DCCorrCrop,1) size(avgcorr,1)/2+2*size(ICS2DCCorrCrop,1) min(avgcorr(:,4)) max(avgcorr(:,4))])

                                        
if filter==1                                       
xlabel(['AC=' num2str(filterwindow/(filterwindow-1)/fiterg2(2)) '   AC_n_w='  num2str(filterwindow/(filterwindow-1)/fiterg3(2))] )                                
else
xlabel(['AC=' num2str(1/fiterg2(2)) '   AC_n_w='  num2str(1/fiterg3(2))] )                                
end                                     
                                        
                                        
                                                      
                                        
 %fitting the A vs. YFP CC curve and plotting figures    
                                        
[fitergs, befrem, corr, ICS2DCCorrCrop]=CCfit_shift(ICS2DCCorr_A_YFP, pixelsize, linetime, pixeltime, shotnoisecorrection, Ch2nonfilteredframe, a1, lowlim, up);
fiterg2=fitergs(1,:);AYFPfiterg3=fitergs(2,:);AYFPfiterg2=fitergs(1,:);

fiterg3=fitergs(2,:);
meann=(mean(ICS2DCCorrCrop,3));
center=(size(ICS2DCCorrCrop,1)+1)/2;
meann(center,center)=0;
subplot(2,6,4)
s=surf(meann);
axis tight
colormap(jet)
zlabel('CC Y-A (\xi,\eta)','FontSize',14)
set(s,'LineStyle','none')
title('CC Y-A')
                                        

subplot(2,6,10)


befremavg=mean(befrem,3);
befremavg(befremavg==max(befremavg))=0;
befremavg(befremavg==min(befremavg))=0;

plot(befremavg,'ko')
avgcorr=mean(corr,3);

hold on

plot(globalsalvo3drics_hor_shift(fiterg3(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'g-')   %not weighted
plot(globalsalvo3drics_hor_shift(fiterg2(1,:),[avgcorr(:,1:3) avgcorr(:,5:6) ]),'b-')

axis([size(avgcorr,1)/2-3*size(ICS2DCCorrCrop,1) size(avgcorr,1)/2+2*size(ICS2DCCorrCrop,1) min(avgcorr(:,4)) max(avgcorr(:,4))])

                                        
if filter==1                                       
xlabel(['AC=' num2str(filterwindow/(filterwindow-1)/fiterg2(2)) '   AC_n_w='  num2str(filterwindow/(filterwindow-1)/fiterg3(2))] )                                
else
xlabel(['AC=' num2str(1/fiterg2(2)) '   AC_n_w='  num2str(1/fiterg3(2))] )                                
end                                     
                                        

 %fitting the A vs. Ch2 CC curve and plotting figures    
[fitergs, befrem, corr, ICS2DCCorrCrop]=CCfit_shift(ICS2DCCorr_A_Ch2, pixelsize, linetime, pixeltime, shotnoisecorrection, Ch2nonfilteredframe, a1, lowlim, up);
scrsz = get(0,'ScreenSize');
fiterg2=fitergs(1,:);ACh2fiterg3=fitergs(2,:);ACh2fiterg2=fitergs(1,:);
fiterg3=fitergs(2,:);
meann=(mean(ICS2DCCorrCrop,3));
center=(size(ICS2DCCorrCrop,1)+1)/2;
meann(center,center)=0;

subplot(2,6,6)
s=surf(meann);
axis tight
colormap(jet)
zlabel('CC C-A (\xi,\eta)','FontSize',14)
set(s,'LineStyle','none')
title('CC Ch2-A')
                                        
                                      

subplot(2,6,12)


befremavg=mean(befrem,3);
befremavg(befremavg==max(befremavg))=0;
befremavg(befremavg==min(befremavg))=0;

plot(befremavg,'ko')
avgcorr=mean(corr,3);

hold on

plot(globalsalvo3drics_hor_shift(fiterg3(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'g-')   %not weighted
plot(globalsalvo3drics_hor_shift(fiterg2(1,:),[avgcorr(:,1:3) avgcorr(:,5:6) ]),'b-')

axis([size(avgcorr,1)/2-3*size(ICS2DCCorrCrop,1) size(avgcorr,1)/2+2*size(ICS2DCCorrCrop,1) min(avgcorr(:,4)) max(avgcorr(:,4))])

                                        
if filter==1                                       
xlabel(['AC=' num2str(filterwindow/(filterwindow-1)/fiterg2(2)) '   AC_n_w='  num2str(filterwindow/(filterwindow-1)/fiterg3(2))] )                                
else
xlabel(['AC=' num2str(1/fiterg2(2)) '   AC_n_w='  num2str(1/fiterg3(2))] )                                
end                                     
                                        
                                        
  
export_fig (pathpdf, '-append' )


try
note{1,1}; %#ok<VUNUS>
catch
note{1,1}=' ';
end
if filter==1
Ch2GFPfiterg2(2)=Ch2GFPfiterg2(2)/filterwindow*(filterwindow-1);
GFPYFPfiterg2(2)=GFPYFPfiterg2(2)/filterwindow*(filterwindow-1);
Ch2YFPfiterg2(2)=Ch2YFPfiterg2(2)/filterwindow*(filterwindow-1);
end
if ftell(fid) <1
fprintf(fid, 'ID\t FileName\t Note\t Ch2_GFP_CC\t qual.control_w\t qual.control_nw\t Ch2_YFP_CC\t qual.control_w\t qual.control_nw\t GFP_A_CC\t qual.control_w\t qual.control_nw\t A_YFP_CC\t qual.control_w\t qual.control_nw\t Ch2_A_CC\t qual.control_w\t qual.control_nw\t GFP_YFP_CC\t qual.control_w\t qual.control_nw\n');

end
fprintf(fid,'%e\t %s\t %s\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',kkk, [pathname filename],note{1,1}, 1/Ch2GFPfiterg2(2),Ch2GFPfiterg2(6),Ch2GFPfiterg3(6),...
    1/Ch2YFPfiterg2(2),Ch2YFPfiterg2(6),Ch2YFPfiterg3(6),...
    1/GFPAfiterg2(2),GFPAfiterg2(6),GFPAfiterg3(6),...
    1/AYFPfiterg2(2),AYFPfiterg2(6),AYFPfiterg3(6),...
    1/ACh2fiterg2(2),ACh2fiterg2(6),ACh2fiterg3(6),...
    1/GFPYFPfiterg2(2), GFPYFPfiterg2(6),GFPYFPfiterg3(6));    



fclose all;
close all
clearvars -except kkk listofbatches listofbatches2 batchesdir



end
 
AjustACamplitude_4colors(batchesdir) %this script will perform AC analysis for the same data
