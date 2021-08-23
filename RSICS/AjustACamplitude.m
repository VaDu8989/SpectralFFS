function AjustACamplitude(batchesdir)  % This function is the analogue of the script ARICCS_Spectraldecomposition.m but it calculates AC instead of CC. For details, please see the comments and explanations in ARICCS_Spectraldecomposition.m
simulationhappening=0;


if simulationhappening==0

listofbatches=dir([batchesdir '/*RICSbatch*']);
kk=1;
listofbatches2=zeros(size(listofbatches,1));
for k=1:size(listofbatches,1)
if listofbatches(k,1).name(end-12:end-10)=='YFP' %#ok<BDSCA>
listofbatches2(kk)=listofbatches(k);
kk=kk+1;
end
end
else
listofbatches=dir([batchesdir '/*total_sim*']);
pathname=batchesdir;
filename='simulation';
kk=1; %#ok<NASGU>
listofbatches2=listofbatches;
end
 

for kkk=1:size(listofbatches2,2)
kkk %#ok<NOPRT>
'out of' %#ok<NOPRT>
size(listofbatches2,2)
fid = fopen([batchesdir '\ACanalisi ARICS' date '.txt'],'a');%where it writes results
pathpdf = [ batchesdir '\ACanalisi ARICS' date '.tiff' ];
load([batchesdir '\' listofbatches2(1,kkk).name]); %#ok<LOAD>
lowlim= [0.17 0 0 -0.00001 3 ];  %%lower limits for fit:  w0   N   D   offset S
up = [0.3 1000 100 0.00001 7];  %upper limits
a1 = [0.20 20 20 0 5]; %starting parameters


%This section is needed only for simulations
% if simulationhappening==0;
%                                                                                                      
% data=load ([pathname  Ch2filename]);
% nameoffield=fieldnames(data);
% namefield=nameoffield{1,1};  
%  
% Ch2simulationresult=data.(namefield);
%     
%                                                                                                          
% data=load ([pathname  GFPfilename]);
% nameoffield=fieldnames(data);
% namefield=nameoffield{1,1};  
%  
%  GFPsimulationresult=data.(namefield);
%     
%                                                                                                          
% data=load ([pathname  YFPfilename]);
% nameoffield=fieldnames(data);
% namefield=nameoffield{1,1};  
%  
% YFPsimulationresult=data.(namefield);
% 
% Ch2frame=double(Ch2simulationresult(:,:,startframe:endframe));
% GFPframe=double(GFPsimulationresult(:,:,startframe:endframe));
% YFPframe=double(YFPsimulationresult(:,:,startframe:endframe));
% end




% if simulationhappening==1
%     try
%  Ch2frame=tempframe3;
% GFPframe=tempframe1;
% YFPframe=tempframe2; 
%     catch
%          load([batchesdir '\Simulation_parameters.mat']);
%     Ch2frame=frame3;
% GFPframe=frame1;
% YFPframe=frame2;   
% endframe=size(frame1,3);
% startframe=1;
% filterwindow=0;
% 
% Use=zeros(size(frame1,1),size(frame1,2),1)+1;
% shotnoisecorrection=1;
%     end
%     filter=0;
% end

    
Usenan=double(Use);
Usenan(Usenan==0)=NaN;
clear frame

%%%%%%%%%%%%%%%%%%filtering slow dynamics, moving average (until line 130)
if endframe-startframe+1<filterwindow+1 && filter==1
   error('Error! You are trying to apply a running average with a window larger than the total amount of frames! \n Remove the filterning by setting filter=0 \n') 
end



Ch2nonfilteredframe=Ch2frame;
GFPnonfilteredframe=GFPframe;
YFPnonfilteredframe=YFPframe;


Ch2firstframe=nanmean(nanmean(Ch2frame(:,:,1).*Usenan));
Ch2lastframe=nanmean(nanmean(Ch2frame(:,:,end).*Usenan));
GFPfirstframe=nanmean(nanmean(GFPframe(:,:,1).*Usenan));
GFPlastframe=nanmean(nanmean(GFPframe(:,:,end).*Usenan));
YFPfirstframe=nanmean(nanmean(YFPframe(:,:,1).*Usenan));
YFPlastframe=nanmean(nanmean(YFPframe(:,:,end).*Usenan));


Ch2totalmean=mean(mean(mean(Ch2frame(:,:,:))));
GFPtotalmean=mean(mean(mean(GFPframe(:,:,:))));
YFPtotalmean=mean(mean(mean(YFPframe(:,:,:))));  
if filter==1 && endframe-startframe+1>filterwindow+1


Ch2frame2=RICSfilter(Ch2frame, filterwindow, startframe, endframe);
GFPframe2=RICSfilter(GFPframe, filterwindow, startframe, endframe);
YFPframe2=RICSfilter(YFPframe, filterwindow, startframe, endframe);
      
clear frame Ch2frame GFPframe YFPframe
Ch2frame=Ch2frame2+Ch2totalmean;
GFPframe=GFPframe2+GFPtotalmean;
YFPframe=YFPframe2+YFPtotalmean;
clear frame2 Ch2frame2 GFPframe2 YFPframe2
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Calculation of autocorrelation between different channels, based on
%ARICS algorithm (DOI: 10.1016/j.bpj.2016.09.012)
ICS2DCCorr_Ch2_AC = arbitrary_salvoV0(Ch2frame, Use);
if rem(size(ICS2DCCorr_Ch2_AC,1),2)==0
   ICS2DCCorr_Ch2_AC=ICS2DCCorr_Ch2_AC(2:end,2:end,:); 
end
        
ICS2DCCorr_GFP_AC = arbitrary_salvoV0(GFPframe, Use);

if rem(size(ICS2DCCorr_GFP_AC,1),2)==0
   ICS2DCCorr_GFP_AC=ICS2DCCorr_GFP_AC(2:end,2:end,:); 
end

ICS2DCCorr_AC_YFP = arbitrary_salvoV0(YFPframe, Use);

if rem(size(ICS2DCCorr_AC_YFP,1),2)==0
   ICS2DCCorr_AC_YFP=ICS2DCCorr_AC_YFP(2:end,2:end,:); 
end


%fitting the Ch2 AC curve and plotting figures


                                        [fitergs, befrem, corr, ICS2DCCorrCrop]=CCfit(ICS2DCCorr_Ch2_AC, pixelsize, linetime, pixeltime, shotnoisecorrection, Ch2nonfilteredframe, a1, lowlim, up);
                                        scrsz = get(0,'ScreenSize');
%fitergs --> first row weighted, second row non-weighted
 fiterg2=fitergs(1,:);%weighted
fiterg3=fitergs(2,:);Ch2fiterg2=fitergs(1,:);
Ch2fitergnw=fitergs(2,:); %#ok<NASGU> %non-weighted
D_Ch2=Ch2fiterg2(3);
fig2=figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]); %#ok<NASGU>
                                        
meann=(mean(ICS2DCCorrCrop,3));
center=(size(ICS2DCCorrCrop,1)+1)/2;

meann(center,center)=0;


subplot(2,3,3)
s=surf(meann);
axis tight
colormap(jet)
zlabel('AC C (a.u.)','FontSize',14)
set(s,'LineStyle','none')
title('AC Ch2')
                                        
                                        
subplot(2,3,6)

befremavg=mean(befrem,3);
befremavg(befremavg==max(befremavg))=0;
befremavg(befremavg==min(befremavg))=0;
plot(befremavg,'ko')
avgcorr=mean(corr,3);

hold on

plot(globalsalvo3drics(fiterg3(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'g-')   %not weighted
plot(globalsalvo3drics(fiterg2(1,:),[avgcorr(:,1:3) avgcorr(:,5:6) ]),'b-')

axis([size(avgcorr,1)/2-3*size(ICS2DCCorrCrop,1) size(avgcorr,1)/2+2*size(ICS2DCCorrCrop,1) min(avgcorr(:,4)) max(avgcorr(:,4))])
t=title( '(blue weighted fit, green not weighted fit)');
set(t, 'Position', get(t,'Position')+ [0 0.001 0 ])                                  
if filter==1                                       
xlabel(['AC=' num2str(filterwindow/(filterwindow-1)/fiterg2(2)) '   AC_n_w='  num2str(filterwindow/(filterwindow-1)/fiterg3(2))] )                                
else
xlabel(['AC=' num2str(1/fiterg2(2)) '   AC_n_w='  num2str(1/fiterg3(2))] )                                
end
                                        
                                        
%fitting the GFP AC curve and plotting figures                                        
                                        
                                        
[fitergs, befrem, corr, ICS2DCCorrCrop]=CCfit(ICS2DCCorr_GFP_AC, pixelsize, linetime, pixeltime, shotnoisecorrection, GFPnonfilteredframe, a1, lowlim, up);
fiterg2=fitergs(1,:);fiterg3=fitergs(2,:);GFPfiterg2=fitergs(1,:);fiterg5=fitergs(1,:); %#ok<NASGU>
GFPfitergnw=fitergs(2,:); %#ok<NASGU> %non-weighted
D_GFP=GFPfiterg2(3);

meann=(mean(ICS2DCCorrCrop,3));
center=(size(ICS2DCCorrCrop,1)+1)/2;

meann(center,center)=0;
subplot(2,3,1)
s=surf(meann);
axis tight
colormap(jet)
zlabel('AC G (a.u.)','FontSize',14)
set(s,'LineStyle','none')
title('AC GFP')

subplot(2,3,4)


befremavg=mean(befrem,3);
befremavg(befremavg==max(befremavg))=0;
befremavg(befremavg==min(befremavg))=0;


plot(befremavg,'ko')
avgcorr=mean(corr,3);

hold on

plot(globalsalvo3drics(fiterg3(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'g-')   %not weighted
plot(globalsalvo3drics(fiterg2(1,:),[avgcorr(:,1:3) avgcorr(:,5:6) ]),'b-')   % weighted

fitergslow=fiterg3;
fitergslow(3)=10^-10;
plot(globalsalvo3drics(fitergslow(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'r-')   %plotting a curve with the slowest possible diffusion coefficent

axis([size(avgcorr,1)/2-3*size(ICS2DCCorrCrop,1) size(avgcorr,1)/2+2*size(ICS2DCCorrCrop,1) min(avgcorr(:,4)) max(avgcorr(:,4))])
title('Red curve is slowest possible model')                           
                           
if filter==1                                       
xlabel(['AC=' num2str(filterwindow/(filterwindow-1)/fiterg2(2)) '   AC_n_w='  num2str(filterwindow/(filterwindow-1)/fiterg3(2))] )                                
else
xlabel(['AC=' num2str(1/fiterg2(2)) '   AC_n_w='  num2str(1/fiterg3(2))] )                                
end                       
                                        
                                        
                                        
                                        
                                        
                                        
                                        
%fitting the YFP AC curve and plotting figures                                                
                                        
                                        
[fitergs, befrem, corr, ICS2DCCorrCrop]=CCfit(ICS2DCCorr_AC_YFP, pixelsize, linetime, pixeltime, shotnoisecorrection, YFPnonfilteredframe, a1, lowlim, up);
fiterg2=fitergs(1,:);fiterg3=fitergs(2,:);YFPfiterg2=fitergs(1,:);  %weighted                                     
YFPfitergnw=fitergs(2,:); %#ok<NASGU> %non-weighted
D_YFP=YFPfiterg2(3);
meann=(mean(ICS2DCCorrCrop,3));
center=(size(ICS2DCCorrCrop,1)+1)/2;

meann(center,center)=0;
subplot(2,3,2)
s=surf(meann);
axis tight
colormap(jet)
zlabel('AC Y (a.u.)','FontSize',14)
set(s,'LineStyle','none')
title('AC YFP')
                                        
                                        
                                        

subplot(2,3,5)


befremavg=mean(befrem,3);
befremavg(befremavg==max(befremavg))=0;
befremavg(befremavg==min(befremavg))=0;

plot(befremavg,'ko')
avgcorr=mean(corr,3);

hold on

plot(globalsalvo3drics(fiterg3(1,:),[avgcorr(:,1:3) avgcorr(:,5:6)]),'g-')   %not weighted
plot(globalsalvo3drics(fiterg2(1,:),[avgcorr(:,1:3) avgcorr(:,5:6) ]),'b-')


axis([size(avgcorr,1)/2-3*size(ICS2DCCorrCrop,1) size(avgcorr,1)/2+2*size(ICS2DCCorrCrop,1) min(avgcorr(:,4)) max(avgcorr(:,4))])
title(['D YFP=' num2str(D_YFP)])
                                        
if filter==1                                       
xlabel(['AC=' num2str(filterwindow/(filterwindow-1)/fiterg2(2)) '   AC_n_w='  num2str(filterwindow/(filterwindow-1)/fiterg3(2))] )                                
else
xlabel(['AC=' num2str(1/fiterg2(2)) '   AC_n_w='  num2str(1/fiterg3(2))] )                                
end             
                                        
                                        

%Exporting figures and writing results

export_fig (pathpdf, '-append' )


try
    note{1,1}; %#ok<VUNUS>
catch
   note{1,1}=' ';
end
if filter==1
Ch2fiterg2(2)=Ch2fiterg2(2)/filterwindow*(filterwindow-1);
GFPfiterg2(2)=GFPfiterg2(2)/filterwindow*(filterwindow-1);
YFPfiterg2(2)=YFPfiterg2(2)/filterwindow*(filterwindow-1);
end
if ftell(fid) <1
fprintf(fid, 'ID\t FileName\t Note\t Ch2_AC\t Ch2_D\t Ch2_B\t Ch2_Mean I\t Ch2_Mean Frame1\t Ch2_Mean last Frame\t GFP_AC\t GFP_D\t GFP_B\t GFP_Mean I\t GFP_Mean Frame1\t GFP_Mean last Frame\t YFP_AC\t YFP_D\t YFP_B\t YFP_Mean I\t YFP_Mean Frame1\t  YFP_Mean last Frame\n');

end
fprintf(fid,'%e\t %s\t %s\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',kkk, [pathname filename],note{1,1}, 1/Ch2fiterg2(2),D_Ch2, 1/Ch2fiterg2(2)* Ch2firstframe/1000/pixeltime ,nanmean(nanmean(nanmean(Ch2frame.*Usenan)))/1000/pixeltime, Ch2firstframe/1000/pixeltime, Ch2lastframe/1000/pixeltime, ...
    1/GFPfiterg2(2),D_GFP, 1/GFPfiterg2(2)* GFPfirstframe/1000/pixeltime ,nanmean(nanmean(nanmean(GFPframe.*Usenan)))/1000/pixeltime, GFPfirstframe/1000/pixeltime, GFPlastframe/1000/pixeltime,...
    1/YFPfiterg2(2),D_YFP, 1/YFPfiterg2(2)* YFPfirstframe/1000/pixeltime ,nanmean(nanmean(nanmean(YFPframe.*Usenan)))/1000/pixeltime, YFPfirstframe/1000/pixeltime, YFPlastframe/1000/pixeltime);    



fclose all;
close all
clearvars -except kkk listofbatches listofbatches2 batchesdir simulationhappening


end
 

