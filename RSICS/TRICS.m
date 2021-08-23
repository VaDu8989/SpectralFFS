clear  
close all
simulationhappening=0;
batchesdir=uigetdir('\\CELLMEMBIOPHYS2\Shared\');

% loading batch files
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
kk=1;

listofbatches2=listofbatches;
end
 

for kkk=1:size(listofbatches2,2)

fid = fopen([batchesdir '\CCCanalisi.txt'],'a');%where it writes results
pathpdf = [ batchesdir '\CCCanalisi.pdf' ];

load([batchesdir '\' listofbatches2(1,kkk).name]);
     
    
lowlim= [0.19 0 0 -0.0001 3 0];  %%lower limit for fit:  w0   N   D   offset S

up = [0.3 100000 100 0.0001 7 13];  %upper limit
a1 = [0.25 500 50 0 5 11]; %starting parameters



filter=1;  % 1 YES, 0 NO   = Filtering slow dynamics away   
filterwindow=4;
%the three channels are referred to as Ch2, GFP, YFP
if simulationhappening==0
                                                                                                     
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
    
    
Ch2frame=double(Ch2acq(:,:,startframe:endframe));
GFPframe=double(GFPacq(:,:,startframe:endframe));
YFPframe=double(YFPacq(:,:,startframe:endframe));

end


clear frame



%%%%%%%%%%%%%%%%%%filtering slow dynamics, moving average 
if endframe-startframe+1<filterwindow+1 && filter==1
error('Error! You are trying to apply a running average with a window larger than the total amount of frames! \n Remove the filterning by setting filter=0 \n') 
end



Ch2nonfilteredframe=Ch2frame;
GFPnonfilteredframe=GFPframe;
YFPnonfilteredframe=YFPframe;

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


Ch2frame=Ch2frame(1:size(Ch2frame,1)-1,1:size(Ch2frame,1)-1,:);
GFPframe=GFPframe(1:size(Ch2frame,1),1:size(Ch2frame,1),:);
YFPframe=YFPframe(1:size(Ch2frame,1),1:size(Ch2frame,1),:);
Use=Use(1:size(Ch2frame,1),1:size(Ch2frame,1));

CCC=manualCCC(Ch2frame, GFPframe, YFPframe, Use); %not all points are calculated to avoid shot noise, see lines 42-44 within manualCCC function

CCC=CCC./1;
CCC_complete=manualCCC_complete(Ch2frame, GFPframe, YFPframe, Use); % here more points are calcolated, see lines 42-42 (line 44 is commented here). This approach make sense for simulations without noise, as a test


%at this point, the triple correlation 4-d surface is calculated. Nothing
%else is needed. 
%Everything that follows is not strictly needed and contains only examples/alternatives for how
%to display this data. Please explore and test for your applications: plotCCC, plotCCC_complete,
%plotCCC_alternative, plotCCC_alternative_complete, 
%plotCCC_alt2 and its smoothed version.
%Any feedback is welcome. 

centerCCC=ceil(size(CCC,1)/2);
plotCCC=zeros(int8(((size(CCC,1)-centerCCC)^2+ (size(CCC,1)-centerCCC)^2)^0.5+1),int8(((size(CCC,1)-centerCCC)^2+ (size(CCC,1)-centerCCC)^2)^0.5+1));
plotCCC_complete=plotCCC;
plotCCC_alternative=plotCCC;
normplotCCC=plotCCC; %for normalization
normplotCCC_alternative=plotCCC;
normplotCCC_complete=plotCCC;
plotCCC_alternative_complete=plotCCC;
normplotCCC_alternative_complete=plotCCC;
for i=1: size(CCC,1)
for j=1:size(CCC,2)
for m=1:size(CCC,3)
for n=1:size(CCC,4)
if   plotCCC(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))==0 && isnan(CCC(i,j,m,n))==0
plotCCC(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= CCC(i,j,m,n);
                
normplotCCC(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= normplotCCC(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+1;
else
if isnan(CCC(i,j,m,n))==0
plotCCC(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))=plotCCC(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+CCC(i,j,m,n);
normplotCCC(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= normplotCCC(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+1;
end
end
end
end
end
end

for i=1: size(CCC,1)
for j=1:size(CCC,2)
for m=1:size(CCC,3)
for n=1:size(CCC,4)
if   plotCCC_alternative(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))==0 && isnan(CCC(i,j,m,n))==0
plotCCC_alternative(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= CCC(i,j,m,n);
normplotCCC_alternative(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= normplotCCC_alternative(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+1;
else
if isnan(CCC(i,j,m,n))==0
plotCCC_alternative(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))=plotCCC_alternative(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+CCC(i,j,m,n);
normplotCCC_alternative(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= normplotCCC_alternative(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+1;
end
end
end
end
end
end

for i=1: size(CCC,1)
for j=1:size(CCC,2)
for m=1:size(CCC,3)
for n=1:size(CCC,4)
if   plotCCC_alternative_complete(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))==0 && isnan(CCC_complete(i,j,m,n))==0
plotCCC_alternative_complete(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= CCC_complete(i,j,m,n);
                
normplotCCC_alternative_complete(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= normplotCCC_alternative_complete(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+1;
else
if isnan(CCC_complete(i,j,m,n))==0
plotCCC_alternative_complete(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))=plotCCC_alternative_complete(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+CCC_complete(i,j,m,n);
normplotCCC_alternative_complete(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= normplotCCC_alternative_complete(     int8(((i-centerCCC)^2+ (m-centerCCC)^2)^0.5+1)      ,int8(((j-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+1;
end
end
end
end
end
end



for i=1: size(CCC_complete,1)
for j=1:size(CCC_complete,2)
for m=1:size(CCC_complete,3)
for n=1:size(CCC_complete,4)
if   plotCCC_complete(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))==0 && isnan(CCC_complete(i,j,m,n))==0
plotCCC_complete(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= CCC_complete(i,j,m,n);
                
normplotCCC_complete(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= normplotCCC_complete(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+1;
else
if isnan(CCC_complete(i,j,m,n))==0
plotCCC_complete(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))=plotCCC_complete(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+CCC_complete(i,j,m,n);
normplotCCC_complete(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))= normplotCCC_complete(     int8(((i-centerCCC)^2+ (j-centerCCC)^2)^0.5+1)      ,int8(((m-centerCCC)^2+(n-centerCCC)^2)^0.5+1))+1;
end
end
end
end
end
end


finalplotCCC=plotCCC./normplotCCC;
finalplotCCC_complete=plotCCC_complete./normplotCCC_complete;
finalplotCCC_alternative=plotCCC_alternative./normplotCCC_alternative;
finalplotCCC_alternative_complete=plotCCC_alternative_complete./normplotCCC_alternative_complete;
surf(finalplotCCC);
savefig([ batchesdir '\' GFPfilename '_filter=' num2str(filterwindow*filter)  '_CCC.fig' ])
figure;
surf(finalplotCCC_complete);
savefig([ batchesdir '\' GFPfilename '_filter=' num2str(filterwindow*filter)  '_CCC_complete.fig' ])
figure;
surf(finalplotCCC_alternative);

savefig([ batchesdir '\' GFPfilename '_filter=' num2str(filterwindow*filter)  '_CCC_alternative.fig' ])
figure;
surf(finalplotCCC_alternative_complete);
savefig([ batchesdir '\' GFPfilename '_filter=' num2str(filterwindow*filter)  '_CCC_alternative_complete.fig' ])



finalplotCCC_alt2_complete=plotof4d(CCC_complete);
[~,idx] = sort(finalplotCCC_alt2_complete(:,3)); % s
finalplotCCC_alt2_complete = finalplotCCC_alt2_complete(idx,:);                                   
cm = colormap(jet(size(finalplotCCC_alt2_complete,1)));                                         % Define Colormap
figure
scatter3(finalplotCCC_alt2_complete(:,1),finalplotCCC_alt2_complete(:,2), finalplotCCC_alt2_complete(:,3), [], cm, 'filled')
grid on
view(30,5)
savefig([ batchesdir '\' GFPfilename '_filter=' num2str(filterwindow*filter)  '_CCC_alt2_complete.fig' ])



finalplotCCC_alt2=plotof4d(CCC);
[~,idx] = sort(finalplotCCC_alt2(:,3)); % s
finalplotCCC_alt2 = finalplotCCC_alt2(idx,:);                                   
cm = colormap(jet(size(finalplotCCC_alt2,1)));                                         % Define Colormap
figure
scatter3(finalplotCCC_alt2(:,1),finalplotCCC_alt2(:,2), finalplotCCC_alt2(:,3), [], cm, 'filled')
grid on
view(30,5)
savefig([ batchesdir '\' GFPfilename '_filter=' num2str(filterwindow*filter)  '_CCC_alt2.fig' ])


smoothedplot=finalplotCCC_alt2;
smoothedplot(:,1)=ceil(smoothedplot(:,1));
smoothedplot(:,2)=ceil(smoothedplot(:,2));

finalplotsmooth=zeros(max(smoothedplot(:,1)),max(smoothedplot(:,1)));
finalplotsmooth_n=finalplotsmooth;

for i=1:size(finalplotsmooth,1)
for j=1:size(finalplotsmooth,2)
finalplotsmooth(i,j)=nanmean(smoothedplot((smoothedplot(:,1)== i-1 & smoothedplot(:,2)== j-1),3));

end
end


figure;
surf(finalplotsmooth)
savefig([ batchesdir '\' GFPfilename '_filter=' num2str(filterwindow*filter)  '_CCC_alt2_smoothed.fig' ])

save([ batchesdir '\' GFPfilename '_filter=' num2str(filterwindow*filter) '_CCC.mat' ], 'CCC', 'finalplotCCC','CCC_complete', 'finalplotCCC_complete', 'finalplotCCC_alternative',  'finalplotCCC_alternative_complete','finalplotCCC_alt2_complete','finalplotCCC_alt2','finalplotsmooth' )
end
