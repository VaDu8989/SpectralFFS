function [speclinefluorescenceseries,speclinearrayalignedba] = speclinealignfuncbr(speclinearraymasked,specbackgroundmasked)
global  blocksize spatialfilter membranewidth

% Block average
fprintf('Block average...\n');
linearraymasked=sum(speclinearraymasked,3);
membranepositions=zeros(length(linearraymasked(:,1))/blocksize,1);
backgroundblocks=zeros(length(linearraymasked(:,1))/blocksize,size(speclinearraymasked,3));
for i=1:length(membranepositions)
    linearrayblock=sum(linearraymasked((i-1)*blocksize+1:i*blocksize,:),1);
    backgroundblocks(i,:)=squeeze(mean(mean(specbackgroundmasked((i-1)*blocksize+1:i*blocksize,:,:),2),1));
    [maxblock,membranepositionblock]=max(linearrayblock);
    membranepositions(i)=membranepositionblock;
end

% Alignment of scan lines (maximum fluorescence)
fprintf('Aligning lines...\n');
maxdriftba=max(membranepositions)-min(membranepositions);
speclinearrayalignedba=zeros(size(speclinearraymasked,1),size(speclinearraymasked,2)+maxdriftba,size(speclinearraymasked,3));
firstblockindex=membranepositions(1);
maxdeviation=max(membranepositions)-firstblockindex;
for j=1:length(membranepositions)
    blockdrift=membranepositions(j)-firstblockindex;
    shiftindex=maxdeviation-blockdrift;
    speclinearrayalignedba((j-1)*blocksize+1:j*blocksize,1+shiftindex:shiftindex+length(linearraymasked(1,:)),:)=speclinearraymasked((j-1)*blocksize+1:j*blocksize,:,:);
end

% Sum of all lines and fit
linearrayalignedmasum=sum(sum(speclinearrayalignedba,3),1);

% Fit average line profile with Gauss function plus sigmoid function
gsigmoid=fittype(@(a1,b1,c1,a2,c2,d1,x) a1*exp(-((x-b1)./c1).^2)+c2./(1+exp(-a2*(x-b1)))+d1);    
membfit=fit([1:1:length(linearrayalignedmasum)]',linearrayalignedmasum'./sum(linearrayalignedmasum),gsigmoid,'StartPoint',[1 100 1 2 0.05 0.01]);

mu=membfit.b1;
sigma=membfit.c1/sqrt(2);
indexmean=round(membfit.b1);

indexcutlow=round(mu-spatialfilter*sigma);
indexcutup=round(mu+spatialfilter*sigma);
membranewidth=indexcutup-indexcutlow-1;

% Optional: Plot of line profile and fit
scrsz=   get(0,'ScreenSize');
figure('OuterPosition',[2*scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Total Membrane fluorescence peak')
plot(1:length(linearrayalignedmasum),linearrayalignedmasum./sum(linearrayalignedmasum),'-b')
hold on
plot(membfit,'-r')
xlabel('Pixel Number')
ylabel('Value')
% fitprofile=gsigmoid(membfit.a1,membfit.b1,membfit.c1,membfit.a2,membfit.c2,membfit.d1,1:length(linearrayalignedmasum));
% profile=linearrayalignedmasum./sum(linearrayalignedmasum);
% save('fitprofile.mat','fitprofile');
% save('profile.mat','profile');    

% filtering of membrane region
speclinearrayalignedba(:,1:indexcutlow,:)=NaN;
speclinearrayalignedba(:,indexcutup:end,:)=NaN;

% Calculation of time series
speclinefluorescenceseries=nansum(speclinearrayalignedba,2);

% Subtraction of background signal (on one side)
for j=1:length(membranepositions)
    speclinefluorescenceseries((j-1)*blocksize+1:j*blocksize,:)=speclinefluorescenceseries((j-1)*blocksize+1:j*blocksize,:)-0.5*backgroundblocks(i,:);
end