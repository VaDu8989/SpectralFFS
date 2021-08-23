function [speclinefluorescenceseries,speclinearrayalignedba] = speclinealignfunc(speclinearraymasked)
global  blocksize spatialfilter membranewidth

% Block average
fprintf('Block average...\n');
linearraymasked=sum(speclinearraymasked,3);
membranepositions=zeros(length(linearraymasked(:,1))/blocksize,1);
for i=1:length(membranepositions)
    linearrayblock=sum(linearraymasked((i-1)*blocksize+1:i*blocksize,:),1);
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

% Fit average line profile with Gauss function (optionally plus heavyside
% function -> function g)
    %ustep = @(x,range) range>=x;
    %g=fittype(@(a1,a2,a3,b1,c1,x) a1*exp(-((x-b1)./c1).^2)+a2*ustep(b1,x)+a3);
    %coeffnames(g);
    %membfit=fit([1:1:length(linearrayalignedmasum)]',linearrayalignedmasum'./sum(linearrayalignedmasum),g,'StartPoint',[1 0.1 0 64 1]);
membfit=fit([1:1:length(linearrayalignedmasum)]',linearrayalignedmasum'./sum(linearrayalignedmasum),'gauss1');
mu=membfit.b1;
sigma=membfit.c1/sqrt(2);
indexmean=round(membfit.b1);
indexcutlow=round(mu-spatialfilter*sigma);
indexcutup=round(mu+spatialfilter*sigma);
membranewidth=indexcutup-indexcutlow-1;

% Optional: Plot of line profile and fit
    % scrsz=   get(0,'ScreenSize');
    % figure('OuterPosition',[2*scrsz(3)/3 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2],'Name','Total Membrane fluorescence peak')
    % plot(1:length(linearrayalignedmasum),linearrayalignedmasum./sum(linearrayalignedmasum),'-b')
    % hold on
    % plot(membfit,'-r')
    % xlabel('Pixel Number')
    % ylabel('Value')

% filtering of membrane region
speclinearrayalignedba(:,1:indexcutlow,:)=NaN;
speclinearrayalignedba(:,indexcutup:end,:)=NaN;

% Calculation of time series
speclinefluorescenceseries=nansum(speclinearrayalignedba,2);