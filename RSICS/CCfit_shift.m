function [out, befrem, corr, ICS2DCorrCrop]=CCfit_shift(ICS2DCorr, pixelsize, linetime, pixeltime, shotnoisecorrection, nonfilteredframe, a1, lowlim, up) %#ok<INUSL>

if pixelsize*size(ICS2DCorr,1)>2
ICS2DCorrCrop = autocropcenter(ICS2DCorr,double(2*round(2/pixelsize)/2));   %cropping around the center, obtaining a 2x2µm square
else
ICS2DCorrCrop = ICS2DCorr; 
end


ICS2DCorrCrop2=ICS2DCorrCrop;
center=(size(ICS2DCorrCrop,1)+1)/2;

shotnoisecorrection=1;

%builds correlation matrix as a function of x displacement from center, y
%displacement and n (sequence number) 

corr=zeros(size(ICS2DCorrCrop2,1)*size(ICS2DCorrCrop2,1),7,size(ICS2DCorrCrop2,3));
for k=1:size(ICS2DCorrCrop2,3)
    l=1;
    for i=1:size(ICS2DCorrCrop2,1) %lines

        for j=1:size(ICS2DCorrCrop2,1) %pixels in one line

            corr(l,1,k)=(j-center);  %x
            corr(l,2,k)=abs(i-center);  %y
            corr(l,3,k)=linetime;           
            corr(l,4,k)=ICS2DCorrCrop2(i,j,k);
            corr(l,5,k)=pixelsize;
            corr(l,6,k)=pixeltime;
            corr(l,7,k)=1;
            l=l+1;
        end
    end
end


%Remove shot noise
befrem=corr(:,4,:);
avgcorr=mean(corr,3);

for k=1:shotnoisecorrection
corr(avgcorr(:,4)==(max(max(avgcorr(:,4)))),7,:)=0;  %setting max and min values to 0. Same for their weight. This is not "Shot noise". This unwanted values happen when small arbitrary regions are chosen, due e.g. to noise
corr(find(avgcorr(:,4)==(min(min(avgcorr(:,4))))),7,:)=0; 
corr(avgcorr(:,4)==(max(max(avgcorr(:,4)))),4,:)=0;
corr(avgcorr(:,4)==(min(min(avgcorr(:,4)))),4,:)=0;
corr(avgcorr(:,1)==0 & corr(:,2)==0,7,:)=0; % weight of central point set to 0
corr(avgcorr(:,1)==0 & corr(:,2)==0,4,:)=corr(find(avgcorr(:,1)==0 & corr(:,2)==0)-1,4,:)/2+corr(find(avgcorr(:,1)==0 & corr(:,2)==0)+1,4,:)/2; %interpolation for central point
end
%%%%%%%%%%%

 


avgcorr=mean(corr,3);
avgcorr(:,7)=(((abs(avgcorr(:,1))+1)).^-2).*(((avgcorr(:,2)+1)).^-1); %weight for fitting. This can be adapted
% avgcorr(:,7)=(std(corr(:,4,:),1,3)./mean(corr(:,4,:),3)).^-1; %Other example for weigths
avgcorr(find(avgcorr(:,1)==0 & corr(:,2)==0),7)=0; %#ok<*FNDSB>; again, this implies ignoring the central point in the AC/CC surface


        
                               

%%%%%%%%%%%%%%%%%%%% Fitting the average correlation (weighted and non
%%%%%%%%%%%%%%%%%%%% weighted fit)


%note that these fitting models contain an additional parameter (shift)
%that should converge to 0 if curves are not too noisy
[fiterg2,~,~,~,~,~,~]=lsqcurvefit(@global3dricsweighted_hor_offset,a1,[avgcorr(:,1:3) avgcorr(:,5:6) avgcorr(:,7)],(avgcorr(:,4).*avgcorr(:,7)), lowlim, up,options); 
[fiterg3,~,~,~,~,~,~]=lsqcurvefit(@globals3drics_hor_shift,a1,[avgcorr(:,1:3) avgcorr(:,5:6)],(avgcorr(:,4)), lowlim, up,options); 

out=[fiterg2 ; fiterg3];