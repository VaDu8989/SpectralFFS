% Signal to noise FCCS
clear all
close all

normalize2single=0;

path= uigetdir; 
files_curveincl=dir([path '/*final_curveincls.txt']);
segnumbs=zeros(size(files_curveincl,1),2);
for i=1:size(files_curveincl,1)
    [namedata,remain]=strtok(files_curveincl(i).name,'.');
    namefile=files_curveincl(i).name;
    curveincl=load([path '/' namefile]);
    %SNCh(i)=dataCFs(1,9);
    segnumbs(i,:)=sum(curveincl,1);
end
figure('Name','Number of segments')
hist(segnumbs)

files_SNR=dir([path '/*final_ACFseg.txt']);
SNRs=zeros(size(files_SNR,1),2);
for i=1:size(files_SNR,1)
    [namedata,remain]=strtok(files_SNR(i).name,'.');
    namefile=files_SNR(i).name;
    ACFinfo=load([path '/' namefile]);
    %SNCh(i)=dataCFs(1,9);
    SNRs(i,1)=mean(ACFinfo(:,3));
    SNRs(i,2)=mean(ACFinfo(:,6));
end

figure('Name','SNRs')
plot(SNRs(:,1),'g')
hold on
plot(SNRs(:,2),'y')

files_Is=dir([path '/*final_fitparameters_seg.txt']);
Is=zeros(size(files_Is,1),2);
for i=1:size(files_Is,1)
    [namedata,remain]=strtok(files_Is(i).name,'.');
    namefile=files_Is(i).name;
    fitparameters=load([path '/' namefile]);
    %SNCh(i)=dataCFs(1,9);
    Is(i,1)=fitparameters(12,1);
    Is(i,2)=fitparameters(13,1);
end

GYratio=Is(:,1)./Is(:,2);
YGratio=1./GYratio;

figure('Name','SNR vs. I ratio')
semilogx(GYratio,SNRs(:,1),'gd')
hold
semilogx(GYratio,SNRs(:,2),'ys')
ylabel('SNR')
xlabel('I(G)/I(Y)')

if normalize2single==1
    pathsingle=uigetdir;
    files_SNR=dir([pathsingle '/*final_ACF.txt']);
    SNRssingle=zeros(size(files_SNR,1),1);
    for i=1:size(files_SNR,1)
        [namedata,remain]=strtok(files_SNR(i).name,'.');
        namefile=files_SNR(i).name;
        ACFinfo=load([pathsingle '/' namefile]);
        %SNCh(i)=dataCFs(1,9);
        SNRssingle(i)=mean(ACFinfo(:,3));
    end
    SNRavg=mean(SNRssingle);

    SNRsnorm=SNRs/SNRavg;

    figure('Name','Normalized SNR vs. I ratio')
    semilogx(GYratio,SNRsnorm(:,1),'gd')
    hold
    semilogx(YGratio,SNRsnorm(:,2),'ys')
end
%xlim([0.1 10])




    