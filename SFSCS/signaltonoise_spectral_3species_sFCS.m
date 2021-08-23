% Signal to noise FCCS
clear all
close all

normalize2single=0;

path= uigetdir; 
files_curveincl=dir([path '/*final_curveincls.txt']);
segnumbs=zeros(size(files_curveincl,1),3);
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
SNRs=zeros(size(files_SNR,1),3);
for i=1:size(files_SNR,1)
    [namedata,remain]=strtok(files_SNR(i).name,'.');
    namefile=files_SNR(i).name;
    ACFinfo=load([path '/' namefile]);
    %SNCh(i)=dataCFs(1,9);
    SNRs(i,1)=mean(ACFinfo(:,3));
    SNRs(i,2)=mean(ACFinfo(:,5));
    SNRs(i,3)=mean(ACFinfo(:,7));
end

figure('Name','SNRs')
plot(SNRs(:,1),'g')
hold on
plot(SNRs(:,2),'y')
hold on
plot(SNRs(:,3),'r')

files_Is=dir([path '/*final_fitparameters_seg.txt']);
Is=zeros(size(files_Is,1),3);
for i=1:size(files_Is,1)
    [namedata,remain]=strtok(files_Is(i).name,'.');
    namefile=files_Is(i).name;
    fitparameters=load([path '/' namefile]);
    %SNCh(i)=dataCFs(1,9);
    Is(i,1)=fitparameters(22,1);
    Is(i,2)=fitparameters(23,1);
    Is(i,3)=fitparameters(24,1);
end

GYratio=Is(:,1)./Is(:,2);
GRratio=Is(:,1)./Is(:,3);
YGratio=1./GYratio;
RGratio=1./GRratio;
RGYratio=Is(:,3)./(Is(:,1)+Is(:,2));
RYratio=Is(:,3)./Is(:,2);
YRratio=1./RYratio;

figure('Name','SNR vs. I ratio')
plot3(GYratio,RGYratio,SNRs(:,1),'gd')
hold on
plot3(GYratio,RGYratio,SNRs(:,2),'ys')
hold on
plot3(GYratio,RGYratio,SNRs(:,3),'ro')
zlabel('SNR')
xlabel('I(G)/I(Y)')
ylabel('I(Ch2)/(I(Y)+I(G))')
set(gca,'XScale','log')
set(gca,'YScale','log')

figure('Name','SNR vs. I ratio 2d')
pointsize=20;
subplot(1,3,1),scatter(GYratio,GRratio,pointsize,SNRs(:,1),'filled','d')
xlabel('I(G)/I(Y)')
ylabel('I(G)/I(Ch2)')
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([0.1 10]);
ylim([0.1 10]);
hold on
subplot(1,3,2),scatter(YGratio,YRratio,pointsize,SNRs(:,2),'filled','s')
xlabel('I(Y)/I(G)')
ylabel('I(Y)/I(Ch2)')
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([0.1 10]);
ylim([0.1 10]);
subplot(1,3,3),scatter(RGratio,RYratio,pointsize,SNRs(:,3),'filled','o')
%zlabel('SNR')
xlabel('I(Ch2)/I(G)')
ylabel('I(Ch2)/I(Y)')
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([0.1 10]);
ylim([0.1 10]);


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




    