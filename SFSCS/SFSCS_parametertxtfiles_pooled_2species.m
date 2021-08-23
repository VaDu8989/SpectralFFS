clear all
close all

% Specify directory in which three-species SFSFCS results files are located
path= uigetdir; 

% Specify directory and filename in which pooled three-species SFSFCS
% results shall be saved
path2= uigetdir; % directory
filename=sprintf('/2020-03-30_mp-G_mp-Y.txt'); % filename

files=dir([path '/*final_fitparameters_seg.txt']);
Ints=zeros(size(files,1),2);
bleaching_fractions=Ints;
taus=zeros(size(files,1),3);
Ns=taus;
for i=1:size(files,1)
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    parametersi=load([path '/' namefile]);
    for j=1:3
    Ns(i,j)=parametersi((j-1)*3+1,1);
    taus(i,j)=parametersi((j-1)*3+2,1);
    end
    Ints(i,:)=parametersi(12:13,1)';
    bleaching_fractions(i,:)=parametersi(10:11,1)';
end
Bs=Ints./Ns(:,1:2); % Brightness
relccs=max([Ns(:,1)./Ns(:,3) Ns(:,2)./Ns(:,3)],[],2); % rel.cc
output=[Ns relccs taus Bs bleaching_fractions];

fid1=fopen([path2 filename],'a');
fprintf(fid1,'sample\t N1\t N2\t Ncc12\t relcc\t tau1 [s]\t tau2 [s]\t taucc12 [s]\t B1 [a.u.]\t B2 [a.u.]\t bleaching1\t bleaching2\n');
for i=1:size(files,1)
    fprintf(fid1,files(i).name(1:end-24));
    fprintf(fid1,'\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output(i,:)');
end
fclose all;