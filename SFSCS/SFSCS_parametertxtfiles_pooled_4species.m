clear all
close all

% Specify directory in which three-species SFSFCS results files are located
path= uigetdir; 

% Specify directory and filename in which pooled three-species SFSFCS
% results shall be saved
path2= uigetdir; % directory
filename=sprintf('/2020-03-31_mp-Y-Ch2-G-A.txt'); % filename

files=dir([path '/*final_fitparameters_seg.txt']);
Ints=zeros(size(files,1),4);
bleaching_fractions=Ints;
taus=zeros(size(files,1),10);
Ns=taus;
for i=1:size(files,1)
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    parametersi=load([path '/' namefile]);
    for j=1:10
    Ns(i,j)=parametersi((j-1)*3+1,1);
    taus(i,j)=parametersi((j-1)*3+2,1);
    end
    Ints(i,:)=parametersi(35:38,1)';
    bleaching_fractions(i,:)=parametersi(31:34,1)';
end
Bs=Ints./Ns(:,1:4); % Brightness values
relcc12s=max([Ns(:,1)./Ns(:,5) Ns(:,2)./Ns(:,5)] ,[],2); % rel.cc 1-2
relcc13s=max([Ns(:,1)./Ns(:,6) Ns(:,3)./Ns(:,6)] ,[],2); % rel.cc 1-3
relcc14s=max([Ns(:,1)./Ns(:,7) Ns(:,4)./Ns(:,7)] ,[],2); % rel.cc 1-4
relcc23s=max([Ns(:,2)./Ns(:,8) Ns(:,3)./Ns(:,8)] ,[],2); % rel.cc 2-3
relcc24s=max([Ns(:,2)./Ns(:,9) Ns(:,4)./Ns(:,9)] ,[],2); % rel.cc 2-4
relcc34s=max([Ns(:,3)./Ns(:,10) Ns(:,4)./Ns(:,10)] ,[],2); % rel.cc 3-4
output=[Ns relcc12s relcc13s relcc14s relcc23s relcc24s relcc34s taus Bs bleaching_fractions];

fid1=fopen([path2 filename],'a');
fprintf(fid1,'sample\t N1\t N2\t N3\t N4\t Ncc12\t Ncc13\t Ncc14\t Ncc23\t Ncc24\t Ncc34\t relcc12\t relcc13\t relcc14\t relcc23\t relcc24\t relcc34\t tau1\t tau2\t tau3\t tau4\t taucc12\t taucc13\t taucc14\t taucc23\t taucc24\t taucc34\t B1\t B2\t B3\t B4\t bleaching1\t bleaching2\t bleaching3\t bleaching4\n'); %alpha
for i=1:size(files,1)
    fprintf(fid1,files(i).name(1:end-24));
    fprintf(fid1,'\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output(i,:)');
end
fclose all;