clear all
close all

% Specify directory in which three-species SFSFCS results files are located
path= uigetdir;

% Specify directory and filename in which pooled three-species SFSFCS
% results shall be saved
path2= uigetdir; % directory
filename=sprintf('/2020-03-31_mp-G_mp-Y_mp-Ch2.txt'); % filename

files=dir([path '/*final_fitparameters_seg.txt']);
Ints=zeros(size(files,1),3);
bleaching_fractions=Ints;
taus=zeros(size(files,1),6);
Ns=taus;

for i=1:size(files,1)
    [namedata,remain]=strtok(files(i).name,'.');
    namefile=files(i).name;
    parametersi=load([path '/' namefile]);
    for j=1:6
    Ns(i,j)=parametersi((j-1)*3+1,1);
    taus(i,j)=parametersi((j-1)*3+2,1);
    end
    Ints(i,:)=parametersi(22:24,1)';
    bleaching_fractions(i,:)=parametersi(19:21,1)';

end
Bs=Ints./Ns(:,1:3); % Brightness values
relcc12s=max([Ns(:,1)./Ns(:,4) Ns(:,2)./Ns(:,4)] ,[],2); % rel.cc 1-2
relcc13s=max([Ns(:,1)./Ns(:,5) Ns(:,5)./Ns(:,5)] ,[],2); % rel.cc 1-3
relcc23s=max([Ns(:,2)./Ns(:,6) Ns(:,3)./Ns(:,6)] ,[],2); % rel.cc 2-3
output=[Ns relcc12s relcc13s relcc23s taus Bs bleaching_fractions];

fid1=fopen([path2 filename],'a');
fprintf(fid1,'sample\t N1\t N2\t N3\t Ncc12\t Ncc13\t Ncc23\t relcc12\t relcc13\t relcc23\t tau1\t tau2\t tau3\t taucc12\t taucc13\t taucc23\t B1\t B2\t B3\t bleaching1\t bleaching2\t bleaching3\n'); %alpha
for i=1:size(files,1)
    fprintf(fid1,files(i).name(1:end-24));
    fprintf(fid1,'\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n',output(i,:)');
end
fclose all;