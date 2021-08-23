%just for odd sized squared input, NaN for cropping
function AC__=manualCCC_complete(input_,input2_, input3_, Use)

sizeAC=9;
AC_=zeros(sizeAC,sizeAC,sizeAC,sizeAC,size(input_,3));

for framenum=1:size(input_,3)

AC=zeros(sizeAC,sizeAC,sizeAC,sizeAC);
ACk=zeros(sizeAC,sizeAC,sizeAC,sizeAC);
Use=double(Use);
input=double(input_(:,:,framenum));
input2=double(input2_(:,:,framenum));
input3=double(input3_(:,:,framenum));

Use(Use==0)=NaN;
input=input.*Use;
meaninput=nanmean(nanmean(input));
input=input-meaninput;
input2=input2.*Use;
meaninput2=nanmean(nanmean(input2));
input2=input2-meaninput2;
input3=input3.*Use;
meaninput3=nanmean(nanmean(input3));
input3=input3-meaninput3;

centerAC=ceil(sizeAC/2);
for ii=1:size(AC,1)
for jj=1:size(AC,2)
for kk=1:size(AC,3)
for mm=1:size(AC,4)
       
        shiftY=jj-centerAC;
        shiftX=ii-centerAC;
        shiftZ=kk-centerAC;
        shiftW=mm-centerAC;
k=0;
for i=1:size(input,1)
for j=1:size(input,2)
       
if i+shiftX>0 && i+shiftX<size(input,1)+1 && j+shiftY>0 && j+shiftY<size(input,2)+1 && i+shiftZ>0 && i+shiftZ<size(input,1)+1 && j+shiftW>0 && j+shiftW<size(input,2)+1
if  isnan(input(i,j))==0 && isnan(input2(i+shiftX,j+shiftY))==0  && isnan(input3(i+shiftZ,j+shiftW))==0  
%           if  (shiftX==0 && shiftY==0)==0 && (shiftZ==0 && shiftW==0)==0 && (shiftZ==shiftX && shiftW==shiftY)==0
        
product=input(i,j)*input2(i+shiftX,j+shiftY)*input3(i+shiftZ,j+shiftW);

%             
              
              
              
              
              
AC(ii,jj,kk,mm)=AC(ii,jj,kk,mm)+product;
k=k+1;

end
end



end
end
AC(ii,jj, kk, mm)=AC(ii,jj, kk, mm)/k;
ACk(ii,jj, kk, mm)=k;    

           
end
end
end
end

AC_(:,:,:,:,framenum)=AC./meaninput./meaninput2./meaninput3;

end
AC__=mean(AC_,5);