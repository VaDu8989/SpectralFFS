%%% Performs auto/crosscorrelation (Hendrix et al.)


function ImageCor2=arbitrary_V0CC(input0, input1, Use)

ImageCor2=zeros(size(input0));

for i=1:size(input0,3)
Image=double(input0(:,:,i));
Image1=double(input1(:,:,i));
Size = [2*size(Image,1)-1, 2*size(Image,2)-1];
        
% Calculates normalization for non-selected regions
Norm=fft2(Use,Size(1),Size(2));
Norm=fftshift(ifft2(Norm.*conj(Norm)));
% Calculates fluctuation image
ImageFluct=Image-mean(Image(Use==1));
ImageFluct1=Image1-mean(Image1(Use==1));
% Applies selection and FFT
ImageFluct=fft2(ImageFluct.*Use,Size(1),Size(2));
ImageFluct1=fft2(ImageFluct1.*Use,Size(1),Size(2));
% Actual correlation
ImageCor = fftshift(real(ifft2(ImageFluct.*conj(ImageFluct1))));
% Corrects for shape of selected region
ImageCor = ImageCor./Norm;
ImageCor2(:,:,i) = ImageCor(ceil(Size(1)/4):round(Size(1)*3/4),ceil(Size(2)/4):round(Size(2)*3/4));
ImageCor2(:,:,i)=ImageCor2(:,:,i)./(mean2(Image(Use==1)).*mean2(Image1(Use==1)));
end
