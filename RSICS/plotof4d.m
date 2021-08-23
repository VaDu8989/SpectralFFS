function out=plotof4d(in)

center=(size(in,1)+1)/2;
n=1;
out=zeros(center^4,3);
for i=1:size(in,1)
    for j=1:size(in,1)
        for k=1:size(in,1)
            for m=1:size(in,1)
                
                out(n,1)=((i-center)^2+(k-center)^2)^0.5;
                out(n,2)=((j-center)^2+(m-center)^2)^0.5;
                out(n,3)=in(i,j,k,m);
                n=n+1;
            end
        end
    end
end

