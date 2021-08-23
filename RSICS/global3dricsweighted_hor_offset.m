function out=global3dricsweighted_hor_offset(par,in)
w0=par(1);
N=par(2);
D=par(3);
offset=par(4);
ff=par(5);
x_offset=par(6);

x=abs(in(:,1)-x_offset);
y=in(:,2);
linet=in(:,3);

pix=in(:,4);
pixt=in(:,5);
weights=in(:,6);



out=(offset+(1/N).*exp(-(((x.*pix./w0).^2+(y.*pix./w0).^2)) ./(1+(4.*D.*(x.*pixt+y.*linet)./(w0^2)))) ./(1+4.*D.*(x.*pixt+y.*linet)./w0^2)./(1+4.*D.*(x.*pixt+y.*linet)./(ff*w0)^2).^(1/2)).*weights;


