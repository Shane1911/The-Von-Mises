% von mises _ Hertz
close all
clear all
clc
N=50;
X=linspace(-2,2,N);
Z=linspace(1e-3,1,N);
k=6.9731;
mu=0.005625;
% k=5;
for i=1:N
x=X(i);
for j=1:N
z=Z(j);
ymax1=@(q)sqrt(1-q.^2);
ymax2=@(q)-sqrt(1-q.^2);
fun2=@(q,p)(1-q.^2-p.^2).^(0.5)./((x-q).^2+(p/k).^2+z^2).^(2.5);
FZ=integral2(@(q,p)fun2(q,p),-1,1,ymax2,ymax1);
F_zz(j,i)=-3*z^3/pi/2/k*FZ;
fun3=@(q,p) (1-q.^2-p.^2).^(0.5).*((1-2*0.3)./((x-q).^2+(p/k).^2).*...
    ((1-z./((x-q).^2+(p/k).^2+z^2).^(0.5)).*(((x-q).^2-(p/k).^2)./...
    ((x-q).^2+(p/k).^2))+z*(p/k).^2./((x-q).^2+(p/k).^2+z^2).^(1.5))...
    -3*z*(x-q).^2./((x-q).^2+(p/k).^2+z^2).^(2.5));
Fx=integral2(@(q,p)fun3(q,p),-1,1,ymax2,ymax1);
F_xx(j,i)=1/pi/2/k*Fx;
fun4=@(q,p) (1-q.^2-p.^2).^(0.5).*((1-2*0.3)./((x-q).^2+(p/k).^2).*...
    ((1-z./((x-q).^2+(p/k).^2+z^2).^(0.5)).*(((p/k).^2-(x-q).^2)./...
    ((x-q).^2+(p/k).^2))+z*(x-q).^2./((x-q).^2+(p/k).^2+z^2).^(1.5))...
    -3*z*(p/k).^2./((x-q).^2+(p/k).^2+z^2).^(2.5));
Fy=integral2(@(q,p)fun4(q,p),-1,1,ymax2,ymax1);
F_yy(j,i)=1/pi/2/k*Fy;
fun5=@(q,p) (1-q.^2-p.^2).^(0.5).*(x-q)./((x-q).^2+(p/k).^2+z^2).^(2.5);
Fxz=integral2(@(q,p)fun5(q,p),-1,1,ymax2,ymax1);
F_xz(j,i)=-3*z^2/pi/2/k*Fxz;
fun6=@(q,p) (1-q.^2-p.^2).^(0.5).*((1-2*0.3).*((x-q)./((x-q).^2+(p/k).^2+z^2).^1.5-...
    3*(x-q)./(((x-q).^2+(p/k).^2+z^2).^0.5.*(((x-q).^2+(p/k).^2+z^2).^0.5+z).^2)+...
    (x-q).^3./(((x-q).^2+(p/k).^2+z^2).^1.5.*(((x-q).^2+(p/k).^2+z^2).^0.5+z).^2)+...
    2*(x-q).^3./(((x-q).^2+(p/k).^2+z^2).*(((x-q).^2+(p/k).^2+z^2).^0.5+z).^3))-...
    3*(x-q).^3./((x-q).^2+(p/k).^2+z^2).^2.5);
Txx=integral2(@(q,p)fun6(q,p),-1,1,ymax2,ymax1);
T_xx(j,i)=mu/pi/2/k*Txx;
fun7=@(q,p) (1-q.^2-p.^2).^(0.5).*((1-2*0.3).*((x-q)./((x-q).^2+(p/k).^2+z^2).^1.5-...
    (x-q)./(((x-q).^2+(p/k).^2+z^2).^0.5.*(((x-q).^2+(p/k).^2+z^2).^0.5+z).^2)+...
    (x-q).*(p/k).^2./(((x-q).^2+(p/k).^2+z^2).^1.5.*(((x-q).^2+(p/k).^2+z^2).^0.5+z).^2)+...
    2*(x-q).*(p/k).^2./(((x-q).^2+(p/k).^2+z^2).*(((x-q).^2+(p/k).^2+z^2).^0.5+z).^3))-...
    3*(x-q).*(p/k).^2./((x-q).^2+(p/k).^2+z^2).^2.5);
Tyy=integral2(@(q,p)fun7(q,p),-1,1,ymax2,ymax1);
T_yy(j,i)=mu/pi/2/k*Tyy;
fun8=@(q,p) (1-q.^2-p.^2).^(0.5).*(x-q)./((x-q).^2+(p/k).^2+z^2).^2.5;
Tzz=integral2(@(q,p)fun8(q,p),-1,1,ymax2,ymax1);
T_zz(j,i)=-3*z^2*mu/pi/2/k*Tzz;
fun11=@(q,p) (1-q.^2-p.^2).^(0.5).*(x-q).^2./((x-q).^2+(p/k).^2+z^2).^2.5;
Tzx=integral2(@(q,p)fun11(q,p),-1,1,ymax2,ymax1);
T_zx(j,i)=-3*z*mu/pi/2/k*Tzx;
end
end
FX=F_xx+T_xx;
FY=F_yy+T_yy;
FZ=F_zz+T_zz;
FXZ=F_xz+T_zx;
% VON=sqrt(0.5*(2*(F_xx-F_zz).^2));
VON=1.3102*sqrt(0.5*((FX-FZ).^2+(FX-FY).^2+(FY-FZ).^2+6*FXZ.^2));
[x1,z1]=meshgrid(X,Z);
surf(x1,z1,VON,'LineStyle','none');
colormap('jet');
colorbar;
shading interp;
view(0,90);
[x1,z1]=meshgrid(X,Z);
contourf(x1,z1,VON);
colormap('jet');
colorbar;
% fun5=@(q,p) (1-q.^2-p.^2).^(0.5).*(x-q)./((x-q).^2+(p/k).^2+z^2).^(2.5);