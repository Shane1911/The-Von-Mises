close 
clear
clc
N=50;
X=linspace(-2,2,N);
Z=linspace(1e-3,1,N);
k=6.9731;
% k=1;
mu=0.002;
for i=1:N
x=X(i);
for j=1:N
z=Z(j);
ymax1=@(q)sqrt(1-q.^2);
ymax2=@(q)-sqrt(1-q.^2);
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
% fun9=@(q,p) (1-q.^2-p.^2).^(0.5).*((1-2*0.3)*...
%     (-(x-q).*(p/k)./(((x-q).^2+(p/k).^2+z^2).^1.5.*(((x-q).^2+(p/k).^2+z^2).^0.5+z).^2)+...
%     -2*(x-q).*(p/k)./(((x-q).^2+(p/k).^2+z^2).*(((x-q).^2+(p/k).^2+z^2).^0.5+z).^3)+...
%     (p/k)./(((x-q).^2+(p/k).^2+z^2).^(0.5).*(((x-q).^2+(p/k).^2+z^2).^0.5+z).^2))+...
%     3*(x-q).^2.*(p/k)./((x-q).^2+(p/k).^2+z^2).^2.5);
% Txy=integral2(@(q,p)fun9(q,p),-1,1,ymax2,ymax1);
% T_xy(j,i)=1/pi/2/k*Txy;
% fun10=@(q,p) (1-q.^2-p.^2).^(0.5).*(x-q).*(p/k)./((x-q).^2+(p/k).^2+z^2).^2.5;
% Tyz=integral2(@(q,p)fun10(q,p),-1,1,ymax2,ymax1);
% T_yz(j,i)=3*z*mu/pi/2/k/k*Tyz;
fun11=@(q,p) (1-q.^2-p.^2).^(0.5).*(x-q).^2./((x-q).^2+(p/k).^2+z^2).^2.5;
Tzx=integral2(@(q,p)fun11(q,p),-1,1,ymax2,ymax1);
T_zx(j,i)=-3*z*mu/pi/2/k*Tzx;
end
end
VON=sqrt(0.5*((F_xx-F_zz).^2+(F_xx-F_yy).^2+(F_yy-F_zz).^2+6*F_xy.^2));
[x1,z1]=meshgrid(X,Z);
surf(x1,z1,VON,'LineStyle','none');
colormap('jet');
colorbar;
shading interp;
view(0,90);
[x1,z1]=meshgrid(X,Z);
contourf(x1,z1,VON);
colormap('jet');