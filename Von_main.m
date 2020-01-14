%% SECTION TITLE:Von Mises stress calculation
% # shane:2020/01/03
close all
clear
clc
%% load zone
Q=[100 100];
E1=2.075e11;
E=E1/(1-0.3^2);
k=[6.9731 8.2435];
deta1=[1.0269 1.0174];
R1=1e-3*[3.9320 3.0557];
a=(6*k.^2.*deta1.*Q.*R1/(pi*E)).^(1/3);
b=a./k;
x=linspace(-2*a(1),2*a(1),500);
y_o=linspace(0,1*a(1),500);
y=y_o(50:500);
[x1,y1]=meshgrid(x,y);
% P=Q(1);
% R=abs(x1);e=sqrt(x1.^2+y1.^2);
% F_xx=P/2/pi*((1-y1./e).*(x1.^2./R.^2).*(1-2*0.3)./R.^2-3*y1.*x1.^2./e.^5);
% F_yy=P/2/pi*((1-y1./e).*(-x1.^2./R.^2)+y1.*x1.^2./e.^3).*(1-2*0.3)./R.^2;
% F_zz=-3*P*y1.^3./e.^5/2/pi;
% F_xz=-3*P*x1.*y1.^2./e.^5/2/pi;
% z=sqrt(0.5*((F_xx-F_yy).^2+(F_zz-F_yy).^2+(F_xx-F_zz).^2+6*F_xz.^2));
% contourf(x,y,z);
% colormap('jet');
% colorbar;
% surf(x,y,z,'LineStyle','none');
% xlim([min(x(:)) max(x(:))]);
% ylim([min(y(:)) max(y(:))]);
% colormap('jet');
% colorbar;
% shading interp;
% view(0,90);
F_xx=0;F_yy=0;F_zz=0;F_xz=0;
for yita=linspace(-a(1),a(1),50)
R=abs(x1-yita);e=sqrt((x1-yita).^2+y1.^2);
PH=Q(1);%3*Q(1)/(2*pi*a(1)*b(1));
F1=PH/2/pi*real(sqrt(1-(yita/a(1)).^2))*((1-2*0.3)./R.^2.*(1-y1./e)...
    .*(x1-yita).^2./R.^2-3*y1.*(x1-yita).^2./e.^5);
F2=PH/2/pi*real(sqrt(1-(yita/a(1)).^2))*((1-2*0.3)./R.^2.*(1-y1./e)...
    .*-(x1-yita).^2./R.^2+y1.*(x1-yita).^2./e.^3);
F3=-3*PH/2/pi*real(sqrt(1-(yita/a(1)).^2)).*y1.^3./e.^5;
F4=-3*PH/2/pi*real(sqrt(1-(yita/a(1)).^2)).*(x1-yita).*y1.^2./e.^5;
F_xz=F_xz+F4;
F_zz=F_zz+F3;
F_yy=F_yy+F2;
F_xx=F_xx+F1;
end
T_xx=0;T_yy=0;T_zz=0;T_xz=0;
for yita=linspace(-a(1),a(1),50)
R=abs(x1-yita);e=sqrt((x1-yita).^2+y1.^2);
mu=0.09;
PH=mu*Q(1);%/(2*pi*a(1)*b(1));
F1=PH/2/pi*real(sqrt(1-(yita/a(1)).^2))*((1-2*0.3)*((x1-yita)./e.^3-3*(x1-yita)./...
    (e.*(e+y1).^2)+(x1-yita).^3./(e.^3.*(e+y1).^2)+2*(x1-yita).^3./...
    (e.^2.*(e+y1).^3))-3*(x1-yita).^3./e.^5);
F2=PH/2/pi*real(sqrt(1-(yita/a(1)).^2))*((1-2*0.3)*((x1-yita)./e.^3-(x1-yita)./(e.*(e+y1).^2)));
F3=-3*PH/2/pi*real(sqrt(1-(yita/a(1)).^2)).*3*(x1-yita).*y1.^2./e.^5;
F4=-3*PH/2/pi*real(sqrt(1-(yita/a(1)).^2)).*3*(x1-yita).^2.*y1./e.^5;
T_xz=T_xz+F4;
T_zz=T_zz+F3;
T_yy=T_yy+F2;
T_xx=T_xx+F1;
end
F_x=F_xx+T_xx;
F_y=F_yy+T_yy;
F_z=F_zz+T_zz;
Fxz=F_xz+T_xz;
z=sqrt(0.5*((F_x-F_y).^2+(F_z-F_y).^2+(F_x-F_z).^2+6*Fxz.^2));
surf(x,y,z,'LineStyle','none');
xlim([min(x(:)) max(x(:))]);
ylim([min(y(:)) max(y(:))]);
colormap('jet');
colorbar;
shading interp;
view(0,90);
% contourf(x,y,z);
% colormap('jet');
% colorbar;









