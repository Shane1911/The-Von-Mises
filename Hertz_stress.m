%% Herte stress
% shane:2020/01/04
%% load zone
close
clear
clc
Q=[32.17 42.15 55.30 84.28];
E1=2.075e11;
E=E1/(1-0.3^2);
k=[6.9731 8.2435];
deta1=[1.0269 1.0174];
R1=1e-3*[3.9320 3.0557];
a=(6.*k(2).^2.*deta1(2).*Q.*R1(2)/(pi*E)).^(1/3);
b=a./k(2);
PH=3*Q'./(2*pi*a(2)*b(2));
% stress
t=0:0.1:2;
F_zz=-(1+k(1)*t.^2).^(-1);
F_xx=-((1+0.3)*(1-t.*atan(1./t)-0.5*(1+t.^2).^(-1)));
F_yy=F_xx;
Tao_1=0.5*abs(F_zz-F_xx);
figure(1)
plot(t,F_zz,'r-');
hold on
plot(t,F_xx);
hold on
plot(t,Tao_1);
% figure(2)
% von_stress = sqrt(0.5*((F_xx-F_yy).^2+(F_zz-F_yy).^2+(F_xx-F_zz).^2));
% plot(t,Tao_1);
% i=1;
% for t=0.02:0.02:2
% fun=@(v)(1-v).^(0.5)./(v+t^2).^(2.5);
% FZ(i)=-1.5*t.^3*quad(@(v)fun(v),1e-8,1);
% i=i+1;
% end
% figure(2)
% t1=0.1:0.2:2;
% plot(FZ,'-*');
%% iint calculation
% i=1;
% for t=0.01:0.01:2
% ymax1=@(x)sqrt(1-x.^2);
% % ymax2=@(x)-sqrt(1-x.^2);
% fun2=@(x,y)(1-x.^2-y.^2).^(0.5)./(x.^2+(y/k(1)).^2+t^2).^(2.5);
% FZ(i)=-6*t.^3/pi/k(1)*integral2(@(x,y)fun2(x,y),1e-8,1,1e-8,ymax1);
% i=i+1;
% end
% figure(2)
% t1=0.1:0.2:2;
% plot(FZ,'-*');
% t=0.1
% ymax1=@(x)sqrt(1-x.^2);
% ymax2=@(x)-sqrt(1-x.^2);
% fun2=@(x,y)x.*(1-x.^2-y.^2).^(0.5)./((x.^2+y.^2+t^2).^(0.5).*((x.^2+y.^2+t^2).^(0.5)+t).^2);
% FZ=integral2(@(x,y)fun2(x,y),-1,1,ymax2,ymax1)