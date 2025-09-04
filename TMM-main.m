% transfer_matrix_method

clear all;
close all;

% 二相组分

a=4e-1;a1=a/2;a2=a-a1;
% rho1=10000;rho2=1500;
% rho1=800;rho2=3790;
rho1=1000;rho2=2000;
E1=1e9;E2=1e9;

% 泊松比
v1=0.3;
v2=0.3;

% lame常数+剪切模量
lambda1=E1*v1/(1+v1)/(1-2*v1);
lambda2=E2*v2/(1+v2)/(1-2*v2);
G1=E1/2/(1+v1);
G2=E2/2/(1+v2);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 纵波波速
cL1=sqrt((lambda1+2*G1)/rho1);
cL2=sqrt((lambda2+2*G2)/rho2);

f=linspace(0,1100,500);
f=f';
w=2*pi.*f;
N=length(w);

% 纵波波矢
kL1=w./cL1; 
kL2=w./cL2;

% 纵波矩阵
KL1=zeros(2,2,N);
HL1=zeros(2,2,N);
KL2=zeros(2,2,N);
HL2=zeros(2,2,N);
TL=zeros(2,2,N);

% EVL即是纵波Bloch矩阵特征值
EVL=zeros(2,N);

% 纵波的数值解法
for ii=1:N
    KL1(:,:,ii)=[exp(1i*kL2(ii)*a1),exp(-1i*kL2(ii)*a1);rho2*cL2^2*kL2(ii)*exp(1i*kL2(ii)*a1),-rho2*cL2^2*kL2(ii)*exp(-1i*kL2(ii)*a1)];
    HL1(:,:,ii)=[exp(1i*kL1(ii)*a1),exp(-1i*kL1(ii)*a1);rho1*cL1^2*kL1(ii)*exp(1i*kL1(ii)*a1),-rho1*cL1^2*kL1(ii)*exp(-1i*kL1(ii)*a1)];
    KL2(:,:,ii)=[exp(1i*kL2(ii)*a),exp(-1i*kL2(ii)*a);rho2*cL2^2*kL2(ii)*exp(1i*kL2(ii)*a),-rho2*cL2^2*kL2(ii)*exp(-1i*kL2(ii)*a)];
    HL2(:,:,ii)=[1,1;rho1*cL1^2*kL1(ii),-rho1*cL1^2*kL1(ii)];
    TL(:,:,ii)=inv(KL1(:,:,ii))*HL1(:,:,ii)*inv(HL2(:,:,ii))*KL2(:,:,ii);
    EVL(:,ii)=eig(TL(:,:,ii));
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% 横波波数
cT1=sqrt(G1/rho1);
cT2=sqrt(G2/rho2);

% 横波波矢
kT1=w./cT1; 
kT2=w./cT2;

% 横波矩阵
KT1=zeros(2,2,N);
HT1=zeros(2,2,N);
KT2=zeros(2,2,N);
HT2=zeros(2,2,N);
TT=zeros(2,2,N);

% EVT即是横波Bloch矩阵特征值
EVT=zeros(2,N);

% 横波的数值解法
for ii=1:N
    KT1(:,:,ii)=[exp(1i*kT2(ii)*a1),exp(-1i*kT2(ii)*a1);rho2*cT2^2*kT2(ii)*exp(1i*kT2(ii)*a1),-rho2*cT2^2*kT2(ii)*exp(-1i*kT2(ii)*a1)];
    HT1(:,:,ii)=[exp(1i*kT1(ii)*a1),exp(-1i*kT1(ii)*a1);rho1*cT1^2*kT1(ii)*exp(1i*kT1(ii)*a1),-rho1*cT1^2*kT1(ii)*exp(-1i*kT1(ii)*a1)];
    KT2(:,:,ii)=[exp(1i*kT2(ii)*a),exp(-1i*kT2(ii)*a);rho2*cT2^2*kT2(ii)*exp(1i*kT2(ii)*a),-rho2*cT2^2*kT2(ii)*exp(-1i*kT2(ii)*a)];
    HT2(:,:,ii)=[1,1;rho1*cT1^2*kT1(ii),-rho1*cT1^2*kT1(ii)];
    TT(:,:,ii)=inv(KT1(:,:,ii))*HT1(:,:,ii)*inv(HT2(:,:,ii))*KT2(:,:,ii);
    EVT(:,ii)=eig(TT(:,:,ii));
end

% kL是纵波Bloch波矢
REVL=real(EVL(1,:));
kL=(real(acos(REVL)/a))';

% 求解纵波模式解析解
F=rho1*cL1/rho2/cL2;
coska=cos(kL1.*a1).*cos(kL2.*a2)-(F+1/F)/2.*sin(kL1.*a1).*sin(kL2.*a2);
kk=real(acos(coska)/a);

% kT是横波Bloch波矢
REVT=real(EVT(1,:));
kT=(real(acos(REVT)/a))';

% 求解横波模式解析解
FF=rho1*cT1/rho2/cT2;
Scoska=cos(kT1.*a1).*cos(kT2.*a2)-(FF+1/FF)/2.*sin(kT1.*a1).*sin(kT2.*a2);
Skk=real(acos(Scoska)/a);

hold on
plot(kL*a,f,'-b')
plot(kT*a,f,'-r')
plot(kk*a,f,'o')
plot(Skk*a,f,'o')
hold off
