function Fcyl= ForceOpt2 (x)

% Parameter
Z=19;
L_arm= 1000 +10*Z;
L_load=L_arm;
L_mc1= 0.4*L_arm;
L_mc2= 0.5*L_arm;
L_mc3= 0.45*L_arm;

h1= 180;
h2= 40;

a1=x(1);
b1= x(2);
a2= x(3);
b2= x(4);
a3= x(5);
b3= x(6);

H_liftmax=L_arm;
H_liftmin=400+5*Z;

m_load=70+Z;
m_1= 10;
m_2=0.005*2*L_arm;
m_3=1.1*m_2;

g=9.81;

W_load=m_load*g;
W_1=m_1*g;
W_2=m_2*g;
W_3=m_3*g;

%Calculating Phi (Phi_min and Phi_max)

Phi_min=asin((H_liftmin-h2-b1-h1)/(L_arm));
Phi_max=asin((H_liftmax-h2-b1-h1)/(L_arm));

%Calculating Epsilon (Eps_min and Eps_max)
Eps_min=asin((H_liftmin-h2-b1-h1)/(L_arm));
Eps_max=asin((H_liftmax-h2-b1-h1)/(L_arm));

% phi=linspace(Phi_min,Phi_max,5000);

%Calculating Beta
Cyl_Ymin= sqrt(b2^2+a2^2).*sin(Phi_min+atan(b2/a2))-b3;
Cyl_Xmin= a3-sqrt(b2^2+a2^2).*cos(Phi_min+atan(b2/a2));

Cyl_Ymax= sqrt(b2^2+a2^2).*sin(Phi_max+atan(b2/a2))-b3;
Cyl_Xmax= a3-sqrt(b2^2+a2^2).*cos(Phi_max+atan(b2/a2));

beta_min=atan(Cyl_Ymin./Cyl_Xmin);
beta_max=atan(Cyl_Ymax./Cyl_Xmax);

% Calculating Cylinder Force
c2=sqrt(a2^2+b2^2);
theta=atan(b2/a2);
eps=[Eps_min];
phi=[Phi_min];
beta=[beta_min];
N=length(beta);
F=zeros(5,N);

for i=1:N
M1= [0 1 0 1 0; ...
    1 0 1 0 0;...
    b1 a1 0 0 0;...
    L_arm*sin(eps(i)) -L_arm*cos(eps(i)) 0 0 0;...
    0 0 L_arm*sin(phi(i)) -L_arm*cos(phi(i)) cos(beta(i))*cos(phi(i))*b2+cos(beta(i))*sin(phi(i))*a2+sin(beta(i))*cos(phi(i))*a2-sin(beta(i))*sin(phi(i))*b2];

C1= [W_load+W_1;...
    0;...
    W_load*L_load+W_1*L_mc1;...
    W_2*cos(eps(i))*L_mc2; ...
    W_3*cos(phi(i))*L_mc3];

F(:,i)=M1\C1;
end

% RAx=F(1,:);
% RAy=W_2+F(2,:);
% RA=sqrt(RAx.^2+RAy.^2);
% 
% RBx=F(1,:);
% RBy=F(2,:);
% RB=sqrt(RBx.^2+RBy.^2);
% 
% RCx=F(5,:).*cos(beta)+F(3,:);
% RCy=W_3+F(4,:)-F(5,:).*sin(beta);
% RC=sqrt(RCx.^2+RCy.^2);
% 
% RDx=F(3,:);
% RDy=F(4,:);
% RD=sqrt(RDx.^2+RDy.^2);

% Fcyl=RA+RB+RC+RD+F(5,:);
Fcyl=F(5,:);

% Fcyl(1)=sqrt(Fcyl_Xmin.^2+Fcyl_Y.^2);
% Fcyl(2)=sqrt(Fcyl_Xmax.^2+Fcyl_Y.^2);
end