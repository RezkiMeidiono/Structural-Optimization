%% Home Assignment Optimization
clear all

% Parameter
Z=19;
L_arm= 1000+10*Z;
L_load=L_arm;
L_mc1= 0.4*L_arm;
L_mc2= 0.5*L_arm;
L_mc3= 0.45*L_arm;

h1= 180;
h2= 40;

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
%============================================
GA=false;

if GA==true
    fun=@ForceOpt;
    nvars = 6;
    A   = [];
    Aeq = [];
    b   = [];
    beq = [];
    lb = [0,0,0,0,0,0];
    ub = [L_mc1,H_liftmin-h1-h2,L_mc3,H_liftmin-h1-h2,L_arm,H_liftmin-h1-h2];
    nonlcon=@nonlconfun;


    [x,fval,exitflag,output] = gamultiobj(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon);

    [m,k]=find(x(:,6)==max(x(:,6)));
    j=m(1);
    a1 = x(j,1);
    b1 = x(j,2);
    a2 = x(j,3);
    b2 = x(j,4);
    a3= x(j,5);
    b3 = x(j,6); 
else
    fun=@ForceOpt2;
    x0=[0 0 0 0 0 0];
    A   = [];
    Aeq = [];
    b   = [];
    beq = [];
    lb = [0,0,0,0,0,0];
    ub = [L_mc1,H_liftmin-h1-h2,L_mc3,H_liftmin-h1-h2,L_arm,H_liftmin-h1-h2];
    nonlcon=@nonlconfun;

    x=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
    a1 = x(1);
    b1 = x(2);
    a2 = x(3);
    b2 = x(4);
    a3= x(5);
    b3 = x(6); 
end
%Calculating Phi 


% Phi_min=asin((H_liftmin-h2-b1-h1)/(L_mc3*2));
% Phi_max=asin((H_liftmax-h2-b1-h1)/(L_mc3*2));
% phi=linspace(Phi_min,Phi_max,5000);


H_lift=linspace(H_liftmin,H_liftmax,100);

phi=asin((H_lift-h2-b1-h1)/(L_arm));
eps=asin((H_lift-h2-b1-h1)/(L_arm));



%Calculating Beta
Cyl_Y= sqrt(b2^2+a2^2).*sin(phi+atan(b2/a2))-b3;
Cyl_X= a3-sqrt(b2^2+a2^2).*cos(phi+atan(b2/a2));

beta=atan(Cyl_Y./Cyl_X);

figure(1)
plot(H_lift,phi*180/pi)
ylabel({'$\varphi$ (deg)'},'Interpreter','latex')
xlabel('Lift Height (mm)')

yyaxis right
plot(H_lift,beta*180/pi,"-r");
title('Arm and Actuator Angle vs. Lift Height')
ylabel({'$\beta$ (deg)'},'Interpreter','latex')
legend({'$\varphi$ (deg)','$\beta$ (deg)'},'Interpreter','latex')
grid on



%% Calculating Cylinder Force and Joints
c2=sqrt(a2^2+b2^2);
theta=atan(b2/a2);
N=length(H_lift);
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
RAx=F(1,:);
RAy=W_2+F(2,:);
RA=sqrt(RAx.^2+RAy.^2);

RBx=F(1,:);
RBy=F(2,:);
RB=sqrt(RBx.^2+RBy.^2);

RCx=F(5,:).*cos(beta)+F(3,:);
RCy=W_3+F(4,:)-F(5,:).*sin(beta);
RC=sqrt(RCx.^2+RCy.^2);

RDx=F(3,:);
RDy=F(4,:);
RD=sqrt(RDx.^2+RDy.^2);

Fcyl=F(5,:);

figure(2)
plot(H_lift,RA/2,'.-',H_lift,RB/2,'.-',H_lift,RC/2,'.-',H_lift,RD/2,'.-')
legend('RA','RB','RC','RD')
ylabel('F Joints (N)')
xlabel('Lift Height (mm)')
title('Joints Forces vs. Lift Height')
grid on

% figure(3)
% plot(H_lift,F_m3(1,:),H_lift,F_m3(2,:),H_lift,F_m3(3,:),H_lift,F_m3(4,:))
% legend('FCx','FCy','FDx','FDy')
% ylabel('F Joints (N)')
% xlabel('Lift Height (mm)')
% title('Joints Forces vs. Lift Height')
% grid on
%% Calculating Cylinder Stroke
L_cyl=sqrt(Cyl_X.^2+Cyl_Y.^2);

figure(3)
plot(H_lift,L_cyl)
ylabel('Actuator Length (mm)')
xlabel('Lift Height (mm)')
title('Actuator Length vs. Lift Height')
grid on

%% Verification
% 

figure(4)
plot(H_lift,Fcyl)
ylabel('F actuator (N)')
xlabel('Lift Height (mm)')
title('Actuator Force vs. Lift Height')
grid on

Wr_cyl=trapz(L_cyl,Fcyl);
Wr_table=(W_load+W_1)*(H_liftmax-H_liftmin)+W_2*L_mc2*(sin(eps(end))-sin(eps(1)))+W_3*L_mc3*(sin(phi(end))-sin(phi(1)));

%% Pareto Front

if GA==true
    figure(5)
    plot(fval(:,1),fval(:,2),'*')
    xlabel('F Actuator LiftMin (N)')
    ylabel('F Actuator LiftMax (N)')

    grid on
end
%% Interference Check
figure(6)
plot(H_lift,H_lift,H_lift,(c2*sin(phi+theta)+h1))
% legend('H lift','H Piston Pin','H Piston Base Pin')
yline(h1+b3,'-g')
legend('H lift','H Piston Pin','H Piston Base Pin')
xlabel('Height (mm)')
ylabel('Height (mm)')
xlim([H_liftmin,1200])
ylim([100,H_liftmax])
grid on 