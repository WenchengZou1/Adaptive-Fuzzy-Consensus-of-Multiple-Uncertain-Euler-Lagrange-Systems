clear
clc
tic
T_max = 40; % 时长
T_stepsize = 0.0001; % 仿真步长
Ts=4;
q=[ 3 4
    3 5
    4 4
    1 2
    5 2
    6 0
   ];
[N,p]=size(q);

A=[0 1 0 1 0 0
    1 0 1 0 1 0
    0 1 0 0 0 1
    1 0 0 0 1 0
    0 1 0 1 0 1
    0 0 1 0 1 0];
B=[0 1 0 0 0 0];
D=diag([2,3,2,2,3,2]);
dq=zeros(N,p);


z=q;
count=1;
z1plt=zeros(N,(T_max-0)/T_stepsize+1);
z2plt=zeros(N,(T_max-0)/T_stepsize+1);
q1plt=zeros(N,(T_max-0)/T_stepsize+1);
q2plt=zeros(N,(T_max-0)/T_stepsize+1);
e11plt=zeros(N,(T_max-0)/T_stepsize+1);
e12plt=zeros(N,(T_max-0)/T_stepsize+1);
e21plt=zeros(N,(T_max-0)/T_stepsize+1);
e22plt=zeros(N,(T_max-0)/T_stepsize+1);
tauplt1=zeros(N,(T_max-0)/T_stepsize+1);
tauplt2=zeros(N,(T_max-0)/T_stepsize+1);
k=1;
k1=0.3;
beta=0.7;
r=0.1;
beta2=1.3;
W1=zeros(N,p);
W2=zeros(N,p);
e1=zeros(N,p);
e2=zeros(N,p);
dqr=zeros(N,p);
dthetah=zeros(N,1);
thetah=zeros(N,1);
tau=zeros(N,p);
ddq=zeros(N,p);
alpha1=1;
alpha2=1;
ktheta1=1;
ktheta2=1;
sigma1=1;
sigma2=1;
a=2.5;
for t=0:T_stepsize:T_max % 时间
    if ~mod(t*10000,10000) 
         disp("t="+t+"s/"+T_max+"s");
    end
    if ~mod(t*10000,Ts*10000) %每隔Ts秒更新一次虚拟系统
        for i=1:N
            sumz=zeros(1,p);
            for j=1:N
                if A(i,j)
                    sumz=sumz+q(i,:)-q(j,:);
                end
            end
            z(i,:)=z(i,:)-k1*sumz;  %更新律
        end
    end
    for i=1:N
        if i==3 && t>0.029 
            aaaa=1;
        end
        e1(i,:)=q(i,:)-z(i,:);  % 计算位置跟踪误差
        for g=1:p %每一维单独考虑
            if norm(e1(i,g))>r
                W1(i,g)=sign(e1(i,g))*abs(e1(i,g))^(2*beta-1);
            else
                W1(i,g)=r^(2*beta-2)*(2-beta)*e1(i,g)+r^(2*beta-4)*(beta-1)*e1(i,g)^3;
            end
        end
        dqr(i,:)=-alpha1*W1(i,:)-alpha2*sign(e1(i,:)).*abs(e1(i,:)).^(2*beta2-1);  %相当于反步法的虚拟控制器
        e2(i,:)=dq(i,:)-dqr(i,:);  %计算速度跟踪误差
        for g=1:p
            if norm(e2(i,g))>r
                W2(i,g)=sign(e2(i,g))*abs(e2(i,g))^(2*beta-1);
            else
                W2(i,g)=(r^(2*beta-2))*(2-beta)*e2(i,g)+(r^(2*beta-4))*(beta-1)*e2(i,g)^3;
            end
        end
        dthetah(i)=0.5*S(e1(i,:),e2(i,:))'*S(e1(i,:),e2(i,:))*e2(i,:)*e2(i,:)'-ktheta1*thetah(i)-ktheta2*sign(thetah(i))*abs(thetah(i)).^(2*beta2-1); % 模糊逻辑的梯度计算
        tau(i,:)=-sigma1*W2(i,:)-sigma2*sign(e2(i,:)).*abs(e2(i,:)).^(2*beta2-1)-1/(2*a^2)*thetah(i)*S(e1(i,:),e2(i,:))'*S(e1(i,:),e2(i,:))*e2(i,:)-0.5*e2(i,:)-e1(i,:);
%        tau(i,:)=-0.5*e2(i,:)-e1(i,:);  % 控制器
         %以下为模型

        
    end
    for i=1:N
        g=[0.2*sin(q(i,1)),1*cos(q(i,2))*sin(q(i,2))];
        C=[-2*0.04*dq(i,2)*sin(q(i,2)) -0.04*dq(i,2)*sin(q(i,2))
           -0.04*dq(i,2)*sin(q(i,2)) 0 ];
        M=[2.25+0.1+4*0.04*cos(q(i,2)) 0.1+2*0.04*cos(q(i,2))
           0.1+2*0.04*cos(q(i,2)) 0.1+2.25 ];
        d=[0.2*cos(t)*sin(0.3*t)-0.1*sin(0.25*t),0.2*sin(t)*cos(0.3*t)-0.1*cos(0.25*t)];
        ddq(i,:)=(M\(tau(i,:)-(C*dq(i,:)')'-g+d)')';
%         ddq(i,:)=tau(i,:);
        dq(i,:)=dq(i,:)+T_stepsize*ddq(i,:);
        q(i,:)=q(i,:)+T_stepsize*dq(i,:);
        thetah(i)=thetah(i)+T_stepsize*dthetah(i);
    end
    
    q1plt(:,count)=q(:,1);
    q2plt(:,count)=q(:,2);
    z1plt(:,count)=z(:,1);
    z2plt(:,count)=z(:,2);
    e11plt(:,count)=e1(:,1);
    e12plt(:,count)=e1(:,2);
    e21plt(:,count)=e2(:,1);
    e22plt(:,count)=e2(:,2);
    tauplt1(:,count)=tau(:,1);
    tauplt2(:,count)=tau(:,2);
    count=count+1;
end

save bishedata

linewd=2;
ptzt=14;
latexzt=16;
figure(100)
subplot(2,1,1)
plot(0:T_stepsize:T_max,z1plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(4,:),"LineWidth",linewd,"Color","#77AC30","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$z_{x}$",'Interpreter','latex','Fontsize',latexzt);
% legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)
subplot(2,1,2)
plot(0:T_stepsize:T_max,z2plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(4,:),"LineWidth",linewd,"Color","#AC7730","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$z_{y}$",'Interpreter','latex','Fontsize',latexzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt-1,"NumColumns",3,"Box","on")


figure(1)
plot(0:T_stepsize:T_max,z1plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(4,:),"LineWidth",linewd,"Color","#77AC30","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z1plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$z_{x}$",'Interpreter','latex','Fontsize',latexzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)



figure(2)
plot(0:T_stepsize:T_max,z2plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(4,:),"LineWidth",linewd,"Color","#AC7730","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,z2plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$z_{y}$",'Interpreter','latex','Fontsize',latexzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)


figure(3)
plot(0:T_stepsize:T_max,q1plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q1plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q1plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q1plt(4,:),"LineWidth",linewd,"Color","#77AC30","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q1plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q1plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$q_{x}$",'Interpreter','latex','Fontsize',latexzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)


figure(4)
plot(0:T_stepsize:T_max,q2plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q2plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q2plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q2plt(4,:),"LineWidth",linewd,"Color","#77AC30","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q2plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,q2plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$q_{y}$",'Interpreter','latex','Fontsize',latexzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)


figure(5)
plot(0:T_stepsize:T_max,e11plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e11plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e11plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e11plt(4,:),"LineWidth",linewd,"Color","#77AC30","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e11plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e11plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$e_{1x}$",'Interpreter','latex','Fontsize',latexzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)


figure(6)
plot(0:T_stepsize:T_max,e12plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e12plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e12plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e12plt(4,:),"LineWidth",linewd,"Color","#77AC30","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e12plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e12plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$e_{1y}$",'Interpreter','latex','Fontsize',latexzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)

figure(7)
plot(0:T_stepsize:T_max,e21plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e21plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e21plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e21plt(4,:),"LineWidth",linewd,"Color","#77AC30","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e21plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e21plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$e_{2x}$",'Interpreter','latex','Fontsize',latexzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)


figure(8)
plot(0:T_stepsize:T_max,e22plt(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e22plt(2,:),"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e22plt(3,:),"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e22plt(4,:),"LineWidth",linewd,"Color","#77AC30","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e22plt(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,e22plt(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("time(s)",'Fontsize',ptzt);
ylabel("$e_{2y}$",'Interpreter','latex','Fontsize',latexzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)




figure(9)
plot3(q1plt(1,:),q2plt(1,:),0:T_stepsize:T_max,"LineWidth",linewd,"Color","#D95319","LineStyle","-");
hold on;
plot3(q1plt(2,:),q2plt(2,:),0:T_stepsize:T_max,"LineWidth",linewd,"Color","r","LineStyle","-");
hold on;
plot3(q1plt(3,:),q2plt(3,:),0:T_stepsize:T_max,"LineWidth",linewd,"Color","#0072BD","LineStyle","-");
hold on;
plot3(q1plt(4,:),q2plt(4,:),0:T_stepsize:T_max,"LineWidth",linewd,"Color","#77AC30","LineStyle","-");
hold on;
plot3(q1plt(5,:),q2plt(5,:),0:T_stepsize:T_max,"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot3(q1plt(6,:),q2plt(6,:),0:T_stepsize:T_max,"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
hold on;
grid on;
xlabel("x",'Fontsize',ptzt);
ylabel("y",'Fontsize',ptzt);
zlabel("t",'Fontsize',ptzt);
legend('$i=1$','$i=2$','$i=3$','$i=4$','$i=5$','$i=6$','Interpreter','latex','Fontsize',ptzt)

figure(10)
plot(0:T_stepsize:T_max,tauplt1(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","--");
hold on;
plot(0:T_stepsize:T_max,tauplt1(2,:),"LineWidth",linewd,"Color","r","LineStyle","-.");
hold on;
plot(0:T_stepsize:T_max,tauplt1(3,:),"LineWidth",linewd+1,"Color","#0072BD","LineStyle",":");
hold on;
plot(0:T_stepsize:T_max,tauplt1(4,:),"LineWidth",linewd+0.5,"Color","#77AC30","LineStyle","--");
hold on;
plot(0:T_stepsize:T_max,tauplt1(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,tauplt1(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
set(gca,"Fontsize",14);
xlabel("Time(s)");
ylabel("$tau_{ix}$",'Interpreter','latex','Fontsize',21);
legend('$i=1$','$i=2$','$i=3$','$i=4$','Interpreter','latex','Fontsize',ptzt)
figure(11)
plot(0:T_stepsize:T_max,tauplt2(1,:),"LineWidth",linewd,"Color","#D95319","LineStyle","--");
hold on;
plot(0:T_stepsize:T_max,tauplt2(2,:),"LineWidth",linewd,"Color","r","LineStyle","-.");
hold on;
plot(0:T_stepsize:T_max,tauplt2(3,:),"LineWidth",linewd+1,"Color","#0072BD","LineStyle",":");
hold on;
plot(0:T_stepsize:T_max,tauplt2(4,:),"LineWidth",linewd+0.5,"Color","#77AC30","LineStyle","--");
hold on;
plot(0:T_stepsize:T_max,tauplt2(5,:),"LineWidth",linewd,"Color","#7E2F8E","LineStyle","-");
hold on;
plot(0:T_stepsize:T_max,tauplt2(6,:),"LineWidth",linewd,"Color","#4DBEEE","LineStyle","-");
set(gca,"Fontsize",14);
xlabel("Time(s)");
ylabel("$tau_{iy}$",'Interpreter','latex','Fontsize',21);
legend('$i=1$','$i=2$','$i=3$','$i=4$','Interpreter','latex','Fontsize',ptzt)

toc
save data1.mat
% load data1.mat

[-1 1 
 0 1 
 1 1
 -1 0
 0 0
 1 0]