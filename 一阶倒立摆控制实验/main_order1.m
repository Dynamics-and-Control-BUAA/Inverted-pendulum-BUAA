clear all;
close all;
clc;

m = 0.2;
M = 0.6;
l = 0.5;
g = 9.8;
I = m*l^2/3;
a23 = -m*m*g*l*l/(I*(m+M)+M*m*l*l);
a43 = m*g*l*(M+m)/(I*(m+M)+M*m*l*l);
b2 = (I+m*l*l)/(I*(m+M)+M*m*l*l);
b4 = -m*l/(I*(m+M)+M*m*l*l);

%linearized system

A = [0 1 0 0;
    0 0 a23 0;
    0 0 0 1;
    0 0 a43 0];
B = [0;
    b2;
    0;
    b4]; 

C = [1  0  0  0
     0  0  1  0];

eig(A);
rank(ctrb(A,B));

Ts = 0.05;
TT = 10;
tspan = 0:Ts:TT;
y0 = [0;0;pi/8;0]; % 初始状态：x, \dot{x}, \theta, \dot{theta}


%% 能控能观性判断
S_c = [B A*B A^2*B A^3*B];
fprintf('rank(S_c) = %d\n', rank(S_c))


%% 极点配置
J = [-1 -2 -1+j*sqrt(3) -1-j*sqrt(3)];
K = place(A,B,J);



%% LQR
% Q = [1 0 0 0; 
%     0 1 0 0;
%     0 0 1 0;
%     0 0 0 1]
% R = 0.001;
% K_lqr = lqr(A,B,Q,R)
% K=K_lqr;
%eig(A-B*K)



%%  极点配置和LQR需要运行该部分

% [t,y] = ode45(@(t,y)IP_order1_dynamic(y,I,m,M,l,g,-K*(y -[-1;0;0;0])),tspan,y0);  %运行
% 
% y1 = y - kron(ones(size(y,1),1),[-1;0;0;0]'); %查看控制输入
% y2 = y1';
% y3 = y2(:);
% size(y3);
% u = -kron(eye(size(y,1)),K)*y3;

%[t,y] = ode45(@(t,y)IP_order1_dynamic(y,I,m,M,l,g,0),tspan,y0);
%%无控制输入的自由运动
% 
% 


%% mpc
A = eye(4)+A*Ts;
B = B*Ts;
% 初始状态 [x; dx; theta; dtheta]
ref = [-1;0;0;0];
N = 20;
refs = repmat(ref,N,1);
% 保存数据
xs     = []; % x
thetas = []; % theta
ts     = []; % time
fs     = []; % F
y=[];
% 初始状态
x = y0;
t = 0;
% Q矩阵,x和theta权重稍大一点
Q = [100 0 0 0
     0 1 0 0
     0 0 100 0
     0 0 0 1];
R = 0.1;
low = -100;
hi = 100;
% 仿真循环
for i = 1:length(tspan)
    % mpc求解
    z = SolveLinearMPC(A, B, x*0, Q, R, low, hi, x, refs, N);
    u = z(1);
    x = A*x+B*u;
    % 保存数据
    y = [y,x];
    fs = [fs, u];
    thetas = [thetas; x(3)];
    xs = [xs, x(1)];
    ts = [ts; t];
    t = t + Ts;
end
y=y';
u= fs;


%% 绘图
gif_flag = 0; % 是否保存gif图
if gif_flag == 1
    %filename = 'gif_name_free_1.gif'; % 动画文件的文件名
    %filename = 'gif_name_feedback_1.gif'; % 动画文件的文件名
    %filename = 'gif_name_lqr_1.gif'; % 动画文件的文件名
    filename = 'gif_name_mpc_1.gif'; % 动画文件的文件名
end


figure('color',[1,1,1]);
set(gcf,'unit','centimeter','position',[2,2,40,23])
%grid on;
hold on;
plot(tspan,u,'LineWidth',3);
xlim([0,TT])         
xlabel('t(s)','FontSize',30)
ylabel('input','FontSize',30)
set(gca,'FontSize',30)
legend_hu=legend({'u_{mpc}'},'interpreter','latex');
set(legend_hu,'Orientation','horizon')  %legend横排
set(legend_hu,'Box','off');             %不显示方框  
hold off;

figure('color',[1,1,1]);
set(gcf,'unit','centimeter','position',[2,2,40,23])
%grid on;
hold on;
plot(tspan,y(:,1),'LineWidth',3);
plot(tspan,y(:,2),'LineWidth',3);
plot(tspan,y(:,3),'LineWidth',3);
plot(tspan,y(:,4),'LineWidth',3);
xlim([0,TT])         
xlabel('t(s)','FontSize',30)
ylabel('x','FontSize',30)
set(gca,'FontSize',30)
legend_hc=legend({'$x$','$$ \dot{x}$$','$$\theta$$','$$\dot{\theta}$$'},'interpreter','latex');
set(legend_hc,'Orientation','horizon')  %legend横排
set(legend_hc,'Box','off');             %不显示方框  
hold off;


figure('color',[1,1,1]);
set(gcf,'unit','centimeter','position',[5,5,25,20])
for k=1:length(tspan)
    IP_order1_draw(y(k,:),m,M,l);
    if gif_flag == 1
         % 保存每一帧为 gif 图像
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if k == 1
            imwrite(imind,cm,filename,'gif','DelayTime',0.01,'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',0.01,'WriteMode','append');
        end
    end
end



%% pid控制
% [G H]=c2d(A,B,Ts);
% u=0;
% sumy0=[0 0 0 0]';
% kp1=1500;
% ki1=3;
% kd1=100;
% kp2=8000;
% ki2=1;
% kd2=150;
% for i=1:size(tspan,2)
%     lasty0=y0;
%     y0=G*y0+H*u;
%     Dy0=y0-lasty0;
%     sumy0=sumy0+y0;
%     u=kp1*y0(1)+ki1*sumy0(1)+kd1*Dy0(1)+kp2*y0(3)+ki2*sumy0(3)+kd2*Dy0(3);
%     y3(:,i)=y0;
% end
% y3=[y3(1,:);y3(3,:)];
% plot(t,y3);
% xlabel('t');
% ylabel('x&θ');
% title('PID');
% legend('x','θ');
