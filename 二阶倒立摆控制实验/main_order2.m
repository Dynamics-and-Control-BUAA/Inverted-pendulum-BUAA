clear all;
close all;
clc;

%M=2; m_1=0.5; m_2=0.5; l_1=0.5; l_2=0.5; L=0.4; g=9.8;
M=0.6; m_1=0.2; m_2=0.2; l_1=0.5; l_2=0.5; L=0.4; g=9.8;
I_1 = 1/12*m_1*(2*l_1)^2; I_2 = 1/12*m_2*(2*l_2)^2;

M_11 = M+m_1+m_2; M_12 = m_1*l_1+m_2*L;  M_13 = m_2*l_2;
M_21 = M_12; M_22 = I_1+m_1*l_1*l_1+m_2*L*L;  M_23 = m_2*L*l_2;   
M_31 = M_13; M_32 = M_23; M_33 = I_2+m_2*l_2*l_2;
M = [M_11 M_12 M_13; M_21 M_22 M_23; M_31 M_32 M_33];
G = [0 0 0; 0 (m_1*l_1+m_2*L)*g 0; 0 0 m_2*g*l_2];
U = [1; 0; 0];
A_a = M\G;
B_b = M\U;

A = zeros(6,6);
B = zeros(6,1);
A(1:3, 4:end) = A_a;
A(4:end, 1:3) = eye(3);
B(1:3) = B_b;
C = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

S_c = [B A*B A^2*B A^3*B A^4*B A^5*B];
fprintf('rank(S_c) = %d\n', rank(S_c));

Ts = 0.05;
TT = 15;
tspan = 0:Ts:TT;
y0 = [0;0; 0; 0; -pi/10; -pi/9];  % 初始状态：\dot{x}, \dot{\theta_1}, \dot{theta_2}, x, \theta_1, \theta_2


%% 极点配置
J = [-6+j*2 -6-j*2 -6 -7 -8 -9];
K = place(A,B,J);

%% LQR
Q = eye(6)+diag([0,0,0,0,100,100]);
R = 0.001;
K_lqr = lqr(A,B,Q,R);
K =K_lqr;

eig(A-B*K)
%%


%%  极点配置和 LQR需要运行该部分

[t,y] = ode45(@(t,y)IP_order2_dynamic(y,-K*(y - [0;0;0;-1;0;0])),tspan,y0);  %运行

%[t,y] = ode45(@(t,y)IP_order2_dynamic(y,0),tspan,y0);


%% mpc
% A = eye(6)+A*Ts;
% B = B*Ts;
% % 初始状态 [x; dx; theta; dtheta]
% ref = [0;0;0;-1;0;0];
% N = 30;
% refs = repmat(ref,N,1);
% % 保存数据
% xs     = []; % x
% thetas = []; % theta
% ts     = []; % time
% fs     = []; % F
% y=[];
% % 初始状态
% x = y0;
% t = 0;
% % Q矩阵,x和theta权重稍大一点
% Q = eye(6)+diag([0,0,0,100,100,100]);
% R = 0.01;
% low = -100;
% hi = 100;
% % 仿真循环
% for i = 1:length(tspan)
%     % mpc求解
%     z = SolveLinearMPC(A, B, x*0, Q, R, low, hi, x, refs, N);
%     u = z(1);
%     x = A*x+B*u;
%     % 保存数据
%     y = [y,x];
%     fs = [fs, u];
%     thetas = [thetas; x(6)];
%     xs = [xs, x(4)];
%     ts = [ts; t];
%     t = t + Ts;
% end
% y=y';


%% 绘图
gif_flag = 0; % 是否保存gif图
if gif_flag == 1
    filename = 'gif_name_free_2.gif'; % 动画文件的文件名
    %filename = 'gif_name_feedback_2.gif'; % 动画文件的文件名
    %filename = 'gif_name_lqr_2.gif'; % 动画文件的文件名
    %filename = 'gif_name_mpc_2.gif'; % 动画文件的文件名
end

figure('color',[1,1,1]);
set(gcf,'unit','centimeter','position',[2,2,40,23])
%grid on;
hold on;
plot(tspan,y(:,1),'LineWidth',3);
plot(tspan,y(:,2),'LineWidth',3);
plot(tspan,y(:,3),'LineWidth',3);
plot(tspan,y(:,4),'LineWidth',3);
plot(tspan,y(:,5),'LineWidth',3);
plot(tspan,y(:,6),'LineWidth',3);
xlim([0,TT])         
xlabel('t(s)','FontSize',30)
ylabel('x','FontSize',30)
set(gca,'FontSize',30)
legend_y=legend({'$ \dot{x}$','$ \dot{\theta}_1$','$ \dot{\theta}_2$','$x$','$\theta_1$','$\theta_2$'},'interpreter','latex');
set(legend_y,'Orientation','horizon')  %legend横排
set(legend_y,'Box','off');             %不显示方框
hold off;


figure('color',[1,1,1]);
set(gcf,'unit','centimeter','position',[5,5,25,20])
for k=1:length(tspan)
    IP_order2_draw(y(k,:));
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



%% pid
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
