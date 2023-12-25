%% 针对二阶倒立摆的线性化模型的稳定性分析（可控性分析）、极点配置实验

clear; clc; close all; warning off;
% syms M m_1 m_2 l_1 l_2 L g I_1 I_2    % 求线性化模型，使用符号求逆
% 具体数据求解
M=2; m_1=0.5; m_2=0.5; l_1=0.2; l_2=0.2; L=0.4; g=9.8;
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
 
%% 能控性、能观性判断
S_c = [B A*B A^2*B A^3*B A^4*B A^5*B];
fprintf('rank(S_c) = %d\n', rank(S_c));

% Q_o = [C; C*A; C*A^2; C*A^3; C*A^4; C*A^5];
% fprintf('rank(Q_o) = %d\n', rank(Q_o));

%% 劳斯判据
disp('eig(A)=');
eig(A)

%% 李雅普诺夫稳定性判据
% P = lyap(A',eye(6))
% eig(P)

%% 状态反馈极点配置
J = [-2+j*2 -2-j*2 -6 -7 -8 -9];
K1 = acker(A,B,J); 
K2 = place(A,B,J); 
disp('K1 = ');
disp(K1);
disp('K2 = ');
disp(K2);
K = K2;

Ts = 0.05;
TT = 15;
tspan = 0:Ts:TT;
y0 = [0;0; 0; 0; -pi/10; -pi/9];

%%  极点配置要运行该部分

[t,y] = ode45(@(t,y)IP_order2_dynamic(y,-K*(y - [0;0;0;-1;0;0])),tspan,y0);  %运行

%[t,y] = ode45(@(t,y)IP_order2_dynamic(y,0),tspan,y0);  %无控制输入的自由运动




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
