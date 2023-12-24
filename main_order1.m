%% 针对一阶倒立摆的线性化模型的稳定性分析（可控性分析）、极点配置实验
clear; clc; close all; warning off;

M = 2;
m = 0.1;
l =0.5;
I = 1/3*m*l^2;  % 注意转动惯量的计算
g = 9.8;

a23 = -m*m*g*l*l/(I*(m+M)+M*m*l*l);
a43 = m*g*l*(M+m)/(I*(m+M)+M*m*l*l);
b2 = (I+m*l*l)/(I*(m+M)+M*m*l*l);
b4 = -m*l/(I*(m+M)+M*m*l*l);

A = [0 1 0 0;
    0 0 a23 0;
    0 0 0 1;
    0 0 a43 0];
B = [0;
    b2;
    0;
    b4];
C =[1 0 0 0]; 

%% 能控性判断
S_c = [B A*B A^2*B A^3*B];
fprintf('rank(S_c) = %d\n', rank(S_c));

Q_o = [C;
    C*A;
    C*A^2;
    C*A^3]
fprintf('rank(Q_o) = %d\n', rank(Q_o));

%% 稳定性判断
% 劳斯判据
disp('eig(A)=');
eig(A)

% 李雅普诺夫稳定性判据
% P = lyap(A',eye(4))
% eig(P)

%% 状态反馈的极点配置问题
J = [-1 -2 -1+j*sqrt(3) -1-j*sqrt(3)];
K1 = acker(A,B,J); 
K2 = place(A,B,J); 
disp('K1 = ');
disp(K1);
disp('K2 = ');
disp(K2);
K = K2;

Ts = 0.05;
TT = 20;
tspan = 0:Ts:TT;
y0 = [0;0;pi/8;0]; % 初始状态：x, \dot{x}, \theta, \dot{theta}


%%  极点配置需要运行该部分

%[t,y] = ode45(@(t,y)pend_cart_1(y,I,m,M,l,g,-K*(y -[-1;0;0;0])),tspan,y0);  %运行
 
[t,y] = ode45(@(t,y)IP_order1_dynamic(y,I,m,M,l,g,0),tspan,y0);  %自由运动
% 


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
