function draw_pend_2(y)
%x', theta1', theta2', x, theta1, theta2
x = y(4);
th1 = y(5);
th2 = y(6);

% kinematics
% x = 3;        % cart position
% th = 3*pi/2;   % pendulum angle
M=2; m1=0.5; m2=0.5; l1=0.5; l2=0.5; L=0.4; g=9.8;
I1 = 1/12*m1*(2*l1)^2; I2 = 1/12*m2*(2*l2)^2;

% dimensions
% L = 2;  % pendulum length
l=2;
W = 1*sqrt(M/5);  % cart width
H = .5*sqrt(M/5); % cart height
wr = .2; % wheel radius
%mr = .3*sqrt(m); % mass radius

% positions
% y = wr/2; % cart vertical position
y = wr/2+H/2; % cart vertical position
w1x = x-.9*W/2;
w1y = 0;
w2x = x+.9*W/2-wr;
w2y = 0;

px1 = x + l1*sin(th1);
py1 = y + l1*cos(th1);

px2 = x + l1*sin(th1)+ l2*sin(th2);
py2 = y + l1*cos(th1)+ l2*cos(th2);

plot([-10 10],[0 0],'k','LineWidth',2)
hold on
rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[.5 0.5 1],'LineWidth',1.5)
rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[0 0 0],'LineWidth',1.5)
rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[0 0 0],'LineWidth',1.5)

plot([x px1],[y py1],'color',[100, 100, 100]/255,'LineWidth',4)
plot([px1 px2],[py1 py2],'color',[180, 180, 180]/255,'LineWidth',4)

%rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[1 0.1 .1],'LineWidth',1.5)

 xlim([-4 4]);
 ylim([-4 4]);
 
% set(gca,'YTick',[-3 3])
%  set(gca,'XTick',[-5 5])

xlabel('$x(m)$','FontSize',30,'interpreter','latex')
ylabel('$y(m)$','FontSize',30,'interpreter','latex')
set(gca,'FontSize',30)

%set(gcf,'Position',[100 550 1000 400])
% box off
drawnow
hold off