function dx = pend_cart_1(x,I,m,M,L,g,u)
%X1 = x ; cart position
%X2 = x_dot ; cart velocity
%X3 = theta; pend angular
%X4 = theta_dot; pend angular velocity

theta = x(3);
sx = sin(theta);
cx = cos(theta);
M_1 = [M+m, m*L*cx; m*L*cx, I+m*L^2];
G_1 = [0, -m*L*sx*x(4);0,0];

G_2 = M_1\G_1;
F_1 = M_1\[u;m*g*L*sx];


dx(1,1)= x(2);
dx(2,1) = G_2(1)+F_1(1);
dx(3,1) = x(4);
dx(4,1) =  G_2(2)+F_1(2);

