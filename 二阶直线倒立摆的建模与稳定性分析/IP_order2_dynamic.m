function dx = pend_cart_2(x,u)
%x', theta1', theta2', x, theta1, theta2

M=2; m1=0.5; m2=0.5; l1=0.2; l2=0.2; L=0.4; g=9.8;
I1 = 1/12*m1*(2*l1)^2; I2 = 1/12*m2*(2*l2)^2;
dx_ = x(1); dth1_ = x(2); dth2_ = x(3);
x_ = x(4); th1_ = x(5); th2_ = x(6); 
M11 = M + m1 + m2;
M12 = (m1*l1+m2*L)*cos(th1_);
M13 = m2*l2*cos(th2_);
M21 = M12;
M22 = I1 + m1*l1^2 + m2*L^2;
M23 = m2*L*l2*cos(th2_ - th1_);
M31 = M13;
M32 = M23;
M33 = I2 + m2*l2^2;
C11 = 0; C12 = -(m1*l1+m2*L)*sin(th1_)*dth1_;
C13 = -m2*l2*sin(th2_)*dth2_;
C21 = 0; C22 = 0; C23 = -m2*L*l2*dth2_*sin(th2_ - th1_);
C31 = 0; C32 = m2*L*l2*dth1_*sin(th2_ - th1_); C33 = 0;
G1 = 0; G2 = -(m1*l1+m2*L)*g*sin(th1_); 
G3 = -m2*g*l2*sin(th2_);
A = [M11 M12 M13 C11 C12 C13;
     M21 M22 M23 C21 C22 C23;
     M31 M32 M33 C31 C32 C33;
      0   0   0   1   0   0 ;
      0   0   0   0   1   0 ;
      0   0   0   0   0   1 ];
B = [u-G1;
     -G2;
     -G3;
     dx_;
     dth1_;
     dth2_];
sys = A\B;


dx(1,1)= sys(1);
dx(2,1) = sys(2);
dx(3,1) = sys(3);
dx(4,1) =  sys(4);
dx(5,1) =  sys(5);
dx(6,1) =  sys(6);

end

