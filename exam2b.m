function exam2b
clc; clear; 
p8
end

function p1
syms k_F k_M L real

P11 = [L/3; 2*L/3; 0];
P22 = [L/3; -L/3; 0];
P33 = [-2*L/3; -L/3; 0];

F11 = [0;0;k_F];
F22 = [0;0;k_F];
F33 = [0;0;k_F];

M11 = [0;0;k_M];
M22 = [0;0;-k_M];
M33 = [0;0;k_M];

t1 = skew(P11)*F11 + M11;
t2 = skew(P22)*F22 + M22;
t3 = skew(P33)*F33 + M33;

W = [t1 t2 t3; k_F k_F k_F];
end

function p2
sigmamax = 1000.00000000;
k_F = 0.00001000;
k_M = 0.00000010;
L = 0.80000000;
m = 1.80000000;
g = 9.81000000;
theta1 = -0.35000000;
theta2 = -0.93000000;
theta3 = 0.85000000;
u = [-0.03467122; 0.00000000; 0.06593665; 6.59366467];

R = Rz(theta1)*Ry(theta2)*Rx(theta3);

f = [0;0;m*g]-R*[0;0;u(4)];
tau = [u(1); u(2); u(3)];

mat2str(f)
mat2str(tau)
end

function p3
g = 9.8100000000;
m = 0.5644603441;
J = [0.0034633777 0.0000000000 0.0000000000; 0.0000000000 0.0031136313 0.0000000000; 0.0000000000 0.0000000000 0.0077965536];
o = [2.7932528414; 0.8421109857; -1.9864036471];
th = [-0.0674697331; -0.1283268077; -0.1145463442];
v = [-0.0002398116; 0.0546551743; -0.0713881971];
w = [-0.0222789339; -0.0608772850; -0.0797067668];
u = [-0.0000902619; 0.0151680096; 0.0100285101; 5.5637124596];

x0=[o;th;v;w];
xd = getRates(J,x0,u,m,g)

q= xd(1:3);
th= xd(4:6);
v= xd(7:9);
w = xd(10:12);

mat2str(q)
mat2str(th)
mat2str(v)
mat2str(w)
end

function p4
g = 9.8100000000;
m = 0.6313027451;
J = [0.0048203145 0.0000000000 0.0000000000; 0.0000000000 0.0048872037 0.0000000000; 0.0000000000 0.0000000000 0.0073289165];
o = [0.2517665986; -1.7679079152; 0.1821398749];
th = [0.0097877359; -0.0888319524; -0.0141941454];
v = [-0.0548472723; -0.0395326020; 0.0020399646];
w = [0.0149655780; -0.0060017666; 0.0401919313];
u = [0.0081970471; -0.0118145722; -0.0053849818; 6.7445788695];
t0 = 0.7000000000;
t1 = 1.2000000000;

x0=[o;th;v;w];
[t,x] = ode45(@(t,x) getRates(J,x,u,m,g),[t0 t1],x0);

x = x(end,:)'
end

function p5

g = 9.8100000000;
m = 0.7277349079;
J = [0.0033980172 0.0000000000 0.0000000000; 0.0000000000 0.0035324011 0.0000000000; 0.0000000000 0.0000000000 0.0074783112];
theta1 = 2.5939444800;
dt = 0.1600000000;

syms o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4 real
xc = [o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3]';
uc = [u1 u2 u3 u4]';
xe = [0;0;0;theta1;0;0;0;0;0;0;0;0];
ue = [0;0;0;m*g];
Ac = subs(subs(jacobian(getRates(J,xc,uc,m,g),xc),xc,xe),uc,ue);
Bc = subs(subs(jacobian(getRates(J,xc,uc,m,g),uc),xc,xe),uc,ue);



[Ad, Bd] = ssdata(c2d(ss(double(vpa(Ac)),double(vpa(Bc)),[],[]),dt))

mat2str(Ad)

mat2str(Bd)
end


function p8
syms p q r real
dp = q;
dq = -2*p +3*q + 8 *r + 8*cos(p);

x = [p;q];
u = r;
dx = [dp;dq];

xe = [5;0];
ue = vpa(solve(subs(dx,x,xe)==0,r));

Ac = double(subs(subs(jacobian(dx,x),x, xe),u,ue));
Bc = double(subs(subs(jacobian(dx,u),x, xe),u,ue));

mat2str(double(Ac));
mat2str(double(Bc));

dt = 0.1;

[Ad, Bd] = ssdata(c2d(ss(Ac,Bc,[],[]),dt));
mat2str(double(Ad));
mat2str(double(Bd));

K = dlqr(Ad,Bd,eye(2),eye(1));

x0 = [4.92; -0.18];
u0 = ue - K*(x0-xe);

mat2str(double(u0));

tf = 0.4;
x = x0;
t = 0;

while t<tf
    u = ue - K*(x-xe);
    xdot=matlabFunction(dx);
    [t,x] = ode45(@(t,x) xdot(x(1),x(2),double(u)),[t t+dt],x);
    x = x(end,:)';
    t = t(end,:)';
end
x
end

function xdot = getRates(J,x,u,m,g)
q= x(1:3);
th= x(4:6);
v= x(7:9);
w = x(10:12);

N = [Rx(th(3))'*Ry(th(2))'*[0;0;1] Rx(th(3))'*[0;1;0] eye(3)*[1;0;0]];
R = Rz(th(1))*Ry(th(2))*Rx(th(3));
f = [0;0;m*g]-R*[0;0;u(4)];
tau = [u(1); u(2); u(3)];

dq = v;
dth = N\w;
dv = f/m;
dw = inv(J)*(tau - skew(w)*J*w);

xdot = [dq; dth; dv; dw];
end

function X = skew(x)
X = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
end

function Rx =Rx(th)
Rx = [1 0 0;0 cos(th) -sin(th); 0 sin(th) cos(th)];
end

function Ry =Ry(th)
Ry = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
end

function Rz= Rz(th)
Rz = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
end