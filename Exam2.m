function Exam2
clear; clc;
%P1
%[fo tau] = P2();
%P3()
%P4();
[Ac Bc] = P5()
%P7(Ac, Bc);
end

function P1

end

function [fo tau] = P2()
sigmamax = 1000.00000000;
k_F = 0.00001000;
k_M = 0.00000010;
L = 0.70000000;
m = 0.80000000;
g = 9.81000000;
theta1 = -1.24000000;
theta2 = 1.17000000;
theta3 = -0.19000000;
u = [-0.15128923; -0.91774808; 0.04257421; 6.87955850];

R = eul2rotm([theta1 theta2 theta3],'ZYX');

fo = mat2str(-(R*[0;0;u(4)])+[0; 0; m*g]);
tau = mat2str([u(1); u(2); u(3)]);
end

function P3()
g = 9.8100000000;
m = 0.5337289919;
J_in1 = [0.0047390955 0.0000000000 0.0000000000; 0.0000000000 0.0035037555 0.0000000000; 0.0000000000 0.0000000000 0.0074007020];
o_1in0 = [-0.8327047773; 1.2978815142; -1.7229329408];
theta = [-0.1191927319; -0.1153879458; -0.0601766387];
v_01in0 = [0.0019415826; 0.0239960466; -0.0047569729];
w_01in1 = [-0.0472663578; 0.0335819998; -0.0606895701];
u = [0.0009485319; -0.0090036504; 0.0081998214; 5.5388712816];
t0 = 0.2000000000;
t1 = 0.3000000000;
    
    % Define the Time Interval
    t = linspace(t0, t1);
    
    % Define Initial Condition
    x0 = [o_1in0; theta; v_01in0; w_01in1];
       
    % Integrate Angular Rates and Euler's Equations
    [t,x] = ode45(@(t,x) getRates(J_in1, x, u, m, g),t,x0)
    
    q = mat2str(x(end, 1:3)')
    theta = mat2str(x(end, 4:6)')
    v = mat2str(x(end, 7:9)')
    w = mat2str(x(end, 10:12)')
    
    
end

function P4()
g = 9.8100000000;
m = 0.4482015029;
J_in1 = [0.0040387697 0.0000000000 0.0000000000; 0.0000000000 0.0040422979 0.0000000000; 0.0000000000 0.0000000000 0.0063924303];
o_1in0 = [0.4105298925; -0.2078515803; -0.0892999493];
theta = [-0.0766302246; 0.0750366804; -0.0352032279];
v_01in0 = [-0.0293432681; -0.0037488968; -0.0537753836];
w_01in1 = [0.0363990085; -0.0187535408; -0.0618633961];
u = [0.0002745712; -0.0060953163; -0.0157263598; 5.0761417090];

x = [o_1in0; theta; v_01in0; w_01in1];

dx = getRates(J_in1, x, u, m, g);

dq = mat2str(dx(1:3))
dtheta = mat2str(dx(4:6))
dv = mat2str(dx(7:9))
dw = mat2str(dx(10:12))
end

function [Ac Bc] = P5()
g = 9.8100000000;
m = 0.7748274155;
J = [0.0043750967 0.0000000000 0.0000000000; 0.0000000000 0.0043330862 0.0000000000; 0.0000000000 0.0000000000 0.0065350919];
theta1 = 2.3329674107;

syms o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4 real

x = [o1; o2; o3; th1; th2; th3; v1; v2; v3; w1; w2; w3];
xe = [0; 0; 0; theta1; 0; 0; 0; 0; 0; 0; 0; 0];
u = [u1; u2; u3; u4];
ue = [0;0;0;m*g];

dx = getRates(J, x, u, m, g);

dq = dx(1:3);
dtheta = dx(4:6);
dv = dx(7:9);
dw = dx(10:12);

%Ac = double(vpa(subs(subs(jacobian(dx,x),x,xe),u,ue)));
%Bc = double(vpa(subs(subs(jacobian(dx,u),x,xe),u,ue)));

Ac = double(vpa(subs(jacobian(dx,x),[x;u],[xe;ue])));
Bc = double(vpa(subs(jacobian(dx,u),[x;u],[xe;ue])));
end

function P6(Ac, Bc)
g = 9.8100000000;
m = 0.5869135260;
J = [0.0035404997 0.0000000000 0.0000000000; 0.0000000000 0.0047818889 0.0000000000; 0.0000000000 0.0000000000 0.0077278919];
o = [-0.5468813884; 0.3336393308; 0.0402405889];
theta1 = 0.0301008677;

syms th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4 real
x = [o; theta1; th2; th3; v1; v2; v3; w1; w2; w3];
u = [u1; u2; u3; u4];
dx = [0;0;0;0;0;0;0;0;0;0;0;0];
sol = solve([Ac*x + Bc*u == dx], [th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4]);
xe  = vpa([o; theta1; sol.th2; sol.th3; sol.v1; sol.v2; sol.v3; sol.w1; sol.w2; sol.w3]);

u4 = vpa(solve(0 == g-((u4*cos(sol.th2)*cos(sol.th3))/(m))));

ue = [sol.u1; sol.u2; sol.u3; u4];

mat2str(double(xe))
mat2str(double(ue))
end

function P7(Ac,Bc)
g = 9.8100000000;
m = 0.7748274155;
J = [0.0043750967 0.0000000000 0.0000000000; 0.0000000000 0.0043330862 0.0000000000; 0.0000000000 0.0000000000 0.0065350919];
theta1 = 2.3329674107;
dt = 0.9100000000;

Ad = expm(Ac*dt);
Bd = integral(@(tau) f(Ac, tau, Bc),0,dt,'ArrayValued',true);

mat2str(Ad)
mat2str(Bd)
end


function f = f(Ac, tau, Bc)
    f= expm(Ac*tau)*Bc;
end


function xdot = getRates(J, x, u, m, g)
    q = x(1:3);
    theta = x(4:6);
    v = x(7:9);
    w = x(10:12);
   
    tau = [u(1); u(2); u(3)];
    
    Rzz = Rz(cos(theta(1)),sin(theta(1)));
    Ryy = Ry(cos(theta(2)),sin(theta(2)));
    Rxx = Rx(cos(theta(3)),sin(theta(3)));

    N = [Rxx'*Ryy'*[0;0;1] Rxx'*[0;1;0] eye(3)*[1;0;0]];
    f = Rzz*Ryy*Rxx*[0;0;-u(4)]+[0;0;m*g];
    
    dq = v;
    dtheta = N\w;
    dv = f/m;
    dw = inv(J)*(tau-skew(w)*J*w);
      
    xdot = [dq; dtheta; dv; dw];
end

%Make Skew Matrix
function S = skew(w)
S = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

function X = Rx(c,s)
    X = [1 0 0; 0 c -s; 0 s c];
end

function Y = Ry(c,s)
    Y= [c 0 s; 0 1 0; -s 0 c];
end

function Z = Rz(c,s)
    Z = [c -s 0; s c 0; 0 0 1];
end