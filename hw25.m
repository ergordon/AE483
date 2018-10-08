function hw25
clc; clear; 
format long
P1
end

function P1()
g = 9.8100000000;
m = 0.7137469554;
J = [0.0037192807 0.0000000000 0.0000000000; 0.0000000000 0.0040796282 0.0000000000; 0.0000000000 0.0000000000 0.0069010173];
o = [0.5629225294; 0.7851218705; -0.5549426634];
theta1 = 1.2690168867;

iJ = inv(J);

syms th1 th2 th3 u4 real 
Rzz = Rz(cos(th1),sin(th1));
Ryy = Ry(cos(th2),sin(th2));
Rxx = Rx(cos(th3),sin(th3));
    
R = subs(Rzz*Ryy*Rxx,th1,theta1);
    
A = [0 0 0 0 0 0 1 0 0 0 0 0;...
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;... 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;...
    0, 0, 0, 0, -cos(theta1)*g, -sin(theta1)*g, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, -sin(theta1)*g, cos(theta1)*g, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

B = [0, 0, 0, 0;...
    0, 0, 0, 0; ...
    0, 0, 0, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0;...
    0, 0, 0, 0;...
    0, 0, 0, -1/m;...
    iJ(1,1), 0, 0, 0;...
    0, iJ(2,2), 0, 0;...
    0, 0, iJ(3,3), 0];



syms th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4 real
x = [o; theta1; th2; th3; v1; v2; v3; w1; w2; w3];
u = [u1; u2; u3; u4];
dx = [0;0;0;0;0;0;0;0;0;0;0;0];
sol = solve([A*x + B*u == dx], [th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4]);
xe  = vpa([o; theta1; sol.th2; sol.th3; sol.v1; sol.v2; sol.v3; sol.w1; sol.w2; sol.w3]);

u4 = vpa(solve(0 == g-((u4*cos(sol.th2)*cos(sol.th3))/(m))));

ue = [sol.u1; sol.u2; sol.u3; u4];

mat2str(double(xe))
mat2str(double(ue))
end

%Make Skew Matrix
function S = skew(w)
S = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

%% Rotation Matrix Maker
function X = Rx(c,s)
    X = [1 0 0; 0 c -s; 0 s c];
end

function Y = Ry(c,s)
    Y= [c 0 s; 0 1 0; -s 0 c];
end

function Z = Rz(c,s)
    Z = [c -s 0; s c 0; 0 0 1];
end

function P2
%% Functions
% pdot = q
% qdot = -9*p^2 - 8 * p - 5*q - 9*r

%% Step 1: Find Equillibrium
% Change Equillibirum Point
pe = 0;
qe = 0;

xe = [pe; qe];

%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms r real
re = vpa(solve(@(r) -9*pe^2 - 8 * pe - 5*qe - 9*r == 0, r));
ue = re;

disp('Step 1')
disp(mat2str(xe))
disp(ue)


%% Step 2 Linearize

syms q p r real

%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dp = q;
dq = -9*p^2 - 8 * p - 5*q - 9*r;

xdot = [dp; dq];
x = [p;q];
u = [r];
Ac = subs(jacobian(xdot,x),p,pe);
Bc = double(jacobian(xdot,u));

disp('Step 2')
disp(mat2str(double(Ac)))
disp(mat2str(Bc))

%% Step 3: Discretize

%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.03;

Ad = double(vpa(expm(Ac*dt)));
Bd = integral(@(tau) double(expm(Ac*tau))*Bc, 0, dt, 'ArrayValued', true);

disp('Step 3')
disp(mat2str(Ad))
disp(mat2str(Bd))


%% Step 4: Design Control Policy
Q = eye(2);
R = 1;

K = dlqr(Ad,Bd,Q,R);
disp('Step 4')
disp(mat2str(K))

%% Step 5: Implement Controly Policy (at initial time)
%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [-0.15; 0.13];

u0 = -K*(x0-xe)+ue;

disp('Step 5')
disp(u0)

%% Step 6: Implement Control Policy (at final time)
t0 = 0;
%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf = 0.18; 

t = t0;
xi = x0;
x = x0;

%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-9*p^2 - 8 * p - 5*q - 9*r
gdot = @(t,x,u) [x(2); -9*x(1).^2-8*x(1)-5*x(2)-9*u];
while t<tf            
    for i = 1:tf/dt
            ui = [-K*(x-xe)+ue];
    end
    % integrate ODE's over the i'th sample interval
    [ti_interval, xi_interval] = ode45(@(t,x) gdot(t,x,double(ui)),[t,t+dt], x);
    x = [xi_interval(end,:)'];
    t = [ti_interval(end)];
end
disp(mat2str(x))
end