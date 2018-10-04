function hw25
clc; clear; 
P2
end

function P2
%% Functions
% pdot = q
% qdot = -4*p + 4*q - 7*r - 7*cos(p)

%% Step 1: Find Equillibrium
% Change Equillibirum Point
pe = -2;
qe = 0;

xe = [pe; qe];

%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms r real
re = vpa(solve(@(r) -4*pe + 4*qe - 7*r - 7*cos(pe) == 0, r));
ue = re;

disp('Step 1')
disp(mat2str(xe))
disp(ue)


%% Step 2 Linearize

syms q p r real

%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dp = q;
dq = -4*p + 4*q - 7*r - 7*cos(p);

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
dt = 0.02;

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
x0 = [-2.23;-0.06];

u0 = -K*(x0-xe)+ue;

disp('Step 5')
disp(u0)

%% Step 6: Implement Control Policy (at final time)
t0 = 0;
%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf = 0.06; 

t = t0;
xi = x0;
x = x0;

%%%%%%%%%%%%%%%%%%%%CHANGE THIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-4*p + 4*q - 7*r - 7*cos(p)
gdot = @(t,x,u) [x(2); -4*x(1) + 4*x(2) - 7*u - 7*cos(x(1))];
while t<tf            
    for i = 1:tf/dt
            ui = [-K*(x-xe)+ue];
    end
    % integrate ODE's over the i'th sample interval
    [ti_interval, xi_interval] = ode45(@(t,x) gdot(t,x,double(ui)),[t,t+dt],x);
    x = [xi_interval(end,:)'];
    t = [ti_interval(end)];
end
disp(mat2str(x))
end
