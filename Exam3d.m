function Exam3d
clc; clear;
P3
end

function P2
F1 = [6.0 3.0; 3.0 8.0];
g1 = [-4.0; 1.5];
h1 = [0.5];
F2 = [3.0 1.5; 1.5 -1.0];
g2 = [-3.0; -1.5];
h2 = [1.5];

eig(F1)
eig(F2)

Z1 = -(F1\eye(2))*g1
Z2 = -(F2\eye(2))*g2

min1 = 0.5*Z1'*F1*Z1+g1'*Z1 + h1
min2 = 0.5*Z2'*F2*Z2+g2'*Z2 + h2

mat2str(Z1)
end

function P4
Ad = [-4 0; -2 0];
Bd = [6 -7; -4 -3];
Qf = [2 0; 0 6];
Qd = [4 0; 0 8];
Rd = [5 0; 0 6];
Pp = [18 1; 1 16];


F = Bd'*Pp*Bd+Rd
g = 2*Bd'*Pp*Ad
h = Ad'*Pp*Ad+Qd
K = 2*F\g
P = h + g'*(-K)+K'*F*K

mat2str(K)
mat2str(P)
end

function P5
Ad = [3 3 1 -5; 5 6 8 5; -2 6 -8 -9; -8 -7 7 6];
Bd = [-5 6; -3 2; 4 5; -9 -8];
Qf = [9 0 0 0; 0 2 0 0; 0 0 5 0; 0 0 0 2];
Qd = [3 0 0 0; 0 7 0 0; 0 0 5 0; 0 0 0 3];
Rd = [4 0; 0 9];
n = 5;

P = Qf;
for i=1:n
    F = Bd'*P*Bd+Rd
    g = 2*Bd'*P*Ad
    h = Ad'*P*Ad+Qd
    K = 2*F\g
    P = h + g'*(-K)+K'*F*K
end

mat2str(P)
mat2str(K)
end

function P3




end

function simulation_template

 o_desired = [-0.63; -0.83; 0.02];

m = 0.73574000;
g = 9.81000000;
J = [0.00444000 0.00000000 0.00000000; 0.00000000 0.00327000 0.00000000; 0.00000000 0.00000000 0.00664000];
sigmamax = 1000.00000000;
kF = 0.00001000;
kM = 0.00000010;
L = 0.60000000;
dt = 0.01000000;

xe = zeros(12,1);
xe(1:3) = o_desired;
mat2str(xe)

ue = [0;0;0;m*g];
mat2str(ue)

% Data in MATLAB format
Qd = diag([10, 100, 1, 10, 10, 100, 1, 100, 10, 100, 1, 100]);
Rd = diag([100, 100, 1, 1]);

x0 = [0.00; 0.90; -0.40; 0.90; 0.50; 0.40; 0.00; 0.10; -0.40; 0.90; -0.10; -0.70];
t0 = 0.00;
t1 = 0.08;

    % Control policy

    K = zeros(4, 12);

    % Variables to keep track of time, state, and input
    t = [t0];
    x = [x0];
    u = [];
    
    % Iterate over t1/dt sample intervals.
    for i = 1:round(t1/dt)

        % Get time and state at start of i'th sample interval
        ti = t(:, i);
        xi = x(:, i);
        
        % Get input to apply throughout the i'th sample interval that can
        % be realized by spin rates that are between 0 and sigmamax
        ui = GetBoundedInputs(ue, kF, kM, L, sigmamax);
        u(:, i) = ui;
        
        % Get time and state at start of (i+1)'th sample interval
        [tsol, xsol] = ode45(@(t, x) h(t, x, ui, g, m, J), [ti ti+dt], xi);
        t(:, i+1) = tsol(end, :)';
        x(:, i+1) = xsol(end, :)';

    end
    
end

function xdot = h(t, x, u, g, m, J)

o = x(1:3)
th = x(4:6)
v = x(7:9)
w = x(10:12)

tau = [u(1); u(2); u(3)]
R = GetR_ZYX(th)
f = [0;0;m*g]

od =v;
thd = GetThetaDot_ZYX(th, w)
vd = 
wd = 

xdot = [od; thd; vd; wd]
end


function [u, sigma] = GetBoundedInputs(u_desired, kF, kM, L, sigmamax)
    
    u = u_desired;
    sigma = zeros(4,1);

end

function A = wedge(a)

    A = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
    
end

function R = GetR_ZYX(theta)

    c1 = cos(theta(1));
    s1 = sin(theta(1));
    c2 = cos(theta(2));
    s2 = sin(theta(2));
    c3 = cos(theta(3));
    s3 = sin(theta(3));
    R = [c1*c2 c1*s2*s3-s1*c3 c1*s2*c3+s1*s3;
         s1*c2 s1*s2*s3+c1*c3 s1*s2*c3-c1*s3;
         -s2 c2*s3 c2*c3];

end

function thetadot = GetThetaDot_ZYX(theta, w)

    c2 = cos(theta(2));
    s2 = sin(theta(2));
    c3 = cos(theta(3));
    s3 = sin(theta(3));
    thetadot = (1/c2)*[0 s3 c3; 0 c2*c3 -c2*s3; c2 s2*s3 s2*c3]*w;
    
end