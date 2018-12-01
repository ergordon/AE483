function [x_final,K,xe,ue]=simulationtemp()
o_desired = [0.23; 0.34; 0.07];

m = 0.75677000;
g = 9.81000000;
J = [0.00402000 0.00000000 0.00000000; 0.00000000 0.00472000 0.00000000; 0.00000000 0.00000000 0.00664000];
sigmamax = 1000.00000000;
kF = 0.00001000;
kM = 0.00000010;
L = 0.80000000;
dt = 0.10000000;

% Data in MATLAB format
Qd = diag([100, 100, 100, 1, 100, 100, 10, 1, 1, 100, 10, 10]);
Rd = diag([1, 100, 10, 10]);

x0 = [0.30; 0.60; -1.00; 0.80; 0.70; 0.50; -0.50; 0.20; 0.30; -0.70; -0.90; 0.50];
t0 = 0.00;
t1 = 0.60;

% Control policy
xe = zeros(12, 1);
xe(1:3)=o_desired;
ue = zeros(4, 1);
ue(4)=m*g;

mat2str(xe)
mat2str(ue)
 

syms o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4 real
xc = [o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3]';
uc = [u1 u2 u3 u4]';

Ac = double(subs(subs(jacobian(h([],xc,uc,g,m,J),xc),xc,xe),uc,ue))
Bc = double(subs(subs(jacobian(h([],xc,uc,g,m,J),uc),xc,xe),uc,ue))


[Ad, Bd] = ssdata(c2d(ss(Ac,Bc,[],[]),dt));
K = dlqr(Ad,Bd,Qd,Rd);
mat2str(K)

  %{
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
        [tsol, xsol] = ode45(@(t, x) h(t, xi, u, g, m, J), [ti ti+dt], xi);
        t(:, i+1) = tsol(end, :)';
        x(:, i+1) = xsol(end, :)';
        XX = xsol(end, :)';
    end
    %}


t = 0;
x = x0;

    while t<t1
    ui = GetBoundedInputs(ue, kF, kM, L, sigmamax);
    [tt,xx] = ode45(@(t, x) h([], x, ui, g, m, J), [t t+dt], x);
    x = xx(end,:)';
    t = tt(end,:)+dt;
    end

    mat2str(x)
    % -0.191673894439; 0.795535418425; -0.43072232135; 0.471444110879; -0.307702775584; -0.103689818148; -0.177498072134; 0.583649368812; 1.09408265825; -0.354592403508; 0.510850438834; -0.149092179936]
end

function xdot = h(t,x,u,g,m,J)
q= x(1:3);
th= x(4:6);
v= x(7:9);
w = x(10:12);

R = GetR_ZYX(th);
f = [0;0;m*g]-R*[0;0;u(4)];
tau = [u(1); u(2); u(3)];

dq = v;
dth = GetThetaDot_ZYX(th, w);
dv = f/m;
dw = inv(J)*(tau - w'*J*w);

xdot = [dq; dth; dv; dw];
end


function [u, sigma] = GetBoundedInputs(u_desired, kF, kM, L, sigmamax)
    W = [kF*L 0 -kF*L 0; 0 kF*L 0 -kF*L; kM -kM kM -kM; kF kF kF kF];
    sigma = W\u_desired;
    sigma = ((sigma>(sigmamax^2)).*((sigmamax^2)*ones(size(sigma)))+(sigma<=(sigmamax^2)).*sigma).*(sigma>0);
    u = W*sigma;
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