function [t, o, theta, odes] = lab3_simulate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   MODIFY
    %
    %   - Change parameter values and sample time to match your quadrotor (see
    %     your results from Labs #1 and #2)
    %   - Change initial time, final time, and initial conditions as you like
    %

    % Parameters
    g = 9.81;                 % acceleration of gravity
    m = 0.715;                  % mass
    J = diag([4093e-6, 3944e-6, 7593e-6]);    % moment of inertia matrix in body frame
    l = 0.17;                  % spar length
    kF = 7.46e-6;              % aerodynamic force coefficient
    kM = 1.23e-7;              % aerodynamic torque coefficient
    sigmamax = 1e3;         % maximum spin rate

    % Initial time
    t0 = 0;
    
    o0 = [0; 0; -2];
    theta0 = [0; 0; 0];
    v0 = [0; 0; 0];
    w0 = [0; 0; 0];

    % Sample time
    dt = (1/50);

    % Final time
    t1 = 10;


    xe = [0;0;-1;0;0;0;0;0;0;0;0;0;];

    syms o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4 real;

    xs = [o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3]';
    us = [u1 u2 u3 u4]';

    f = h(0,xs,us,g,m,J);

    
    eqn = h(0,xe,us,g,m,J) == zeros(12,1);
    sol = solve(eqn);
    ue = [double(sol.u1);
        double(sol.u2);
        double(sol.u3);
        double(sol.u4)];

    Ac = double(subs(jacobian(f,xs),[xs;us],[xe;ue]));
    Bc = double(subs(jacobian(f,us),[xs;us],[xe;ue]));


    Ad = expm(Ac*dt)
    fun = @(t) expm(Ac*t)*Bc;
    Bd = integral(fun,0,dt,'ArrayValue',true)

    Q = [500 0 0 0 0 0 0 0 0 0 0 0;
         0 500 0 0 0 0 0 0 0 0 0 0;
         0 0 500 0 0 0 0 0 0 0 0 0;
         0 0 0 250 0 0 0 0 0 0 0 0;
         0 0 0 0 250 0 0 0 0 0 0 0;
         0 0 0 0 0 250 0 0 0 0 0 0;
         0 0 0 0 0 0 1 0 0 0 0 0;
         0 0 0 0 0 0 0 1 0 0 0 0;
         0 0 0 0 0 0 0 0 1 0 0 0;
         0 0 0 0 0 0 0 0 0 1 0 0;
         0 0 0 0 0 0 0 0 0 0 1 0;
         0 0 0 0 0 0 0 0 0 0 0 1];
    R = [50 0 0 0 ;
         0 50 0 0;
         0 0 50 0;
         0 0 0 50];
    K = dlqr(Ad, Bd, Q, R);
   

    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create variables to keep track of time, state, input, and desired position
    t = [t0];
    x = [o0; theta0; v0; w0];
    u = [];
    odes = [];
%      times = t0:dt:t1;
%      for i=1:length(times)
%          % odes = [odes mickymouse((pi*times(i))/30)];
%      end
     times = t0:dt:t1;
     for i=1:length(times)
         odes = [odes [0;0;-1]];
     end     

    % Iterate over t1/dt sample intervals.
    for i = 1:(t1/dt)

        % Get time and state at start of i'th sample interval
        ti = t(:, i);
        xi = x(:, i);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %   MODIFY
        %
        %   - Get desired position at start of i'th sample interval
        %
        %     (You may also need to redefine your equilibrium point!)
        %

        odesi = odes(:,i);

        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         odes(:, i) = odesi;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %   MODIFY
        %
        %   - Get input that will be applied throughout the i'th sample
        %     interval (in other words, implement your control policy)
        %   - Don't forget to make sure that this input can be realized by
        %     spin rates that are between 0 and sigmamax
        %
        
%         u_desired = -K*(xi-xe) + ue;
        xe = [odesi;xe(4:12)];
        u_desired = -K*(xi-xe) + ue;

        ui = GetBoundedInputs( u_desired , kF, kM, l, sigmamax);

        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u(:, i) = ui;

        % Get time and state at start of (i+1)'th sample interval
        [tsol, xsol] = ode45(@(t, x) h(t, x, ui, g, m, J), [ti ti+dt], xi);
        t(:, i+1) = tsol(end, :)';
        x(:, i+1) = xsol(end, :)';

    end

    % Get position and orientation
    o = x(1:3, :);
    theta = x(4:6, :);

end

function xdot = h(t, x, u, g, m, J)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   MODIFY
    %
    %   - Compute xdot given x, u, and other parameters (see HW2.2)

    
    o = x(1:3);
    t = x(4:6);
    v = x(7:9);
    w = x(10:12);
    
    what = [0 -w(3) w(2);
            w(3) 0 -w(1);
            -w(2) w(1) 0];
    %ZYX
    R = Rz(t(1))*Ry(t(2))*Rx(t(3));
    
    N = [Rx(t(3))'*Ry(t(2))'*[0;0;1] Rx(t(3))'*[0;1;0] [1;0;0]]^(-1);
    
    odot = v;
    tdot = N*w;
    vdot = 1/m*([0;0;m*g] + R*[0;0;-u(4)]);
    wdot = J^(-1)*(u(1:3) - what*J*w);
    
    xdot = [odot;tdot;vdot;wdot];

    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


function [u, sigma] = GetBoundedInputs(u_desired, kF, kM, l, sigmamax)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   MODIFY
    %
    %   - Compute an input u (4x1) that is close to the input u_desired but
    %     that can be realized by spin rates that are between 0 and sigmamax
    %   - Compute the spin rates sigma (4x1) - these are not used by the rest
    %     of the simulation code, but will be necessary in your onboard C code,
    %     so it is useful to make sure you know how to compute them here
    %
    %     (See HW2.1.2)

    u = u_desired;
    W = [l*kF -l*kF 0 0; 0 0 l*kF -l*kF; kM kM -kM -kM; kF kF kF kF];
    sigma_d = W^-1*u;
    sigma_d(sigma_d > sigmamax^2) =  sigmamax^2;
    sigma_d(sigma_d < 0) =  0;
    sigma = sigma_d;
    u = W*sigma;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
