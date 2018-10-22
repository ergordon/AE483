function simulation_template

    % Parameters (model)
    g = 10;                 % acceleration of gravity
    m = 1;                  % mass
    J = diag([1, 1, 1]);    % moment of inertia matrix in body frame
    sigmamax = 1e3;         % maximum spin rate of each rotor
    kF = 1e-5;              % aerodynamic force coefficient
    kM = 1e-5;              % aerodynamic torque coefficient
    L = 1;                  % spar length
    dt = 0.1;               % sample time
    
    % Parameters (simulation)
    x0 = zeros(12, 1);      % initial state
    t0 = 0;                 % initial time
    t1 = 1;                 % final time

    % Control policy
    xe = zeros(12, 1);
    ue = zeros(4, 1);
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
        ui = GetBoundedInputs(zeros(4, 1), kF, kM, L, sigmamax);
        u(:, i) = ui;
        
        % Get time and state at start of (i+1)'th sample interval
        [tsol, xsol] = ode45(@(t, x) h(t, x, ui, g, m, J), [ti ti+dt], xi);
        t(:, i+1) = tsol(end, :)';
        x(:, i+1) = xsol(end, :)';

    end
    
end

function xdot = h(t, x, u, g, m, J)

    xdot = zeros(12, 1);

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