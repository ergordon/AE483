function simulation_template

    % Parameters (model)
    m = 0.41097000;
    g = 9.81000000;
    J = [0.00411000 0.00000000 0.00000000; 0.00000000 0.00407000 0.00000000; 0.00000000 0.00000000 0.00623000];
    sigmamax = 1000.00000000;
    kF = 0.00001000;
    kM = 0.00000010;
    L = 0.40000000;
    dt = 0.04000000;
    
    % Parameters (simulation)
    x0 = [-0.80; 0.30; -1.40; -0.70; 0.80; -0.90; -0.20; 0.70; 0.00; 0.70; 0.90; 0.70];
    t0 = 0.00;
    t1 = 0.28;

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