function [data, params] = lab4_simulate
    
    % Parameters
   
    % - gravity
    params.g = 9.81;
    % - mass
    params.m = 0.715;
    % - moment of inertia
    params.J = diag([4093e-6, 3944e-6, 7593e-6]);
    % - spar length
    params.l = 0.17;
    % - aerodynamic force and moment coefficients
    params.kF = 7.46e-6;
    params.kM = 1.23e-7;
    % - maximum spin rate
    params.sigmamax = (1e3);
    % - radius of bounding volume (sphere) around quadrotor
    params.r = 3*params.l;
    % - sample time
    params.dt = (1/50);
    
    % Simulation
    % - initial time
    t0 = 0;
    % - initial state
    o0 = [-1.5; 1; -2];
    theta0 = [0; 0; 0];
    v0 = [0; 0; 0];
    w0 = [0; 0; 0];
    x0 = [o0; theta0; v0; w0];
    % - final time
    t1 = 10;
    
    % Problem
    % - desired position
    o_desired = [-1.5; 1.5; -2];
    % - goal position
    o_goal = [1.5; -1.0; -1.5];
    % - obstacles
    %   * create an empty cell array to hold obstacles
    obst = {};
    %   * uncomment this line to add a spherical obstacle (center, radius)
    % obst = AddObstacle_Sphere(obst, [0; 1; -2], 0.5);
    %   * add n random spherical obstacles
    n = 10;
    obst = AddObstacle_RandomSpheres(obst, n, 0.1, 0.5, 2.5, 2.5, 2.5, ...
                                     o_desired, o_goal, params);
	
    % DESIGN
    %
    %     %%%%%%%%%
    %     % FIXME %
    %     %%%%%%%%%
    %
    % - controller
    %
    xe = [0;0;-1;0;0;0;0;0;0;0;0;0;];

    syms o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4 real;

    xs = [o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3]';
    us = [u1 u2 u3 u4]';

    f = h(0,xs,us,params.g,params.m,params.J);


    eqn = h(0,xe,us,params.g,params.m,params.J) == zeros(12,1);
    sol = solve(eqn);
    ue = [double(sol.u1);
        double(sol.u2);
        double(sol.u3);
        double(sol.u4)];

    Ac = double(subs(jacobian(f,xs),[xs;us],[xe;ue]));
    Bc = double(subs(jacobian(f,us),[xs;us],[xe;ue]));


    Ad = expm(Ac*params.dt);
    fun = @(t) expm(Ac*t)*Bc;
    Bd = integral(fun,0,params.dt,'ArrayValue',true);

    Q = diag([250 250 1150 250 250 250 1 1 1 1 1 1]);
    R = diag([50 50 50 50]);

    K = dlqr(Ad, Bd, Q, R);

    params.xe = xe;
    params.ue = ue;
    params.K = K;
        
    %
    % - planner
    params.k_att = 5;
    params.b_att = 5;
    params.k_rep = 1;
    params.b_rep = 1;
    params.k_des = .01;
    params.b_des = .01;
    
    % Data to store
    data.t = [t0];
    data.x = [x0];
    data.u = [];
    data.o_desired = [o_desired];
    data.o_goal = [o_goal];
    data.obst = {obst};
    
    % Iterate over t1/dt sample intervals.
    for i = 1:round(t1/params.dt)

        % Current time and state
        t = data.t(:, end);
        x = data.x(:, end);
        
        % Planner (update o_desired)
        o_desired = planner(t, x, o_desired, o_goal, obst, params);
        
        % Controller (update u)
        u = controller(t, x, o_desired, params);
        
        % Simulator (update t and x)
        [tsol, xsol] = ode45(@(t, x) h(t, x, u, params.g, params.m, params.J), [t t+params.dt], x);
        
        % Store data
        data.o_desired(:, end+1) = o_desired;
        data.u(:, end+1) = u;
        data.t(:, end+1) = tsol(end, :)';
        data.x(:, end+1) = xsol(end, :)';
        data.obst{:, end+1} = obst;

    end

    % Get position and orientation
    data.o = data.x(1:3, :);
    data.theta = data.x(4:6, :);

end

function obst = AddObstacle_Sphere(obst, p, s)
    obst{end+1} = struct('type', 1, 'p', p, 's', s);
end

function obst = ...
    AddObstacle_RandomSpheres(obst, ...
                              n, smin, smax, dx, dy, dz, ...
                              o_start, o_goal, params)

	% n = number of spheres to add
    % smin = minimum radius
    % smax = maximum radius
    % o_start = start position
    % o_goal = goal position
    % (-dx, dx), (-dy, dy), (-dz, 0) = room dimensions

    % Loop to add n spheres one at a time
    for i=1:n

        % Keep trying until we add a sphere that does not interfere
        % with the start or goal position of the quadrotor.
        while(1)
            % Random center position in the room
            p = [dx;dy;dz].*([-1;-1;-1]+[2;2;1].*rand(3,1));
            % Random radius
            s = smin+(smax-smin)*rand;
            % Check for non-interference
            ds = norm(p-o_start)-(params.r+s);
            dg = norm(p-o_goal)-(params.r+s);
            if ((ds>0)&&(dg>0.5))
                obst = AddObstacle_Sphere(obst, p, s);
                break;
            end
        end

    end
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

function u = controller(t, x, o_desired, params)
    
    %%%%%%%%%
    % FIXME %
    %%%%%%%%%
    
    xe = [o_desired;params.xe(4:12)];
    u_desired = -params.K*(x-xe) + params.ue;

    u = GetBoundedInputs( u_desired , params);
    
end

function [u, sigma] = GetBoundedInputs(u_desired, params)

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
    kF = params.kF;
    kM = params.kM;
    l = params.l;
    sigmamax = params.sigmamax;
    
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

function o_desired = planner(t, x, o_desired, o_goal, obst, params)
    
 % Rename o_desired, o_goal, and the radius of the bounding volume for convenience
 q = o_desired;
 q_goal = o_goal;
 r = params.r;

 %{
 % Get attractive part of gradient
 % if ( ... )
 % gradf = ...
 % else
 % gradf = ...
 
% params.k_att = 1;
% params.b_att = 1;
% params.k_rep = 1;
% params.b_rep = 1;
% params.k_des = 1;
% params.b_des = 1;
 %}
 
  % Get attractive part of gradient
 if ( norm(q-q_goal) <= params.b_att )
    gradf = params.k_att*norm(q-q_goal);
 else
    gradf = params.k_att*params.b_att*norm(q-q_goal);
 end
 
 % Get repulsive part of gradient
 for i=1:length(obst)
     % Get center and radius of i'th spherical obstacle
     p = obst{i}.p;
     s = obst{i}.s;
     %{
     % Compute distance and gradient of distance to this obstacle
     % d = ...
     % dgrad = ...
     % Compute repulsive gradient for this obstacle and add it to total gradient
     % if ( ... )
     % gradf = gradf + ...
     %}
     %Compute distance and gradient of distance to this obstacle
     d = norm(q-p)-(s+r);
     dgrad = 1;
     %Compute repulsive gradient for this obstacle and add it to total gradient
    if ( d <= params.b_rep )
        gradf = gradf + params.k_rep*(1/d - 1/params.b_rep)*(-1)/d^2;
    end
 end

 % Take a step
 % if ( ... )

 % q = ...
 % else
 % q = ...

  % Take a step
 if ( params.k_des*gradf <= params.b_des )
    q = q - params.k_des*gradf;
 else
    q = q - params.b_des*(gradf/norm(gradf));
 end
 
 % Recover o_desired

    o_desired = q;

end







