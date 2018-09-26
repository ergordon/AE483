function hw22
clc;
%p221;
%p222;
%p223;
%p224'
end

%% Problem 2.2.1
function p221
    g = 9.8100000000;
    m = 0.7884295359;
    J_in1 = [0.0033016970 0.0000000000 0.0000000000; 0.0000000000 0.0045141760 0.0000000000; 0.0000000000 0.0000000000 0.0070462466];
    o_1in0 = [-0.2681102596; -0.4454069188; -1.0110801983];
    theta = [-0.1072574566; -0.0167450674; 0.0736743339];
    v_01in0 = [0.0424034614; 0.0691002554; 0.0221973885];
    w_01in1 = [-0.0273341945; 0.1348373512; -0.0267022083];
    u = [-0.0062083254; 0.0056255452; -0.0069389030; 5.2760775180];

    tau = [u(1); u(2); u(3)];
    
    Rzz = Rz(cos(theta(1)),sin(theta(1)));
    Ryy = Ry(cos(theta(2)),sin(theta(2)));
    Rxx = Rx(cos(theta(3)),sin(theta(3)));


    R1 = Rxx'*Ryy';
    R2 = Rxx';
    R3 = eye(3)';
    
    N = [R1*[0;0;1] R2*[0;1;0] R3*[1;0;0]];
    f = Rzz*Ryy*Rxx*[0;0;-u(4)]+[0;0;m*g];
    
    %%Answers
    odot = mat2str(v_01in0)
    thetadot = mat2str(N\w_01in1)
    vdot = mat2str(f/m)
    wdot = mat2str(inv(J_in1)*(tau-skew(w_01in1)*J_in1*w_01in1))
end

%% Problem 2
function p222
g = 9.8100000000;
m = 0.4740181293;
J_in1 = [0.0031807557 0.0000000000 0.0000000000; 0.0000000000 0.0039471198 0.0000000000; 0.0000000000 0.0000000000 0.0072035444];
o_1in0 = [1.5288527023; -1.6938940209; 0.1480760158];
theta = [0.0658567281; -0.1302761989; 0.0373092197];
v_01in0 = [0.1088819380; -0.0446576492; 0.0273892616];
w_01in1 = [0.0111772510; 0.1131011919; 0.0075007393];
u = [-0.0091671459; -0.0004433811; -0.0028581238; 6.2127181537];
t0 = 0.9000000000;
t1 = 1.3000000000;

    tau = [u(1); u(2); u(3)];
    
    Rzz = Rz(cos(theta(1)),sin(theta(1)));
    Ryy = Ry(cos(theta(2)),sin(theta(2)));
    Rxx = Rx(cos(theta(3)),sin(theta(3)));


    R1 = Rxx'*Ryy';
    R2 = Rxx';
    R3 = eye(3)';
    
    N = [R1*[0;0;1] R2*[0;1;0] R3*[1;0;0]];
    f = Rzz*Ryy*Rxx*[0;0;-u(4)]+[0;0;m*g];

    % Define the Time Interval
    t = linspace(t0, t1);
    
    % Define Initial Condition
    x0 = [o_1in0; theta; v_01in0; w_01in1]
    
   
    
    % Integrate Angular Rates and Euler's Equations
    [t,x] = ode45(@(t,x) xdot1(J_in1, x, u, m, g),t,x0);
    
    q = mat2str(x(end, 1:3)')
    theta = mat2str(x(end, 4:6)')
    v = mat2str(x(end, 7:9)')
    w = mat2str(x(end, 10:12)')
end

%ODE 45 Function
function xdot = xdot1(J, x, u, m, g)
    q = x(1:3);
    theta = x(4:6);
    v = x(7:9);
    w = x(10:12);
    
    
    tau = [u(1); u(2); u(3)];
    
    Rzz = Rz(cos(theta(1)),sin(theta(1)));
    Ryy = Ry(cos(theta(2)),sin(theta(2)));
    Rxx = Rx(cos(theta(3)),sin(theta(3)));


    R1 = Rxx'*Ryy';
    R2 = Rxx';
    R3 = eye(3)';
    
    N = [R1*[0;0;1] R2*[0;1;0] R3*[1;0;0]];
    f = Rzz*Ryy*Rxx*[0;0;-u(4)]+[0;0;m*g];
    
        dq = v;
        dtheta = N\w;
        dv = f/m;
        dw = inv(J)*(tau-skew(w)*J*w);
      
    xdot = [dq; dtheta; dv; dw];
end

%% Problem 3
function p223
g = 9.8100000000;
m = 0.6548233901;
J = [0.0033023256 0.0000000000 0.0000000000; 0.0000000000 0.0038182521 0.0000000000; 0.0000000000 0.0000000000 0.0077020080];
o = [1.6832810573; 0.4664217887; -0.9376381731];
theta = [-0.0927721369; 0.0268658940; -0.1404935697];
v = [-0.0074800033; 0.0042389697; -0.0110375549];
w = [0.0121601878; -0.0221588857; -0.0238367838];
kp1 = 2.0593057641;
kv1 = 2.8858475996;
kp2 = 0.7337656408;
kv2 = 2.7490681609;
kp3 = 2.3095506788;
kv3 = 1.6948042837;
kp4 = -2.3025848171;
kv4 = -0.8947561143;

u = [-kp1*theta(3)-kv1*w(1);...
    -kp2*theta(2)-kv2*w(2);...
    -kp3*theta(1)-kv3*w(3);...
    m*g-kp4*o(3)-kv4*v(3);]

    tau = [u(1); u(2); u(3)];
    
    Rzz = Rz(cos(theta(1)),sin(theta(1)));
    Ryy = Ry(cos(theta(2)),sin(theta(2)));
    Rxx = Rx(cos(theta(3)),sin(theta(3)));


    R1 = Rxx'*Ryy';
    R2 = Rxx';
    R3 = eye(3)';
    
    N = [R1*[0;0;1] R2*[0;1;0] R3*[1;0;0]];
    f = Rzz*Ryy*Rxx*[0;0;-u(4)]+[0;0;m*g];
    
    %%Answers
    odot = mat2str(v)
    thetadot = mat2str(N\w)
    vdot = mat2str(f/m)
    wdot = mat2str(inv(J)*(tau-skew(w)*J*w))
end

%% Problem 4
function p224
g = 9.8100000000;
m = 0.6443725805;
J = [0.0031339102 0.0000000000 0.0000000000; 0.0000000000 0.0039908001 0.0000000000; 0.0000000000 0.0000000000 0.0072814037];
o = [-1.0733231753; 0.6947016552; -0.5046010705];
theta = [0.0107339694; 0.0074542375; 0.0039804912];
v = [-0.0388085809; 0.0688175008; -0.0218117937];
w = [0.0002586294; -0.0207095141; 0.0524144159];
t0 = 0.3000000000;
t1 = 0.4000000000;
kp1 = 1.1753060167;
kv1 = 2.5714897045;
kp2 = 2.5545979482;
kv2 = 0.8486010536;
kp3 = 0.6228952611;
kv3 = 1.7975392795;
kp4 = -1.0537905769;
kv4 = -2.4523115606;

u0 = [-kp1*theta(3)-kv1*w(1);...
    -kp2*theta(2)-kv2*w(2);...
    -kp3*theta(1)-kv3*w(3);...
    m*g-kp4*o(3)-kv4*v(3);]

    tau = [u0(1); u0(2); u0(3)];
    
    Rzz = Rz(cos(theta(1)),sin(theta(1)));
    Ryy = Ry(cos(theta(2)),sin(theta(2)));
    Rxx = Rx(cos(theta(3)),sin(theta(3)));


    R1 = Rxx'*Ryy';
    R2 = Rxx';
    R3 = eye(3)';
    
    N = [R1*[0;0;1] R2*[0;1;0] R3*[1;0;0]];
    f = Rzz*Ryy*Rxx*[0;0;-u0(4)]+[0;0;m*g];

    % Define the Time Interval
    t = linspace(t0, t1);
    
    % Define Initial Condition
    x0 = [o; theta; v; w]
    
   
    
    % Integrate Angular Rates and Euler's Equations
    [t,x] = ode45(@(t,x) xdot2(J, x, m, g,kp1,kv1,kp2,kv2,kp3,kv3,kp4,kv4),t,x0);
    
    %answer
    q = mat2str(x(end, 1:3)')
    theta = mat2str(x(end, 4:6)')
    v = mat2str(x(end, 7:9)')
    w = mat2str(x(end, 10:12)')
end

%ODE 45 Function
function xdot = xdot2(J, x, m, g,kp1,kv1,kp2,kv2,kp3,kv3,kp4,kv4)
    q = x(1:3);
    theta = x(4:6);
    v = x(7:9);
    w = x(10:12);
    
    u = [-kp1*theta(3)-kv1*w(1);...
    -kp2*theta(2)-kv2*w(2);...
    -kp3*theta(1)-kv3*w(3);...
    m*g-kp4*q(3)-kv4*v(3);]
    
    tau = [u(1); u(2); u(3)];
    
    Rzz = Rz(cos(theta(1)),sin(theta(1)));
    Ryy = Ry(cos(theta(2)),sin(theta(2)));
    Rxx = Rx(cos(theta(3)),sin(theta(3)));


    R1 = Rxx'*Ryy';
    R2 = Rxx';
    R3 = eye(3)';
    
    N = [R1*[0;0;1] R2*[0;1;0] R3*[1;0;0]];
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