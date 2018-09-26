function hw23
clc;

%p231;
p232;
end

function p231
g = 9.81;
m = 0.44;
theta2 = -0.64;
theta3 = -0.73;
o3 = -4.62;

A = mat2str([0 1; 0 0])

B = mat2str([0; -((cos(theta2)*cos(theta3))/(m))])

end

function p232
g = 9.8100000000;
m = 0.7738987843;
J = [0.0048001031 0.0000000000 0.0000000000; 0.0000000000 0.0041859126 0.0000000000; 0.0000000000 0.0000000000 0.0077481019];
theta1 = 1.3167895686;

iJ = inv(J);

syms th1 th2 th3 u4 real 
Rzz = Rz(cos(th1),sin(th1));
Ryy = Ry(cos(th2),sin(th2));
Rxx = Rx(cos(th3),sin(th3));
    
R = Rzz*Ryy*Rxx;
    
f = [0;0;m*g]+R*[0;0;-u4];

v = vpa(f/m);

jacobian(v,[th1 th2 th3 u4]);

vpa(subs(subs(jacobian(v,[th1 th2 th3 u4]),th1,theta1),[th2 th3 u4],[0 0 (m*g)]));

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

mat2str(A)
mat2str(B)
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