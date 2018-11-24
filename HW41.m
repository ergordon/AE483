function HW41
clc; clear; 
P2
end

function P2
q = [-1.18; -3.34; 4.51];
r = 0.24;
p = [-0.46; 3.44; 4.56];
s = 0.22;

dq = norm(p-q) - (s+r)
fgrad = -((p-q)/norm(p-q))

mat2str(dq)
mat2str(fgrad)
end

function P3
qi = [-1.72; 1.37; 0.07];
k_des = 0.62;
fgrad = [3.48; 0.96; -4.38];

qi = qi-k_des*fgrad
mat2str(qi)
end

function P4
k_att = 1.26;
k_rep = 1.49;
k_des = 0.01;
q_goal = [4.91; -0.10; -3.72];
q = [-4.84; 0.15; -2.02];
r = 1.50;
p = [-1.26; -0.14; -3.40];
s = 1.38;
t0 = 0.00;
t1 = 1.00;
dt = 0.05;

dq = norm(p-q) - (s+r)
fgrad = k_att * (q-q_goal) + k_rep*(1/dq^3)*((p-q)/norm(p-q))
qi = q-k_des*fgrad

mat2str(qi)
end
