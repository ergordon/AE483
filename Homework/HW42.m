function HW42
clc; clear;
P2
end

function P1

end

function P2
k_rep = 0.12;
b_rep = 1.03;
d = 1.67;
dgrad = [0.23; 0.22; 4.01];

if (d <= b_rep)
    dqFrepq = -k_rep*((1/d)-(1/b_rep))*(1/d)^2*dgrad;
else
    dqFrepq = zeros(size(dgrad));
end
mat2str(dqFrepq)
end

function P3
k_att = 0.42;
b_att = 1.09;
q_goal = [0.38; 0.12; 0.56];
k_rep = 0.36;
b_rep = 1.76;
q = [-1.30; 0.64; 0.47];
d = 1.93;
dgrad = [0.05; 4.32; 2.00];

if (norm(q-q_goal) <= b_att)
    dqFattq = k_att*(q-q_goal);
else
    m = b_att;
    dqFattq = m*k_att*((q-q_goal)/norm(q-q_goal));
end

if (d <= b_rep)
    dqFrepq = -k_rep*((1/d)-(1/b_rep))*(1/d)^2*dgrad;
else
    dqFrepq = 0;
end

dqFq =dqFattq + dqFrepq;
mat2str(dqFq)
end