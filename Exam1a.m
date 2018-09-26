syms c1 c2 c3 s1 s2 s3 real

theta1 = 1.09;
theta2 = 1.14;
theta3 = -0.89;
theta1dot = -0.94;
theta2dot = -0.25;
theta3dot = -0.36;

c3 = cos(theta3)
c2 = cos(theta2)
c1 = cos(theta1)
s1 = sin(theta1)
s2 = sin(theta2)
s3 = sin(theta3)

Ry = [c1 0 s1; 0 1 0; -s1 0 c1]
Rz = [c2 -s2 0; s2 c2 0; 0 0 1]
Ry1 = [c3 0 s3; 0 1 0; -s3 0 c3]

a = Ry1'*Rz'*[0;1;0]
b = Ry1'*[0;0;1]
c = [0;1;0]

w = [a b c]

mat2str(w*[theta1dot;theta2dot;theta3dot])