function hello
    % CODE TO CONVERT FROM XYZ EULER ANGLES TO A ROTATION MATRIX
    theta1 = -2.63;
    theta2 = 1.73;
    theta3 = 0.17;

    c1 = cos(theta1);
    c2 = cos(theta2);
    c3 = cos(theta3);

    s1 = sin(theta1);
    s2 = sin(theta2);
    s3 = sin(theta3);

    %XZX
    %R =  z(c1,s1) * y(c2, s2) * z(c3,s3)
    R =  x(c1,s1)
end

function X = x(c,s)
    X = [1 0 0; 0 c -s; 0 s c];
end

function Y = y(c,s)
    Y= [c 0 s; 0 1 0; -s 0 c];
end

function Z = z(c,s)
    Z = [c -s 0; s c 0; 0 0 1];
end