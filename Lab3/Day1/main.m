clc; clear; close all;

[t, o, theta, odes] = lab3_simulate();

%%
%moviefile = 'mickeyMouse.mp4';
%lab3_visualize(t, o, theta, odes, moviefile)

%%
figure()
plot(o(1,:),o(2,:),'b-');
hold on
plot(odes(1,:),odes(2,:),'r--');
title('xy')

figure()
plot(t,o(1,:),'r-',t,o(2,:),'b-',t,o(3,:),'g-');
hold on
plot(t(2:end),odes(1,:),'r--',t(2:end),odes(2,:),'b--',t(2:end),odes(3,:),'g--');
legend('x','y','z','x_des','y_des','z_des');
title('position')


figure()
plot(t,theta(1,:),t,theta(2,:),t,theta(3,:));
legend('t1','t2','t3');
title('euler angles')
