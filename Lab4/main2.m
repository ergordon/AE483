%%
clc; clear; close all;

%plot the actual/desired/obstacle for experiment


load('ExpResults1.mat')

t = table2array(ExpResults1(4999:end,1));
x = table2array(ExpResults1(4999:end,4));
y = table2array(ExpResults1(4999:end,5));
z = table2array(ExpResults1(4999:end,6));

x_desired = table2array(ExpResults1(4999:end,10));
y_desired = table2array(ExpResults1(4999:end,11));
z_desired = table2array(ExpResults1(4999:end,12));

x_obst = table2array(ExpResults1(4999:end,13));
y_obst = table2array(ExpResults1(4999:end,14));
z_obst = table2array(ExpResults1(4999:end,15));

figure()
hold on;
subplot(311)
plot(t,x,t,x_desired,'linewidth',2)
grid on
    title('X Position V. Time')
    xlabel('Time [s]','fontsize',16)
    ylabel('X Position [m]','fontsize',16)
    lgd = legend('Actual X Position', 'Desired X Position')
    set(gca,'fontsize',16)
subplot(312)
plot(t,y,t,y_desired,'linewidth',2)
grid on
    title('Y Position V. Time')
    xlabel('Time [s]','fontsize',16)
    ylabel('Y Position [m]','fontsize',16)
    set(gca,'fontsize',16)
subplot(313)
plot(t,z,t,z_desired,'linewidth',2)
grid on
    title('Z Position V. Time')
    xlabel('Time [s]','fontsize',16)
    ylabel('Z Position [m]','fontsize',16)
    set(gca,'fontsize',16)
    
    set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','ExperimentPositions.pdf');
    
figure()
plot(t,x_obst,t,y_obst,t,z_obst, 'linewidth',2)
grid on
lgd = legend('Obstacle X', 'Obstacle Y', 'Obstacle Z')
    title('Obstacle Position V. Time')
    xlabel('Time [s]','fontsize',16)
    ylabel('Z Position [m]','fontsize',16)
    set(gca,'fontsize',16)

