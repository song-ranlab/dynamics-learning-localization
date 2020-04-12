clc, clear

dt = 0.01;
tf = 10;

g=9.81;
l=1;

theta0 = pi/10;

for ii = 1:tf/dt

t=ii*dt;
theta = theta0*cos(sqrt(g/l)*t);

traj(ii,1)= t;
traj(ii,2)= theta;

end

traj(:,3) = l*sin(traj(:,2));
traj(:,4) = -l*cos(traj(:,2));

plot(traj(:,1),traj(:,2),'o')
legend('Pendulum 1')
title('Theta(t)')
xlabel('time')
ylabel('theta')

plot3(traj(:,1),traj(:,3),traj(:,4))
legend('Pendulum 1')
title('X(t),Y(t)')
xlabel('time')
ylabel('x position')
zlabel('y position')