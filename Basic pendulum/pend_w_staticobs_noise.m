clc
clear all
hold off

% Code exploring DLL with a Feature X and an Observer Y.
% The obsever can be held at a fixed point, or given its own trajectory.
% Values/Measurements/parameters to consider:
% Range
% Heading
% Velocity of observer/feature (constatnt?)
% Relative velocity between observer/feature
% Rcceleration of observer/feature (constatnt?)
% Relative acceleration between observer/feature
% Subtracting the observer from the feature gives the motion as if the
% observer is the origin.

% add noise
% rebuild feature trajectory from sensor measurement

%observer noise
%extract dynamics
%Mass spring damper
%
%% Initialize
%t = -100:0.1:100;
g = 9.81;                   %acceration of gravity
l=1;                        % pendulum length
N = 1;                      %number of oscilations
dt = 0.01;                  % Time step
tf = 0.5*N*(2*pi*sqrt(l/g));% length of simulation, multiple of half period
tf = round(tf);             %final time of simulation, rounded for ease of indexing
n = tf/dt;                  %number of iterations


%observer Y=[x1,y1]
    Y=zeros(2,n); % initialize observer

    Y(1,:) = 0.4;
    Y(2,:) = -0.98;

%feature X=[x2,y2]
    X=zeros(2,n); % initialize feature trajectory

    theta0 = pi/10; %initial pedulum theta

    for ii = 1:tf/dt
        t=ii*dt;
        theta = theta0*cos(sqrt(g/l)*t); % Pendulum EOM (small order approx)

        traj(ii,1)= t;
        traj(ii,2)= theta;

    end

    traj(:,3) = l*sin(traj(:,2));
    traj(:,4) = -l*cos(traj(:,2));

    X(1,:) = traj(:,3);
    X(2,:) = traj(:,4);
    T = traj(:,1).'; %time vector
    
    %X(1,:) = t;
    %X(2,:) = 0;
    % 
    %X(1,:) = t;
    %X(2,:) = t;
    
    figure(1)
    clf
    
    plot(traj(:,1),traj(:,2))
    legend('Pendulum 1')
    title('Theta(t)')
    xlabel('time')
    ylabel('theta')
    
    figure(2)
    clf
    plot3(traj(:,1),traj(:,3),traj(:,4),traj(:,1),Y(1,:),Y(2,:),'linewidth',1.5 )
    legend('Pendulum 1', 'Observer')
    title('X(t),Y(t)')
    xlabel('time')
    ylabel('x position')
    zlabel('y position')
    grid on
    
    
%% Analytical Range & Heading
D = X-Y;  %Relative Dynamics
del_x = D(1,:); % x delta
del_y = D(2,:); % y delta
range = sqrt((D(1,:).^2)+(D(2,:).^2)); %total range
heading = atan2(D(2,:),D(1,:)); % heading

for ii = 1:tf/dt
if heading(ii) <= 0
    heading(ii) = (2*pi)-abs(heading(ii));
else
end
end
%% Observer noise
obs_sigma = 0.001;
obs_noise = obs_sigma * randn(1,n);
obs_sensed = Y + obs_noise;

%% Noisy measurments

D_noise = X-obs_sensed;
range_obs_noise = sqrt((D_noise(1,:).^2)+(D_noise(2,:).^2)); %total range
heading_obs_noise = atan2(D_noise(2,:),D_noise(1,:)); % heading

range_sigma = 0.01;
heading_sigma = 0.01;

range_noise = range_sigma * randn(1,n); 
heading_noise = heading_sigma * randn(1,n); 

range_sensed = range_obs_noise + range_noise;
heading_sensed = heading_obs_noise + heading_noise;
for ii = 1:tf/dt
if heading_sensed(ii) <= 0
    heading_sensed(ii) = (2*pi)-abs(heading_sensed(ii));
else
end
end
%% Calculate Feature Trajectory

del_y_sensed = range_sensed.*sin(heading_sensed);
del_x_sensed = range_sensed.*cos(heading_sensed);

%% Fit points
p = polyfit(del_x_sensed,del_y_sensed,5);
xspace = del_x_sensed(1,n):dt:del_x_sensed(1,1);
Polycurve = polyval(p,xspace);

pX = polyfit(T,del_x_sensed,5);
PolyX = polyval(pX,T);

pY = polyfit(T,del_y_sensed,5);
PolyY = polyval(pY,T);

%% Velocity and Acceleraiton profiles
for jj = 1:tf/dt-1
Yvel(jj) = (PolyY(jj+1) - PolyY(jj))/dt;
end

for kk = 1:tf/dt-1
Xvel(kk) = (PolyX(kk+1) - PolyX(kk))/dt;
end

Vel = sqrt(Yvel.^2 + Xvel.^2);

for mm = 1:tf/dt-2
Yacel(mm) = (Yvel(mm+1) - Yvel(mm))/dt;
end

for nn = 1:tf/dt-2
Xacel(nn) = (Xvel(nn+1) - Xvel(nn))/dt;
end

Acel = sqrt(Yacel.^2 + Xacel.^2);


%% PLOTs
figure(3)
clf

nexttile
plot(Y(1,:), Y(2,:), 'o', X(1,:), X(2,:), 'Linewidth', 1)
%axis([-10 100 -1 2])
legend('observer','feature')
title('Trajectories')
xlabel('x')
ylabel('y')

nexttile
plot(D(1,:), D(2,:), 'k', 'Linewidth', 1)
legend('Relative motion of X')
title('Relative Motion of X viewed from Y')
xlabel('x')
ylabel('y')

nexttile
plot (T,del_x, 'b', T,del_y, 'r', T,range, 'g', 'Linewidth', 1)
legend('x delta','y delta','range')
title('Distance Measurements')
xlabel('time')
ylabel('distance')

nexttile
plot(T, heading, 'k', 'Linewidth', 1)
legend('Heading')
title('Heading Measurements')
xlabel('time')
ylabel('radians')

figure(4)
clf

nexttile
plot (T,range_sensed, 'b', 'Linewidth', 1)
legend('Measured range')
title('Distance Measurement')
xlabel('time')
ylabel('distance')

nexttile
plot (T,heading_sensed, 'b', 'Linewidth', 1)
legend('Measured heading')
title('Heading Measurement')
xlabel('time')
ylabel('radians')

nexttile
plot (T,del_x_sensed, 'b', T,del_y_sensed, 'g', T, PolyX, 'r', T, PolyY, 'K', 'Linewidth', 1)
legend('x delta','y delta')
title('Reconstructed Component Measurements')
xlabel('time')
ylabel('distance')

nexttile
plot(xspace, Polycurve, 'Linewidth', 1)
hold on
scatter(del_x_sensed, del_y_sensed, 'k')
%lsline
legend('Relative motion of X')
title('Reconstructed Relative Motion of X viewed from Y')
xlabel('x')
ylabel('y')

figure(5)
clf

plot (T(2:n), Yvel,'r', T(2:n), Xvel, 'b', 'Linewidth', 1)
legend('Y Velocity', 'Xvelocity')
title('Velocity Profiles')
xlabel('time')
ylabel('Velocity')

figure(6)
clf

plot (T(3:n), Yacel,'r', T(3:n), Xacel, 'b', T(3:n), Acel, 'g',  'Linewidth', 1)
legend('Y Acceleration', 'X Acceleration', 'Magnitude')
title('Acceleration Profiles')
xlabel('time')
ylabel('Acceleration')
