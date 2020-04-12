clc

clear


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
t = -100:0.1:100;
n = length(t);

%observer Y=[x2,y2]
    Y=zeros(2,n);
%     Y(1,:) = 0;
%     Y(2,:) = 0;

    Y(1,:) = 5;
    Y(2,:) = 5;

%     Y(1,:) = t;
%     Y(2,:) = 0;

%feature X=[x1,y1]
    X=zeros(2,n);
%     X(1,:) = 0;
%     X(2,:) = t;

    X(1,:) = t;
    X(2,:) = 0;
    % 
    % X(1,:) = t;
    % X(2,:) = t;
    % 
    % X(1,:) = 0;
    % X(2,:) = 0.5.*t;


%% Analytical Range & Heading
D = X-Y;  %Relative Dynamics
del_x = D(1,:); % x delta
del_y = D(2,:); % y delta
range = sqrt((D(1,:).^2)+(D(2,:).^2)); %total range
heading = atan2(D(2,:),D(1,:)); % heading

%% Observer noise
obs_sigma = 0.01;
obs_noise = obs_sigma * randn(1,n);
obs_sensed = Y + obs_noise;

%% Noisy measurments

D_noise = X-obs_sensed;
range_obs_noise = sqrt((D_noise(1,:).^2)+(D_noise(2,:).^2)); %total range
heading_obs_noise = atan2(D_noise(2,:),D_noise(1,:)); % heading

range_sigma = 0.001;
heading_sigma = 0.01;

range_noise = range_sigma * randn(1,n); 
heading_noise = heading_sigma * randn(1,n); 

range_sensed = range_obs_noise + range_noise;
heading_sensed = heading_obs_noise + heading_noise;

%% Calculate Feature Trajectory

del_y_sensed = range_sensed.*sin(heading_sensed);
del_x_sensed = range_sensed.*cos(heading_sensed);

%% Fit points
p = polyfit(del_x_sensed,del_y_sensed,4);

%% PLOTs
figure(1)
clf

nexttile
plot(Y(1,:), Y(2,:), 'o', X(1,:), X(2,:), 'Linewidth', 2)
%axis([-10 100 -1 2])
legend('observer','feature')
title('Trajectories')
xlabel('x')
ylabel('y')

nexttile
plot(D(1,:), D(2,:), 'k', 'Linewidth', 2)
legend('Relative motion of X')
title('Relative Motion of X viewed from Y')
xlabel('x')
ylabel('y')

nexttile
plot (t,del_x, 'b', t,del_y, 'r', t,range, 'g', 'Linewidth', 2)
legend('x delta','y delta','range')
title('Distance Measurements')
xlabel('time')
ylabel('distance')

nexttile
plot(t, heading, 'k', 'Linewidth', 2)
legend('Heading')
title('Heading Measurements')
xlabel('time')
ylabel('radians')

figure(2)
clf

nexttile
plot (t,range_sensed, 'b', 'Linewidth', 2)
legend('Measured range')
title('Distance Measurement')
xlabel('time')
ylabel('distance')

nexttile
plot (t,heading_sensed, 'b', 'Linewidth', 2)
legend('Measured heading')
title('Heading Measurement')
xlabel('time')
ylabel('radians')

nexttile
plot (t,del_x_sensed, 'b', t,del_y_sensed, 'g', 'Linewidth', 2)
legend('x delta','y delta')
title('Reconstructed Component Measurements')
xlabel('time')
ylabel('distance')

nexttile
%plot(del_x_sensed, del_y_sensed, 'k', 'Linewidth', 2)
scatter(del_x_sensed, del_y_sensed, 'k')
lsline
legend('Relative motion of X')
title('Reconstructed Relative Motion of X viewed from Y')
xlabel('x')
ylabel('y')