clc
clear
set(0, 'DefaultLineLineWidth', 1.75);

% pendulum parameters 
m = 1;  
g = 1;
l = 1;

% simulation time parameters
tf = 10;   
dt = 0.05;
n = tf/dt + 1;

% noise parameters
omega_sigma = 0.3; %process sensor noise
theta_sigma = 0.01; % Position measurement noise
z_sigma = 0.1; %DLL measurement noise
S = 0; %kalman filter initial covariance

poly_order = 6; % polynomial regression order

x_0 = [pi/3,0,0];    %Initial state of observer in the form [theta_0,omega_0,z]
y_0 = [pi/4,0,2];    %Initial state of Feature in the form [theta_0,omega_0,z]

%initialize data matrices
T = [0:dt:tf];  
X = zeros(n,3);
Y = zeros(n,3);

%% Simulate Ground Truth
X(1,:)=x_0; % initial state x
Y(1,:)=y_0; % initial state y

for ii=2:n %propogate analytical state vectors for the whole time vector
   %Observer
   X(ii,:) = [X(ii-1,1) + X(ii-1,2)*dt, X(ii-1,2)-(g/l^2)*sin(X(ii-1,1))*dt,X(ii-1,3)]; 
   %Feature
   Y(ii,:) = [Y(ii-1,1) + Y(ii-1,2)*dt, Y(ii-1,2)-(g/l^2)*sin(Y(ii-1,1))*dt,Y(ii-1,3)]; 
    
end    

%plot observer analytical solution
figure(1)
clf
plot(T,X(:,1),'b',T,X(:,2),'r--')
legend('$\theta(t)$','$\omega(t)$','fontsize', 18, 'interpreter', 'latex')
title('Observer Analytical Solution','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('radians,$\frac{radians}{s}$','fontsize', 18, 'interpreter', 'latex')

%Plot feature analytical solution
figure(2)
clf
plot(T,Y(:,1),'b',T,Y(:,2),'r--')
legend('$\theta(t)$','$\omega(t)$','fontsize', 18, 'interpreter', 'latex')
title('Feature Analytical Solution','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('radians,$\frac{radians}{s}$','fontsize', 18, 'interpreter', 'latex')

%Convert to cartesian coordinates
X_cart = [l*sin(X(:,1)), l*cos(X(:,1)),X(:,3)];
Y_cart = [l*sin(Y(:,1)), l*cos(Y(:,1)),Y(:,3)];

% Plot Observer and feature trajectory in 3-space
figure(12)
clf
plot3(X_cart(:,1),X_cart(:,3),-X_cart(:,2),'b',Y_cart(:,1),Y_cart(:,3),-Y_cart(:,2),'r:')
hold on
plot3(0,0,0,'b*',0,2,0,'r*')
legend('Observer', 'Dynamic Feature','fontsize', 18, 'interpreter', 'latex')
title('Plot of Pendulum Trajectories in Cartesian Coordinates','fontsize', 20, 'interpreter', 'latex')
xlabel('X position(m)','fontsize', 18, 'interpreter', 'latex')
ylabel('Z position (m)','fontsize', 18, 'interpreter', 'latex')
zlabel('Y position (m)','fontsize', 18, 'interpreter', 'latex')
xlim([-1 1])
ylim([-0.2 2])
zlim([-1 0.1])
grid on
%% simulate noisy sensor data

omega_sensed_x = X(:,2)+(omega_sigma*randn(n,1)); %observer omega sensed
theta_sensed_x = X(:,1)+(theta_sigma*randn(n,1)); %observer theta sensed

z_sensed = [Y(:,1) - X(:,1), Y(:,3)- X(:,3)]+(z_sigma*randn(n,2));

% figure(3)
% clf
% plot(T,z_sensed(:,1),'b',T,omega_sensed_x,'r')
% legend('Delta theta','Omega(t)')
% title('Sensor Values')
% xlabel('time')
% ylabel('radians')

%% Kalman Filter
%Basic Kalman filter with theta and omega measurement

f =@(X)X(ii-1,1) + X(ii-1,2)*dt; %motion model
h =@(theta_sensed_x)theta_sensed_x; %measurement model

Q = omega_sigma; 
R = theta_sigma;
X_noise(1) = x_0(1); %noisy motion model estimate
theta_est(1) = x_0(1); %noisy position measuremnt

for ii=2:n
    X_noise(ii) = X_noise(ii-1) + omega_sensed_x(ii-1)*dt;
    %z = h(X);
    z = theta_sensed_x(ii);
    [theta_est(ii),S]=ekf(f,X,S,h,z,Q,R);
%   mu = f(mu)+Q*randn();

    Sigma(ii) = S;
end

E = [(X(:,1)-X_noise.').^2,(X(:,1)-theta_est.').^2];
error_values = sum(E)

%Plot EKF results
figure(4)
clf
plot(T,theta_est,'b-.',T,X_noise,'k--',T,X(:,1),'r')
legend('$\hat{\theta}$','$\theta$ Noisy','$\theta$ Analytical','fontsize', 18, 'interpreter', 'latex')
title('Extended Kalman Filter Results','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('radians','fontsize', 18, 'interpreter', 'latex')

%Plot EKF error
figure(5)
clf
plot(T,E(:,1),'b',T,E(:,2),'r--')
legend('Process $\tilde\theta^2$','Kalman Filter $\tilde\theta^2$','fontsize', 18, 'interpreter', 'latex')
title('Kalman Filter vs. Process Noise $\tilde\theta^2$ Comparison','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('$radians^2$','fontsize', 18, 'interpreter', 'latex')


%% Baseline feature dynamics Regression

theta_sensed_y = Y(:,1) + theta_sigma*randn(n,1); %feature theta sensed
omega_sensed_y = Y(:,2) + omega_sigma*randn(n,1); %Feature omega sensed

B_0 = [1; 1; 1; 1]; %Initial coeefiecient guess
B_omega = gauss_newton_sin(T,omega_sensed_y,B_0); %sinusoidal omega regression
y_omega=-B_omega(1)*sin(B_omega(2)*Y(:,2)+B_omega(3))+B_omega(4); % evaluate regressed omega

B_theta = polyfit(T.',theta_sensed_y,poly_order); %polynomial regression for theta
y_theta = polyval(B_theta,T); %evaluate regressed theta

%plot baseline omega regression
figure(6)
clf
plot(T,Y(:,2),'r',T,y_omega,'--b',T,omega_sensed_y,'k-.')
legend('$\omega(t)$','$\omega$ Regression','$\omega$ Noisy','fontsize', 18, 'interpreter', 'latex')
title('Baseline Feature $\omega$ Regression','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('radians/s','fontsize', 18, 'interpreter', 'latex')
ylim([-8 10])

%plot baseline theta regression
figure(7)
clf
plot(T,Y(:,1),'r',T,y_theta,'--b',T,theta_sensed_y,'k-.')
legend('$\theta(t)$','$\theta$ Regression','$\theta$ Noisy','fontsize', 18, 'interpreter', 'latex')
title('Baseline Feature $\theta$ Regression','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('radians','fontsize', 18, 'interpreter', 'latex')
%% DLL Recurrursion with Kalman Filter and Nonlinear Regression

delta_theta = Y(:,1)-theta_est.' + z_sigma*randn(n,1); %noisy measurement
theta_est_y = theta_est.'+delta_theta; %noisy global reconstruction
 
p = polyfit(T.',theta_est_y,poly_order); %regression on theta
y_theta_DLL = polyval(p,T); %evaluate

figure(8)
clf
plot(T,Y(:,1),'r',T,y_theta_DLL,'--b',T,theta_est_y,'-.k')
legend('$\theta(t)$','$\hat\theta$ Regression','$\hat\theta$ Noisy','fontsize', 18, 'interpreter', 'latex')
title('DLL Feature $\theta$ Regression','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('radians','fontsize', 18, 'interpreter', 'latex')

omega_est_y(1)=0;
for ii=2:n %numerical differentiation
    omega_est_y(ii) = (theta_est_y(ii)-theta_est_y(ii-1))/dt;
end

%omega regression
B_0 = [1; 1; 1; 1];
B_omega_DLL = gauss_newton_sin(T,omega_sensed_y,B_0); %regression on omega
y_omega_DLL=-B_omega_DLL(1)*sin(B_omega_DLL(2)*Y(:,2)+B_omega_DLL(3))+B_omega_DLL(4);%evaluate

figure(9)
clf
plot(T,Y(:,2),'r',T,y_omega_DLL,'--b',T,omega_est_y,'-.k')
legend('$\omega(t)$','$\hat\omega$ Regression','$\hat\omega$ Noisy','fontsize', 18, 'interpreter', 'latex')
title('DLL Feature $\omega$ Regression','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('radians/s','fontsize', 18, 'interpreter', 'latex')

%Calculate DLL error
DLL_theta_error = [(Y(:,1)-y_theta.').^2,(Y(:,1)-y_theta_DLL.').^2];
DLL_theta_error_sum = sum(DLL_theta_error)

figure(10)
clf
plot(T,DLL_theta_error(:,1),'b',T,DLL_theta_error(:,2),'--r')
legend('Baseline $\tilde\theta^2$ Regression','DLL $\tilde\theta^2$ Regression','fontsize', 18, 'interpreter', 'latex')
title('DLL $\tilde\theta^2$ Regression vs. Baseline $\tilde\theta^2$ Regression Comparison','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('$radians^2$','fontsize', 18, 'interpreter', 'latex')

DLL_omega_error = [(Y(:,2)-y_omega).^2,(Y(:,2)-y_omega_DLL).^2];
DLL_omega_error_sum = sum(DLL_omega_error)

figure(11)
clf
plot(T,DLL_omega_error(:,1),'*',T,DLL_omega_error(:,2),'--')
legend('Baseline $\tilde\omega^2$ Regression','DLL $\tilde\omega^2$ Regression','fontsize', 18, 'interpreter', 'latex')
title('DLL $\tilde\omega^2$ Regression vs. Baseline $\tilde\omega^2$ Regression Comparison','fontsize', 20, 'interpreter', 'latex')
xlabel('time $(s)$','fontsize', 18, 'interpreter', 'latex')
ylabel('$(radians/s)^2$','fontsize', 18, 'interpreter', 'latex')


