%% Create trajectory for the Aircraft 
x = trajectory();
xi = x(1,:); xi_d = x(2,:); eta = x(3,:); eta_d = x(4,:);

% Measurement indices
T = 10; k = 1:T:501;
% Measurement Noise variances
sigma_t = 1; sigma_r = 2500;

%% Create Sensor measurements
% Sensor Location
xi0 = -1e+4; eta0 = 0;
% Range measurement
r = sqrt((xi(k)-xi0).^2 + (eta(k)-eta0).^2);
% Angle measurement
theta = atan2((eta(k)-eta0),(xi(k)-xi0));
% Measurement + noise --> z(k) = h(x(k))+w(k)
z = [r + randn(size(k)).*sqrt(sigma_r); theta + randn(size(k)).*deg2rad(sqrt(sigma_t))];

%% Unbiased Conversion of Sensor Measurements
% Multiplicative Bias
b = exp(-(deg2rad(sigma_t)^2)/2);
% Translate trajectory in xi dirn as sensor is not at [0,0]
xi_m = b^-1.*z(1,:).*cos(z(2,:))-1e+4;
% Compute trajectory in eta dirn
eta_m = b^-1.*z(1,:).*sin(z(2,:));

%% Filtering
% Run Kalman Filter
[x_kf, P] = Kalman_filter(xi_m,eta_m,T,sigma_r,(deg2rad(sigma_t))^2,b,z);

% Run IMM-L

% Run IMM-CT

%% Plot trajectory and metrics
figure
plot(xi, eta, 'LineWidth', 2); hold on;
plot(xi_m,eta_m,'-*');
plot(x_kf(1,:),x_kf(3,:),'-.o','LineWidth',1.25);
title('Trajectory of the Aircraft');
xlabel({'$\xi$'},'Interpreter','latex','FontSize',14); 
ylabel({'$\eta$'},'Interpreter','latex','FontSize',14); 