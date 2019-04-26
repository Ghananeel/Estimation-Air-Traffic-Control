% Create trajectory for the Aircraft 
x = trajectory();
xi = x(1,:); xi_d = x(2,:); eta = x(3,:); eta_d = x(4,:);
plot(xi, eta, 'LineWidth', 2); title('Trajectory of the Aircraft');
xlabel({'$\xi$'},'Interpreter','latex','FontSize',14); 
ylabel({'$\eta$'},'Interpreter','latex','FontSize',14); 

% Measurement indices
k = 1:10:500;

% Create Sensor measurements
xi0 = -1e+4; eta0 = 0;
r = sqrt((xi(k)-xi0).^2 + (eta(k)-eta0).^2);
theta = atan((eta(k)-eta0)./(xi(k)-xi0));
z = [r + randn(size(k)).*sqrt(2500); theta + randn(size(k)).*sqrt(1)];
figure; plot(theta); hold on
plot(z(2,:))

