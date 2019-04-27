function [x,P] = Kalman_filter(xi,eta,T,sigma_r,sigma_t,b,z)
    % White Noise Accleration Model
    F = [1, T, 0, 0;
         0, 1, 0, 0;
         0, 0, 1, T;
         0, 0, 0, 1];
     
    % Process noise Gain
    G = [T^2/2, 0;
         T,   0;
         0, T^2/2;
         0,   T];
    % Measurement matrix
    H = [1, 0, 0, 0;
         0, 0, 1, 0];
    
    % Process Noise Covariance
    Q = [T^4/4, T^3/2, 0,    0;
         T^3/2, T^2,   0,    0;
          0,   0,  T^4/4, T^3/2;
          0,   0,  T^3/2,  T^2];
    
    % Measurement Covariance   
    syms r theta;   
    r11 = (b^-2 - 2)*r^2*(cos(theta))^2 + (1/2)*(r^2 + sigma_r )*(1 + b^4*cos(2*theta));   
    r22 = (b^-2 - 2)*r^2*(sin(theta))^2 + (1/2)*(r^2 + sigma_r )*(1 - b^4*cos(2*theta));
    r12 = (((b^-2)*r^2)/2 + (r^2 + sigma_r )*(b^4)/2 - r^2)*sin(2*theta);
    
    % Measurement Covariance
    R_cov = [r11,r12;
             r12,r22];
    R = matlabFunction(R_cov);
      
    % State Covariance Initialization
    p_0 = double(subs(r11, [r,theta], [z(1,2),z(2,2)]));
    p_1 = double(subs(r22, [r,theta], [z(1,2),z(2,2)]));
    P(:,:,1) = [p_0,   p_0/T,      0,     0;
                p_0/T, p_0/(T^2),  0,     0;
                0,       0,     p_1,   p_1/T;
                0,       0,    p_1/T, p_1/(T^2)];
    
    % State Initialization
    x = zeros(4,size(xi,2)-2);
    x(1,1) = xi(2); x(3,1) = eta(2);
    x(2,1) = (xi(2)-xi(1))/T;
    x(4,1) = (eta(2)-eta(1))/T;
    
    % Run Kalman Filter 
    for i = 2:size(xi,2)-1
        % State Covariance and Measurement prediction Prediction
        x(:,i) = F*x(:,i-1) + G*(randn(2,1)*sqrt(4));
        P(:,:,i) = F*P(:,:,i-1)*F' + Q.*4;
        z_pred = H*x(:,i);
        
        % Measurement Residual, Innovation Covariance, Filter Gain
        nu = [xi(i+1); eta(i+1)] - z_pred;
        S = feval(R,z(1,i+1),z(2,i+1)) + H*P(:,:,i)*H';
        W = P(:,:,i)*H'*inv(S);
        
        % State and Covariance update
        x(:,i) = x(:,i) + W*nu;
        P(:,:,i) = P(:,:,i) - W*S*W';
    end
end
