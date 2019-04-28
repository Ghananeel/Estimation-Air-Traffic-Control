classdef KFilter
    properties
        theta, r
        
        % State Matrix
        F = [1, T, 0, 0;
             0, 1, 0, 0;
             0, 0, 1, T;
             0, 0, 0, 1];
         
        % Process Noise Gain
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
         
         r11 = (b^-2 - 2)*r^2*(cos(theta))^2 + (1/2)*(r^2 + sigma_r )*(1 + b^4*cos(2*theta));   
         r22 = (b^-2 - 2)*r^2*(sin(theta))^2 + (1/2)*(r^2 + sigma_r )*(1 - b^4*cos(2*theta));
         r12 = (((b^-2)*r^2)/2 + (r^2 + sigma_r )*(b^4)/2 - r^2)*sin(2*theta);
         
         % Measurement Covariance
         R_cov = [r11,r12;
                  r12,r22];
              
         % State Covariance Initialization
         P = [r11,   r11/T,      0,     0;
              r11/T, r11/(T^2),  0,     0;
                0,       0,     r22,   r22/T;
                0,       0,    r22/T, r22/(T^2)];
    end
    methods
        function 
