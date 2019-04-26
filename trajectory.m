function x = trajectory()
    dt = 1;
    % Starting Const Vel State
    state = [0;0;0;250];
    
    % Create each maneuver from start time to end
    x1 = const_vel(1,100,state,dt);
    x2 = ct(2, 101, 130, x1(:,end), dt);
    x3 = const_vel(131 ,200 ,x2(:,end),dt);
    x4 = ct(-1, 201, 245, x3(:,end), dt);
    x5 = ct(1, 246, 335, x4(:,end), dt);
    x6 = ct(-1, 336, 380, x5(:,end), dt);
    x7 = const_vel(381 ,500 ,x6(:,end),dt);
    
    % Concatenate all maneuvers
    x = [x1,x2,x3,x4,x5,x6,x7];
 end

function x = const_vel(t_start, t_end, x0, dt)
    x = []; x(:,1) = x0;
    % Discrete State Matrix for Const Vel Model
    A = [1, dt, 0, 0;
          0, 1, 0, 0;
          0, 0, 1, dt;
          0, 0, 0, 1];
    % Run for total time
    for i = 1:(t_end - t_start)
        x(:,i+1) = A*x(:,i);
    end
end
    
function x = ct(omega, t_start, t_end, x0, dt)
    omega = deg2rad(omega);
    % Cont State Matrix for a CT Model 
    f = @(t,x)[x(2); -omega*x(4); x(4); omega*x(2)];
    % Runge-Kutta integration
    [t,x] = ode45(f,t_start:dt:t_end,x0);
    x = transpose(x);
end

