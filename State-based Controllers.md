# Desired Path Design

The path is divided into five segments under fixed simulation time 63 seconds, and sampled at 100 Hz. The path starts at a 3.15 m straight-line section with a constant surge speed at 0.5 m/s, followed by turning into a three-fourth fraction of a circle that has a radius curvature of 3.15 m at the constant surge speed 0.785 m/s and angular rate 0.2494 rad/s. 

The third segment is a 6.3 m vertical line section at constant surge speed 0.5 m/s, while the fourth segment is another three-fourth fraction of a circle similar to the physical parameters from first circle but negative angular rate. 

Finally, the path is composed of another 3.15 m straight line with the same surge speed at 0.5 m/s. An additional section piece may be needed for the MPC scheme to look ahead future desired states, and this extra segment revert back to the first segment of the infinite loop until ```p``` horizon points are met.

<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Images/Desired%20path.jpg width=60% height=60%>
</p><p align="center">
  
# PD+FBL
 
There are some extra equations used to help us to figure out the inputs and outputs of the controller and transformation of the states from controller to the USV:
  
Let's take derivatives of the vessel kinematic equation to get the acceleration vector in the global frame

<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PD1.svg>
</p><p align="center">

recover ```z``` specified as the FBL state vector ```eta``` the FBL input vector equal to the acceleration vector in the global frame
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/12.svg>
</p><p align="center">
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/equation_FBL.svg>
</p><p align="center">
  
Using basic feedback control law ```u=-Kx```, we have
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PD2.svg>
</p><p align="center">
  
Since the desired accleration of the routes are assumed to be zero, let's rewrite the equation:
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PD3.svg>
</p><p align="center">  

The next step is to transform the acceleration vector in the global frame to the acceleration in the body fixed frame, followed by inverse dynamics to get the resultant force and moment vector.
  
 ```matlab
%% Controller part initialization  
q_FSF = zeros(6,N); % x y psi u v r
tao_FSF = zeros(3,N);
f_FSF = zeros(2,N);
q_FSF(1:3,1) = [0.5 1 pi/6]';% Set up initial conditions
q_dot_FSF = zeros(6,N);
a_n = zeros(3,N);
a_b = zeros(3,N);
G = zeros(6,3);
G(4,1) = 1/M(1,1);
G(5,2) = 1/M(2,2);
G(6,3) = 1/M(3,3);

%% Estimator part initialization
q_hat_FSF = zeros(6,N);
aa = [0.3 -0.3 0.1 0.1 0.1 0.1]; 
P_hat_FSF = zeros(6,6,N); % Define covariance matrix
P_hat_FSF(:,:,1) = P_guess;
q_hat_FSF(:,1) = q_guess + transpose(aa); % set the initial guess for estimator
  
tic 
for i = 2:N
    %% simulate USV
    
    [q_dot_FSF(1:3,i-1),q_dot_FSF(4:6,i-1)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_FSF(4:6,i-1),tao_FSF(:,i-1),q_FSF(1:3,i-1),t);
    [q_FSF(1:3,i),q_FSF(4:6,i)]= Vehicle(q_dot_FSF(1:3,i-1),q_dot_FSF(4:6,i-1),q_FSF(1:3,i-1),q_FSF(4:6,i-1),t);
    
    %% Model propriocdeptive sensors with processing noise
    u_m(1,i-1) = tao_FSF(1,i-1) + sqrt(Q_EKF(1,1))*randn(1); 
    u_m(2,i-1) = tao_FSF(2,i-1) + sqrt(Q_EKF(2,2))*randn(1);
    u_m(3,i-1) = tao_FSF(3,i-1) + sqrt(Q_EKF(3,3))*randn(1);
        
    %% Model GPS,COMPASS and IMU with noise    
    z_m(1,i) = q_FSF(1,i) + sqrt(R_EKF(1,1))*randn(1);
    z_m(2,i) = q_FSF(2,i) + sqrt(R_EKF(2,2))*randn(1);
    z_m(3,i) = q_FSF(3,i) + sqrt(R_EKF(3,3))*randn(1);
    z_m(4,i) = q_FSF(4,i) + sqrt(R_EKF(4,4))*randn(1);
    z_m(5,i) = q_FSF(5,i) + sqrt(R_EKF(5,5))*randn(1);
    z_m(6,i) = q_FSF(6,i) + sqrt(R_EKF(6,6))*randn(1);
    
    %% Simulate EKF fusion of GPS 
    [q_hat_FSF(:,i),P_hat_FSF(:,:,i)] = Sensor_EKF(D,M,X_u_dot,Y_v_dot,t,q_hat_FSF(:,i-1),u_m(:,i-1),z_m(:,i),R_EKF,Q_EKF,P_hat_FSF(:,:,i-1));
   
    % define g_0
    
    g_0(1,4) = cos(q_hat_FSF(3,i)); 
    g_0(1,5) = -sin(q_hat_FSF(3,i)); 
    g_0(2,4) = sin(q_hat_FSF(3,i)); 
    g_0(2,5) = cos(q_hat_FSF(3,i)); 
    g_0(3,6) = 1;
    g_0(4,4) = -D(1,1)/M(1,1);
    g_0(4,5) = M(2,2)*q_hat_FSF(6,i)/M(1,1);
    g_0(5,5) = -D(2,2)/M(2,2);
    g_0(5,6) = M(1,1) * q_hat_FSF(4,i)/M(2,2);
    g_0(6,4) = (M(1,1)-M(2,2)) * q_hat_FSF(5,i) / M(3,3);
    g_0(6,6) = -D(3,3) / M(3,3);
    
    %% I/O linearization PD Control
    
    err_FSF(1,i) = q_hat_FSF(1,i-1) - q_d(1,i);
    err_FSF(2,i) = q_dot_FSF(1,i-1) - q_d_dot(1,i);
    err_FSF(3,i) = q_hat_FSF(2,i-1) - q_d(2,i);
    err_FSF(4,i) = q_dot_FSF(2,i-1) - q_d_dot(2,i);
    err_FSF(5,i) = q_hat_FSF(3,i-1) - q_d(3,i);
    err_FSF(6,i) = q_dot_FSF(3,i-1) - q_d_dot(3,i);
    Coriolis_FSF = [0 0 -M(2,2)*q_hat_FSF(5,i);
                0 0 M(1,1)*q_hat_FSF(4,i);
                M(2,2)*q_hat_FSF(5,i) -M(1,1)*q_hat_FSF(4,i) 0];
    a_n(:,i) = -K_FSF_MATRIX * err_FSF(:,i);
    J_dot = [-q_hat_FSF(6,i)*sin(q_hat_FSF(3,i)) -q_hat_FSF(6,i)*cos(q_hat_FSF(3,i)) 0;
             q_hat_FSF(6,i)*cos(q_hat_FSF(3,i)) -q_hat_FSF(6,i)*sin(q_hat_FSF(3,i)) 0;
             0 0 0];
    a_b(:,i) = g_0(1:3,4:6) \ (a_n(:,i) - J_dot * q_hat_FSF(4:6,i));
    tao_FSF(:,i) = M * a_b(:,i) + Coriolis_FSF * q_hat_FSF(4:6,i) + D * q_hat_FSF(4:6,i);
    f_FSF(1,i) = tao_FSF(1,i)*0.5 + tao_FSF(3,i)/B;
    f_FSF(2,i) = tao_FSF(1,i)*0.5 - tao_FSF(3,i)/B;
    
    %% Simulate external disturbance that exerts to the body
    tao_wave(1,i) = ExwaveForce(PSD_wave(1,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(1,i));
    tao_wave(2,i) = ExwaveForce(PSD_wave(2,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(2,i));
    tao_wave(3,i) = ExwaveForce(PSD_wave(3,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(3,i));
    tao_wind(:,i) = ExwindForce(beta_w(1,i),V_w(1,i),q_FSF(4:6,i),q_FSF(3,i),rho_air,L);
    %tao_FSF(:,i) = tao_FSF(:,i) - tao_wind(:,i) - tao_wave(:,i);    
end
toc
```
  
# Nonlinear PID
  
A PID Proportional-Integration-Derivative (PID) controller with the approximate linearization technique is selected as the second feedback control algorithm to minimize the error and operate towards the desired states by monitoring the forces and moment as the control input.
  
```q_d``` is defined as the desired pose of Heron M300 in the global frame and the desired velocity state in the body-fixed frame, respectively. The PID controller is extended with an extra vector
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PID1.svg>
</p><p align="center">
  
The term is expressed as an integral term control vector of the PID controller. To find the integral term, the global frame positions are integrated as a function of the time step and the position vectors, and ```T``` is the current step of the simulation.
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PID2.svg>
</p><p align="center">

Let ```q_t``` redefined to be an error vector between the actual pose and the desired pose in the global frame, and ```V_t``` be an error vector between actual state and the desired state of velocity in the body-fixed frame shown below. ```τ_d``` is the desired input force and resulting moment vector and measured by desired waypoints simulation. A nonlinear PID control law depending on the state feedback control law is defined as:    
  
 <p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PID3.svg>
</p><p align="center">

```K_I```,```K_D```,and ```K_P``` are expressed as a positive integration gain matrix, derivative gain matrix, and proportional gain matrix for each vessel vectors. The gain matrixes values are stabilized to pursue the equilibrium point at each time step by placing left-hand side poles for each vessel state. Let the error state X derive as 

<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PID4.svg>
</p><p align="center">
  
More transformation matrixes are defined below to express the error state
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PID5.svg>
</p><p align="center">
  
A simplified version of the continuous error state space equation is derived as follows.
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PID6.svg>
</p><p align="center">  
  
all ```X``` states are modified to be ```X = X_hat + X_d``` before taking Jocobian partial derivatives with respect to ```X_hat```, followed by putting the equilibrium point (0,0) into the stability matrix. Pole vector should be a 1×9 matrix, and output of the function pole placement is a PID gain 3×9 matrix
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/PID7.svg>
</p><p align="center"> 
  
```matlab
%% Control part initialization
q_PID = zeros(6,N); % x y psi u v r
tao_PID = zeros(3,N);
f_PID = zeros(2,N);
q_PID(1:3,1) = [0.5 1 pi/6]';% Set up initial conditions
q_dot_PID = zeros(6,N);
K_gain = zeros(3,9,N);
uI = zeros(3,N);
Integral_term = zeros(3,N);

%% Estimator part initialization
q_hat_PID = zeros(6,N);
aa = [0.3 -0.3 0.1 0.1 0.1 0.1]; 
q_guess = q_PID(:,1);

P_guess = diag([5 -5 1 1 1 1].^2);  % Set the initial pose covariance estimate as a diagonal matrix

q_hat_PID(:,1) = q_guess + transpose(aa); % set the initial guess for estimator

z_m = zeros(6,N);% measurement odemetry and GPS signal, as well as IMU unit
u_m = zeros(3,N); % input value with noise

R_EKF = diag([0.5 0.5 0.5 0.05 0.05 0.05].^2);% Set the true process and measurement noise covariances (output measurement noise)
Q_EKF = diag([0.1 0.1 0.1].^2); %(input processing noise)

P_hat_PID = zeros(6,6,N); % Define covariance matrix
P_hat_PID(:,:,1) = P_guess;
  
tic
for i = 2:N
    %% Simulate PID Controller
    [K_gain(:,:,i-1), tao_PID(:,i-1),Integral_term(:,i-1)] = FSF_PID(q_hat_PID(:,i-1), q_d(:,i), M, D, X_u_dot, Y_v_dot, tao_d(:,i), t, uI(:,i-1), i);
    tao_wave(1,i-1) = ExwaveForce(PSD_wave(1,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_PID(:,i-1),Drift_wave(1,i));
    tao_wave(2,i-1) = ExwaveForce(PSD_wave(2,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_PID(:,i-1),Drift_wave(2,i));
    tao_wave(3,i-1) = ExwaveForce(PSD_wave(3,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_PID(:,i-1),Drift_wave(3,i));
    tao_wind(:,i-1) = ExwindForce(beta_w(1,i),V_w(1,i),q_FSF(4:6,i-1),q_PID(3,i-1),rho_air,L);   
%    [tao_PID(1,i-1),tao_PID(3,i-1)] = sat(tao_PID(1,i-1),tao_PID(3,i-1));
    

    f_PID(1,i-1) = tao_PID(1,i-1)*0.5 + tao_PID(3,i-1)/B;
    f_PID(2,i-1) = tao_PID(1,i-1)*0.5 - tao_PID(3,i-1)/B;
   % tao_PID(:,i-1) = tao_PID(:,i-1) - tao_wind(:,i-1) - tao_wave(:,i-1);
    [q_dot_PID(1:3,i-1),q_dot_PID(4:6,i-1)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_PID(4:6,i-1),tao_PID(:,i-1),q_PID(1:3,i-1),t);
    [q_PID(1:3,i),q_PID(4:6,i)]= Vehicle(q_dot_PID(1:3,i-1),q_dot_PID(4:6,i-1),q_PID(1:3,i-1),q_PID(4:6,i-1), t);    
    
    %% Model propriocdeptive sensors with processing noise
    u_m(1,i-1) = tao_PID(1,i-1) + sqrt(Q_EKF(1,1))*randn(1); 
    u_m(2,i-1) = tao_PID(2,i-1) + sqrt(Q_EKF(2,2))*randn(1);
    u_m(3,i-1) = tao_PID(3,i-1) + sqrt(Q_EKF(3,3))*randn(1);
        
    %% Model GPS,COMPASS and IMU with noise    
    z_m(1,i) = q_PID(1,i) + sqrt(R_EKF(1,1))*randn(1);
    z_m(2,i) = q_PID(2,i) + sqrt(R_EKF(2,2))*randn(1);
    z_m(3,i) = q_PID(3,i) + sqrt(R_EKF(3,3))*randn(1);
    z_m(4,i) = q_PID(4,i) + sqrt(R_EKF(4,4))*randn(1);
    z_m(5,i) = q_PID(5,i) + sqrt(R_EKF(5,5))*randn(1);
    z_m(6,i) = q_PID(6,i) + sqrt(R_EKF(6,6))*randn(1);
    
    %% Simulate EKF fusion of GPS 
    [q_hat_PID(:,i),P_hat_PID(:,:,i)] = Sensor_EKF(D,M,X_u_dot,Y_v_dot,t,q_hat_PID(:,i-1),u_m(:,i-1),z_m(:,i),R_EKF,Q_EKF,P_hat_PID(:,:,i-1));
end
toc
  
  function [K_gain, new_tao, Integral_term] = FSF_PID(q_b, q_d, M, D, X_u_dot, Y_v_dot, tao_ref, t, Integral_term_a, i)

    %poleX = [-13.25 -13.11 -25.12 -34.2 -34.25 -34.1 -14.8 -8.6 -10.43]; % Define our pole
    %poleX = [-9.25 -10.11 -9.5 -27.2 -30 -28.5 -9.8 -10.1 -7.15];
    poleX = [-9.25 -10.11 -9.5 -24.2 -30 -25.25 -5.8 -5.1 -6.97];
    
    %poleX = [-3.5 -2.9 -0.5 -10.11 -15 -20.25 -3 -3.25 -0.3];
    %poleX = [-3.25 -2.85 -1.7 -80.2 -90 -70.25 -2.9 -2.5 -1.7];
    A_Matrix = zeros(9);
    B_Matrix = cat(1,zeros(6,3),M);
    
    A_Matrix(1,2) = q_d(6);
    A_Matrix(1,4) = 1;
  

    A_Matrix(2,1) = -q_d(6);
    A_Matrix(2,5) = 1;

    A_Matrix(3,6) = 1;

    A_Matrix(4,5) = q_d(6);    
    A_Matrix(4,7) = 1;  
    
    A_Matrix(5,4) = -q_d(6);
    A_Matrix(5,8) = 1;

    A_Matrix(6,9) = 1;

    A_Matrix(7,7) = -D(1,1)*q_d(4)/M(1,1);
    A_Matrix(7,8) = (M(2,2)/M(1,1))*q_d(5)*q_d(6);
    A_Matrix(7,9) = (M(2,2)/M(1,1))*q_d(6)*q_d(5);

    A_Matrix(8,7) = -(M(1,1)/M(2,2))*q_d(6)*q_d(4);
    A_Matrix(8,8) = -D(2,2)*q_d(5)/M(2,2);
    A_Matrix(8,9) = -(M(1,1)/M(2,2))*q_d(4)*q_d(6);

    A_Matrix(9,7) = -((X_u_dot - Y_v_dot)/M(3,3))*q_d(4)*q_d(5);
    A_Matrix(9,8) = -((X_u_dot - Y_v_dot)/M(3,3))*q_d(4)*q_d(5);
    A_Matrix(9,9) = -D(3,3)*q_d(6)/M(3,3);
    
    %% Single place pole design
    T = 0:t:6301;
    K_gain = place(A_Matrix,B_Matrix,poleX);
    %% compute PID Integral Term   

    Integral_term_a(1,i-1) = (cos(q_b(3))*(q_b(1) - q_d(1)) + sin(q_b(3))*(q_b(2) - q_d(2)));
    Integral_term_a(2,i-1) = (-sin(q_b(3))*(q_b(1) - q_d(1)) + cos(q_b(3))*(q_b(2) - q_d(2)));
    Integral_term_a(3,i-1) = (q_b(3) - q_d(3));
    Integral_term = trapz(T(1,1:i-1),Integral_term_a(:,1:i-1),2);
    
    %% apply PID control law
    new_tao(1) = -K_gain(1,4) * (cos(q_b(3))*(q_b(1) - q_d(1)) + sin(q_b(3))*(q_b(2) - q_d(2))) - K_gain(1,7) * (q_b(4)-q_d(4)) - K_gain(1,1)*Integral_term(1) + tao_ref(1); 
    new_tao(2) = -K_gain(2,5) * (-sin(q_b(3))*(q_b(1)-q_d(1)) + cos(q_b(3))*(q_b(2) - q_d(2))) - K_gain(2,8) * (q_b(5)-q_d(5)) - K_gain(2,2)*Integral_term(2) + tao_ref(2); 
    new_tao(3) = -K_gain(3,6) * (q_b(3) - q_d(3)) - K_gain(3,9) * (q_b(6) - q_d(6)) - K_gain(3,3)*Integral_term(3) + tao_ref(3);
    
end
```
# MPC+FBL
# Nonlinear MPC
# HMPC+FBL
# GPR+PD+FBL
# GPR+MPC+FBL

