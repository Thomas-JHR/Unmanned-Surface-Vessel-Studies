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
# MPC+FBL
# Nonlinear MPC
# HMPC+FBL
# GPR+PD+FBL
# GPR+MPC+FBL

