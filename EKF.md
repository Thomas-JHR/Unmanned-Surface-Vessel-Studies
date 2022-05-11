# Linear Extended Kalman Filter

 - This function describes filtration and poeses estimation in vessel dynamics and kinematics systems. 
 
 - Nonlinear observers and estimators have to be initially modeled such that the vessel model is observable and controllable. 
 
 - The discrete-time Extended Kalman filter for the multivariable system of vessel motion control is designed based on the **Nonlinear model**  estimate the poses
 
First of all, let's define the EKF-localization states as follows:

<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/ekf1.svg>
</p><p align="center">

```f(q_k,τ_k)``` is a nonlinear vector field. ```ν_k``` and ```w_k``` are white input noise processes and measurement white senor noise, respectively. The filter system involves the process noise covariance weight matrix ```Q_EKF``` and pose measurement noise covariance matrix ```R_EKF``` described as
 
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/variance_ekf.svg>
</p><p align="center">

The basic extended Kalman filter algorithm including the Priori estimate and the Posterior estimate is shown below to estimate the measured poses of the vessel. Initial values of estimated poses and covariance matrixes are predicted before the progress of filtering:
 
 Prediction Poriori Estimate:
 
 <p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/ekf_11.svg>
</p><p align="center">
 
 
 Correction Posteriori Estimate:
 
<p align="center">
<img src= https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/ekf_12.svg>
</p><p align="center">
 
```F_k```, ```L_k```, and ```M_k``` are Jacobian matrixes determined using the Euler forward method.
 
  <p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/ekf_13.svg>
</p><p align="center">
 
 
 ```matlab
 
 function [eta_new,P_new] = Sensor_EKF(D,M,X_u_dot,Y_v_dot,t,eta,u_m,z,R,Q,P)

A = [0 0 0 cos(eta(3)) -sin(eta(3)) 0;...
     0 0 0 sin(eta(3)) cos(eta(3)) 0;...
     0 0 0 0 0 1;...
     0 0 0 -D(1,1)/M(1,1) M(2,2)*eta(6)/M(1,1) 0;...
     0 0 0 0 -D(2,2)/M(2,2) -M(1,1)*eta(4)/M(2,2);...
     0 0 0 -(X_u_dot-Y_v_dot)*eta(5)/M(3,3) 0 -D(3,3)/M(3,3)];
 
B = [0 0 0;0 0 0; 0 0 0;1/M(1,1) 0 0;0 1/M(2,2) 0;0 0 1/M(3,3)]; 
% Define system matrices

% C = [1 0 0 0 0 0;...
%     0 1 0 0 0 0;...
%     0 0 1 0 0 0];
C = eye(6);

F = [1 0 -t*(eta(4)*sin(eta(3))+eta(5)*cos(eta(3))) t*cos(eta(3)) -t*sin(eta(3)) 0;...
     0 1 t*(eta(4)*cos(eta(3))-eta(5)*sin(eta(3))) t*sin(eta(3)) t*cos(eta(3)) 0;...
     0 0 1 0 0 t;...
     0 0 0 1-D(1,1)*(t/M(1,1)) M(2,2)*t*eta(6)/M(1,1) M(2,2)*t*eta(5)/M(1,1);...
     0 0 0 -M(1,1)*eta(6)*(t/M(2,2)) 1-D(2,2)*(t/M(2,2)) -M(1,1)*t*eta(4)/M(2,2);...
     0 0 0 -t*(X_u_dot-Y_v_dot)*eta(5)/M(3,3) -t*(X_u_dot-Y_v_dot)*eta(4)/M(3,3) 1-D(3,3)*(t/M(3,3))];
L = t*B;
% Compute priori estimate
eta_new = eta + t*(B*u_m + A*eta); % use our nonlinear model
P_new = F*P*transpose(F) + L*Q*transpose(L);
P_new = (P_new + transpose(P_new))/2;

% % Define more system matrices
%  
H = C;
M_matrix = eye(6);

% Compute posteriori estimate
K = P_new * transpose(H) * inv(H * P_new * transpose(H) + M_matrix * R * transpose(M_matrix)); 
eta_new = eta_new + K*(z - C*eta_new);
P_new = (eye(6) - K*H)*P_new;
P_new = (P_new + transpose(P_new))/2;

end
 
 ```
 
