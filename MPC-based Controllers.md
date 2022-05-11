# Linear MPC+FBL

MPC-based controllers for the Heron vessel were developed based on the the controllers development of the hydraulic manipulator in [^1]. MPC is a methdology attrempts to optimize the system behaviour over same prediction horizon ```p``` as a discrete number of time steps in the future from the current state ```k```. Suppose that ```p``` is an array ```[1 2 ... p]``` and define a cost function as

  <p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/MPC1.svg>
</p><p align="center">

The vectors inside the cost function are defined as 
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/MPC2.svg>
</p><p align="center">  

```eta_vector``` is future input accleration vector in the global frame. ```z_vector``` is prediction output state vector. ```z``` is the desired prediction output state vector
  
Let's recover state-space equation for the nonlinear vessel. The linear system as a function of the FBL control state with respect to the FBL input acceleration vector, previously defined in the FBL section, is used to estimate the next FBL control state by generating sequential mathematical expression as a discrete-time linear system by discretizing the continuous-time linear state space model,
 
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/MPC3.svg>
</p><p align="center">    
  
The first discrete time step linear system is successfully derived as a function of the previous FBL control state and the previous FBL input state. To estimate the FBL control state at ```k+2```, we use lookhead (prediction horizon timesteps)  and compute the MPC-based controller design matrixces ```M``` and ```L```, and positive weighting state cost ```Q``` and input cost matrices ```R``` that are tunable to change path following behaviours
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/MPC4.svg>
</p><p align="center">      
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/MPC5.svg>
</p><p align="center">     
  
Next we are gonna to derive the cost function as
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/MPC6.svg>
</p><p align="center">    
  
Optimize ```eta_vector``` by taking partial derivatives wrt ```eta_vector``` and rearange the equation as
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/MPC7.svg>
</p><p align="center">    
  
``` matlab
  %% MPC variables define
F_MPC = eye(6) + t * A_FSF;
G_MPC = t * B_FSF; 
[nr,nc] = size(G_MPC);
q_MPC = zeros(nr,N);
q_MPC(1:3,1) = [0.5 1 pi/6]';
z_MPC = zeros(nr,N);

q_dot_MPC = zeros(nr,N);
tao_MPC = zeros(nc,N);
f_MPC = zeros(2,N);
eta_MPC = zeros(nc,N);

L_MPC = zeros(nr*p,nr);
M_MPC = zeros(nr*p,nc*p);
for i =1:p
    L_MPC(nr*i-nr+1:nr*i, 1:nr) = mpower(F_MPC, i);
    for j =1:p-i+1
        M_MPC(nr*(p-i+1)-nr+1:nr*(p-i+1), nc*j-nc+1:nc*j) = mpower(F_MPC, p-i-j+1) * G_MPC;
    end
end

Q_MPC = 1 * kron(eye(p),eye(nr)); % State Cost
R_MPC = 0.36 * kron(eye(p),eye(nc)); % Input Cost 

K_MPC = (transpose(M_MPC)*Q_MPC*M_MPC + R_MPC) \ transpose(M_MPC) * Q_MPC;

Zd_MPC = zeros(nr*p,1); % For only MPC use
Zd_NMPC = zeros(nr*p,1); % For only NMPC use

%% Estimator part initialization
q_hat_MPC = zeros(6,N);
aa = [0.3 -0.3 0.1 0.1 0.1 0.1]; 
P_guess = diag([5 -5 1 1 1 1].^2);  % Set the initial pose covariance estimate as a diagonal matrix
P_hat_MPC = zeros(6,6,N); % Define covariance matrix
P_hat_MPC(:,:,1) = P_guess;
q_hat_MPC(:,1) = q_guess + transpose(aa); % set the initial guess for estimator
  
tic
for i = 2:N
    [q_dot_MPC(1:3,i-1),q_dot_MPC(4:6,i-1)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_MPC(4:6,i-1),tao_MPC(:,i-1),q_MPC(1:3,i-1),t);
    [q_MPC(1:3,i),q_MPC(4:6,i)]= Vehicle(q_dot_MPC(1:3,i-1),q_dot_MPC(4:6,i-1),q_MPC(1:3,i-1),q_MPC(4:6,i-1),t);
    
    %% Model propriocdeptive sensors with processing noise
    u_m(1,i-1) = tao_MPC(1,i-1) + sqrt(Q_EKF(1,1))*randn(1); 
    u_m(2,i-1) = tao_MPC(2,i-1) + sqrt(Q_EKF(2,2))*randn(1);
    u_m(3,i-1) = tao_MPC(3,i-1) + sqrt(Q_EKF(3,3))*randn(1);
        
    %% Model GPS,COMPASS and IMU with noise    
    z_m(1,i) = q_MPC(1,i) + sqrt(R_EKF(1,1))*randn(1);
    z_m(2,i) = q_MPC(2,i) + sqrt(R_EKF(2,2))*randn(1);
    z_m(3,i) = q_MPC(3,i) + sqrt(R_EKF(3,3))*randn(1);
    z_m(4,i) = q_MPC(4,i) + sqrt(R_EKF(4,4))*randn(1);
    z_m(5,i) = q_MPC(5,i) + sqrt(R_EKF(5,5))*randn(1);
    z_m(6,i) = q_MPC(6,i) + sqrt(R_EKF(6,6))*randn(1);
    
    %% Simulate EKF fusion of GPS 
    [q_hat_MPC(:,i),P_hat_MPC(:,:,i)] = Sensor_EKF(D,M,X_u_dot,Y_v_dot,t,q_hat_MPC(:,i-1),u_m(:,i-1),z_m(:,i),R_EKF,Q_EKF,P_hat_MPC(:,:,i-1));
    
    %% MPC Control
    z_MPC(1,i-1) = q_hat_MPC(1,i-1);
    z_MPC(2,i-1) = q_dot_MPC(1,i-1);
    z_MPC(3,i-1) = q_hat_MPC(2,i-1);
    z_MPC(4,i-1) = q_dot_MPC(2,i-1);
    z_MPC(5,i-1) = q_hat_MPC(3,i-1);
    z_MPC(6,i-1) = q_dot_MPC(3,i-1);
    for j = 1:p
        Zd_MPC(nr*j-nr+1,1) = q_d_new(1,i+j-1);
        Zd_MPC(nr*j-nr+2,1) = q_d_dot_new(1,i+j-1);
        Zd_MPC(nr*j-nr+3,1) = q_d_new(2,i+j-1);
        Zd_MPC(nr*j-nr+4,1) = q_d_dot_new(2,i+j-1);
        Zd_MPC(nr*j-nr+5,1) = q_d_new(3,i+j-1);
        Zd_MPC(nr*j-nr+6,1) = q_d_dot_new(3,i+j-1);
    end
    eta = -K_MPC * (L_MPC*z_MPC(:,i-1) -Zd_MPC);
    eta_MPC(:,i) = eta(1:3);
    g_0(1,4) = cos(q_MPC(3,i)); 
    g_0(1,5) = -sin(q_MPC(3,i)); 
    g_0(2,4) = sin(q_MPC(3,i)); 
    g_0(2,5) = cos(q_MPC(3,i)); 
    g_0(3,6) = 1;
    g_0(4,4) = -D(1,1)/M(1,1);
    g_0(4,5) = M(2,2)*q_MPC(6,i)/M(1,1);
    g_0(5,5) = -D(2,2)/M(2,2);
    g_0(5,6) = M(1,1) * q_MPC(4,i)/M(2,2);
    g_0(6,4) = (M(1,1)-M(2,2)) * q_MPC(5,i) / M(3,3);
    g_0(6,6) = -D(3,3) / M(3,3);
    Coriolis_MPC = [0 0 -M(2,2)*q_MPC(5,i);
                0 0 M(1,1)*q_MPC(4,i);
                M(2,2)*q_MPC(5,i) -M(1,1)*q_MPC(4,i) 0];
    J_dot = [-q_MPC(6,i)*sin(q_MPC(3,i)) -q_MPC(6,i)*cos(q_MPC(3,i)) 0;
             q_MPC(6,i)*cos(q_MPC(3,i)) -q_MPC(6,i)*sin(q_MPC(3,i)) 0;
             0 0 0];
    eta_MPC_b = g_0(1:3,4:6) \ (eta_MPC(:,i) - J_dot * q_MPC(4:6,i));
    tao_MPC(:,i) = M * eta_MPC_b + Coriolis_MPC * q_hat_MPC(4:6,i) + D * q_MPC(4:6,i);    

    f_MPC(1,i) = tao_MPC(1,i)*0.5 + tao_MPC(3,i)/B;
    f_MPC(2,i) = tao_MPC(1,i)*0.5 - tao_MPC(3,i)/B;
    
    %% Simulate disturbance
    tao_wave(1,i) = ExwaveForce(PSD_wave(1,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(1,i));
    tao_wave(2,i) = ExwaveForce(PSD_wave(2,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(2,i));
    tao_wave(3,i) = ExwaveForce(PSD_wave(3,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(3,i));
    tao_wind(:,i) = ExwindForce(beta_w(1,i),V_w(1,i),q_FSF(4:6,i),q_FSF(3,i),rho_air,L);
    %tao_MPC(:,i) = tao_MPC(:,i) - tao_wind(:,i) - tao_wave(:,i);
end
toc  
  
```
  
# Nonlinear MPC
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/NMPC.JPG>
</p><p align="center">      
 
```matlab
%% NMPC variables define
q_NMPC = zeros(nr,N+p);
q_NMPC(1:3,1) = [0.5 1 pi/6]';
z_NMPC = zeros(nr*p,1);
q_dot_NMPC = zeros(nr,N+p);
tao_NMPC = zeros(nc,N+p);
f_NMPC = zeros(2,N);
F_NMPC = zeros(6,6,N+p);
G_NMPC = t*cat(1,zeros(3),inv(M));
Hu_NMPC = kron(eye(p),G_NMPC); % wrt to M_NMPC = zeros(nr*p,nc*p);
Hz_NMPC = zeros(nr*p+nr,nr*p); % wrt to M_NMPC = zeros(nr*p,nc*p);

M_NMPC = zeros(nr*p,nc*p); 
Q_NMPC = 2500 * kron(eye(p),eye(nr)); % 2500
L_NMPC = zeros(nr*p,nr);
R_NMPC = 0.4 * kron(eye(p),eye(nc)); % 0.5
eta_NMPC = zeros(nc,N+p);
new_eta = zeros(nc*p,1);

%% Estimator part initialization
q_hat_NMPC = zeros(6,N+p); 
P_guess = diag([5 -5 1 1 1 1].^2);  % Set the initial pose covariance estimate as a diagonal matrix
P_hat_NMPC = zeros(6,6,N); % Define covariance matrix
P_hat_NMPC(:,:,1) = P_guess;
q_hat_NMPC(:,1) = q_guess + transpose(aa); % set the initial guess for estimator
u_m_NMPC = zeros(3,N+p);
  
 tic
for i = 2:N
    [q_dot_NMPC(1:3,i-1),q_dot_NMPC(4:6,i-1)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_NMPC(4:6,i-1),tao_NMPC(:,i-1),q_NMPC(1:3,i-1),t);
    [q_NMPC(1:3,i),q_NMPC(4:6,i)]= Vehicle(q_dot_NMPC(1:3,i-1),q_dot_NMPC(4:6,i-1),q_NMPC(1:3,i-1),q_NMPC(4:6,i-1),t); 
     
    %% Model propriocdeptive sensors with processing noise
    u_m_NMPC(1,i-1) = tao_NMPC(1,i-1) + sqrt(Q_EKF(1,1))*randn(1); 
    u_m_NMPC(2,i-1) = tao_NMPC(2,i-1) + sqrt(Q_EKF(2,2))*randn(1);
    u_m_NMPC(3,i-1) = tao_NMPC(3,i-1) + sqrt(Q_EKF(3,3))*randn(1);
        
    %% Model GPS,COMPASS and IMU with noise    
    z_m(1,i) = q_NMPC(1,i) + sqrt(R_EKF(1,1))*randn(1);
    z_m(2,i) = q_NMPC(2,i) + sqrt(R_EKF(2,2))*randn(1);
    z_m(3,i) = q_NMPC(3,i) + sqrt(R_EKF(3,3))*randn(1);
    z_m(4,i) = q_NMPC(4,i) + sqrt(R_EKF(4,4))*randn(1);
    z_m(5,i) = q_NMPC(5,i) + sqrt(R_EKF(5,5))*randn(1);
    z_m(6,i) = q_NMPC(6,i) + sqrt(R_EKF(6,6))*randn(1);
    
    %% Simulate EKF fusion of GPS 
    [q_hat_NMPC(:,i),P_hat_NMPC(:,:,i)] = Sensor_EKF(D,M,X_u_dot,Y_v_dot,t,q_hat_NMPC(:,i-1),u_m_NMPC(:,i-1),z_m(:,i),R_EKF,Q_EKF,P_hat_NMPC(:,:,i-1));
    % Update L and M
     for j =1:p 
             [q_dot_NMPC(1:3,i-1+j),q_dot_NMPC(4:6,i-1+j)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_hat_NMPC(4:6,i-1+j),u_m_NMPC(:,i-1+j),q_hat_NMPC(1:3,i-1+j),t);
             [q_hat_NMPC(1:3,i+j),q_hat_NMPC(4:6,i+j)]= Vehicle(q_dot_NMPC(1:3,i-1+j),q_dot_NMPC(4:6,i-1+j),q_hat_NMPC(1:3,i-1+j),q_hat_NMPC(4:6,i-1+j),t);
             
             F_NMPC(:,:,i-1+j) = [1 0 -t*(q_hat_NMPC(4,i-1+j)*sin(q_hat_NMPC(3,i-1+j))+q_hat_NMPC(5,i-1+j)*cos(q_hat_NMPC(3,i-1+j))) t*cos(q_hat_NMPC(3,i-1+j)) -t*sin(q_hat_NMPC(3,i-1+j)) 0;...
             0 1 t*(q_hat_NMPC(4,i-1+j)*cos(q_hat_NMPC(3,i-1+j))-q_hat_NMPC(5,i-1+j)*sin(q_hat_NMPC(3,i-1+j))) t*sin(q_hat_NMPC(3,i-1+j)) t*cos(q_hat_NMPC(3,i-1+j)) 0;...
             0 0 1 0 0 t;...
             0 0 0 1-D(1,1)*(t/M(1,1)) M(2,2)*t*q_hat_NMPC(6,i-1+j)/M(1,1) M(2,2)*t*q_hat_NMPC(5,i-1+j)/M(1,1);...
             0 0 0 -M(1,1)*q_hat_NMPC(6,i-1+j)*(t/M(2,2)) 1-D(2,2)*(t/M(2,2)) -M(1,1)*t*q_hat_NMPC(4,i-1+j)/M(2,2);...
             0 0 0 -t*(X_u_dot-Y_v_dot)*q_hat_NMPC(5,i-1+j)/M(3,3) -t*(X_u_dot-Y_v_dot)*q_hat_NMPC(4,i-1+j)/M(3,3) 1-D(3,3)*(t/M(3,3))];  
             Hz_NMPC(nr*(j+1)-5:nr*(j+1),nr*j-5:nr*j) = F_NMPC(:,:,i-1+j);
             Zd_NMPC(nr*j-nr+1:nr*j,1) = q_d_new(:,i+j);
             z_NMPC(nr*j-nr+1:nr*j,1) = q_hat_NMPC(:,i-1+j);
             tao_lookahead(nc*j-nc+1:nc*j,1) = tao_NMPC(:,i-1+j);
      end
    % Update M matrix
    M_NMPC = (eye(length(Hz_NMPC(1:nr*p,1:nr*p)))-Hz_NMPC(1:nr*p,1:nr*p)) \ Hu_NMPC; % We call it H*
    
    % Update Nonlinear gain
    K_NMPC = transpose(M_NMPC)*Q_NMPC*M_NMPC + R_NMPC;
    tao_delta = K_NMPC \ (transpose(M_NMPC)*Q_NMPC*(Zd_NMPC - z_NMPC) - R_NMPC*tao_lookahead);
    new_eta = tao_lookahead + tao_delta;
    tao_NMPC(:,i) =  new_eta(1:3);
    %[tao_NMPC(1,i),tao_NMPC(3,i)] = sat(tao_NMPC(1,i),tao_NMPC(3,i));
    f_NMPC(1,i) = tao_NMPC(1,i)*0.5 + tao_NMPC(3,i)/B;
    f_NMPC(2,i) = tao_NMPC(1,i)*0.5 - tao_NMPC(3,i)/B;
    
    % Simulate disturbance
    tao_wave(1,i) = ExwaveForce(PSD_wave(1,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(1,i));
    tao_wave(2,i) = ExwaveForce(PSD_wave(2,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(2,i));
    tao_wave(3,i) = ExwaveForce(PSD_wave(3,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(3,i));
    tao_wind(:,i) = ExwindForce(beta_w(1,i),V_w(1,i),q_FSF(4:6,i),q_FSF(3,i),rho_air,L);   
    %tao_NMPC(:,i) = tao_NMPC(:,i) - tao_wind(:,i) - tao_wave(:,i);
end   
```

# Hybrid MPC+FBL
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/HMPC.JPG>
</p><p align="center">   

```matlab
%% HMPC variables define
q_HMPC = zeros(nr,N+p);
q_dot_HMPC = zeros(nr,N+p);
q_HMPC(1:3,1) = [0.5 1 pi/6]';
z_HMPC = zeros(nr*p,N);
eta_HMPC = zeros(nc,N+p);
tao_HMPC = zeros(nc,N);
f_HMPC = zeros(2,N);
delta_z = zeros(nr,N);
eta_lookbehind = zeros(nc*p,1);

F_HMPC = eye(6) + t * A_FSF;
G_HMPC = t * B_FSF; 
Q_HMPC = 0.3 * kron(eye(p),eye(nr));% 0.3 1 //0.3
R_HMPC = 0.2 * kron(eye(p),eye(nc)); % 0.095 0.3135 for the second //0.2

for i =1:p
    L_HMPC(nr*i-nr+1:nr*i, 1:nr) = mpower(F_HMPC, i);
    for j =1:p-i+1
        M_HMPC(nr*(p-i+1)-nr+1:nr*(p-i+1), nc*j-nc+1:nc*j) = mpower(F_HMPC, p-i-j+1) * G_HMPC;
    end
end
K_HMPC = (transpose(M_HMPC)*Q_HMPC*M_HMPC + R_HMPC);

%% Estimator part initialization
q_hat_HMPC = zeros(6,N+p); 
P_guess = diag([5 -5 1 1 1 1].^2);  % Set the initial pose covariance estimate as a diagonal matrix
P_hat_HMPC = zeros(6,6,N); % Define covariance matrix
P_hat_HMPC(:,:,1) = P_guess;
q_hat_HMPC(:,1) = q_guess + transpose(aa); % set the initial guess for estimator 
  
tic
for i = 2:N  
    [q_dot_HMPC(1:3,i-1),q_dot_HMPC(4:6,i-1)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_HMPC(4:6,i-1),tao_HMPC(:,i-1),q_HMPC(1:3,i-1),t);
    [q_HMPC(1:3,i),q_HMPC(4:6,i)]= Vehicle(q_dot_HMPC(1:3,i-1),q_dot_HMPC(4:6,i-1),q_HMPC(1:3,i-1),q_HMPC(4:6,i-1),t); 
     
     %% Model propriocdeptive sensors with processing noise
    u_m(1,i-1) = tao_HMPC(1,i-1) + sqrt(Q_EKF(1,1))*randn(1); 
    u_m(2,i-1) = tao_HMPC(2,i-1) + sqrt(Q_EKF(2,2))*randn(1);
    u_m(3,i-1) = tao_HMPC(3,i-1) + sqrt(Q_EKF(3,3))*randn(1);
        
    %% Model GPS,COMPASS and IMU with noise    
    z_m(1,i) = q_HMPC(1,i) + sqrt(R_EKF(1,1))*randn(1);
    z_m(2,i) = q_HMPC(2,i) + sqrt(R_EKF(2,2))*randn(1);
    z_m(3,i) = q_HMPC(3,i) + sqrt(R_EKF(3,3))*randn(1);
    z_m(4,i) = q_HMPC(4,i) + sqrt(R_EKF(4,4))*randn(1);
    z_m(5,i) = q_HMPC(5,i) + sqrt(R_EKF(5,5))*randn(1);
    z_m(6,i) = q_HMPC(6,i) + sqrt(R_EKF(6,6))*randn(1);
    
    %% Simulate EKF fusion of GPS 
    [q_hat_HMPC(:,i),P_hat_HMPC(:,:,i)] = Sensor_EKF(D,M,X_u_dot,Y_v_dot,t,q_hat_HMPC(:,i-1),u_m(:,i-1),z_m(:,i),R_EKF,Q_EKF,P_hat_HMPC(:,:,i-1));
    
    for j = 1:p
            g_0(1,4) = cos(q_hat_HMPC(3,i-1+j)); 
            g_0(1,5) = -sin(q_hat_HMPC(3,i-1+j)); 
            g_0(2,4) = sin(q_hat_HMPC(3,i-1+j)); 
            g_0(2,5) = cos(q_hat_HMPC(3,i-1+j)); 
            g_0(3,6) = 1;
            g_0(4,4) = -D(1,1)/M(1,1);
            g_0(4,5) = M(2,2)*q_hat_HMPC(6,i-1+j)/M(1,1);
            g_0(5,5) = -D(2,2)/M(2,2);
            g_0(5,6) = M(1,1) * q_hat_HMPC(4,i-1+j)/M(2,2);
            g_0(6,4) = (M(1,1)-M(2,2)) * q_hat_HMPC(5,i-1+j) / M(3,3);
            g_0(6,6) = -D(3,3) / M(3,3);
            Coriolis_HMPC = [0 0 -M(2,2)*q_hat_HMPC(5,i-1+j);
                        0 0 M(1,1)*q_hat_HMPC(4,i-1+j);
                        M(2,2)*q_hat_HMPC(5,i-1+j) -M(1,1)*q_hat_HMPC(4,i-1+j) 0];
            J_dot = [-q_hat_HMPC(6,i-1+j)*sin(q_hat_HMPC(3,i-1+j)) -q_hat_HMPC(6,i-1+j)*cos(q_hat_HMPC(3,i-1+j)) 0;
                    q_hat_HMPC(6,i-1+j)*cos(q_hat_HMPC(3,i-1+j)) -q_hat_HMPC(6,i-1+j)*sin(q_hat_HMPC(3,i-1+j)) 0;
                    0 0 0];           
            tao_head = M*(g_0(1:3,4:6)^(-1))*(eta_HMPC(:,i-1+j) - J_dot * q_hat_HMPC(4:6,i-1+j)) + Coriolis_HMPC * q_hat_HMPC(4:6,i-1+j) + D * q_hat_HMPC(4:6,i-1+j);               
            
            [q_dot_HMPC(1:3,i-1+j),q_dot_HMPC(4:6,i-1+j)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_hat_HMPC(4:6,i-1+j),tao_head,q_hat_HMPC(1:3,i-1+j),t);
            [q_hat_HMPC(1:3,i+j),q_hat_HMPC(4:6,i+j)]= Vehicle(q_dot_HMPC(1:3,i-1+j),q_dot_HMPC(4:6,i-1+j),q_hat_HMPC(1:3,i-1+j),q_hat_HMPC(4:6,i-1+j),t); 
            eta_lookbehind(nc*j-nc+1:nc*j,1) = eta_HMPC(:,i-1+j);            
            z_HMPC(nr*j-nr+1,i) = q_hat_HMPC(1,i+j) - q_d_new(1,i+j);
            z_HMPC(nr*j-nr+2,i) = q_dot_HMPC(1,i+j) - q_d_dot_new(1,i+j);
            z_HMPC(nr*j-nr+3,i) = q_hat_HMPC(2,i+j) - q_d_new(2,i+j);
            z_HMPC(nr*j-nr+4,i) = q_dot_HMPC(2,i+j) - q_d_dot_new(2,i+j);
            z_HMPC(nr*j-nr+5,i) = q_hat_HMPC(3,i+j) - q_d_new(3,i+j);
            z_HMPC(nr*j-nr+6,i) = q_dot_HMPC(3,i+j) - q_d_dot_new(3,i+j);  
    end
    delta_z(:,i) = z_HMPC(1:6,i) - z_HMPC(1:6,i-1);
    new_eta = -K_HMPC \ (transpose(M_HMPC)*Q_HMPC*((z_HMPC(:,i) + L_HMPC * delta_z(:,i))) + R_HMPC * eta_lookbehind); % delta eta
    eta = new_eta + eta_lookbehind;
    eta_HMPC(:,i) = eta(1:3);
    g_0(1,4) = cos(q_HMPC(3,i)); 
    g_0(1,5) = -sin(q_HMPC(3,i)); 
    g_0(2,4) = sin(q_HMPC(3,i)); 
    g_0(2,5) = cos(q_HMPC(3,i)); 
    g_0(3,6) = 1;
    g_0(4,4) = -D(1,1)/M(1,1);
    g_0(4,5) = M(2,2)*q_HMPC(6,i)/M(1,1);
    g_0(5,5) = -D(2,2)/M(2,2);
    g_0(5,6) = M(1,1) * q_HMPC(4,i)/M(2,2);
    g_0(6,4) = (M(1,1)-M(2,2)) * q_HMPC(5,i) / M(3,3);
    g_0(6,6) = -D(3,3) / M(3,3);
    Coriolis_HMPC = [0 0 -M(2,2)*q_HMPC(5,i);
                0 0 M(1,1)*q_HMPC(4,i);
                M(2,2)*q_HMPC(5,i) -M(1,1)*q_HMPC(4,i) 0];
    J_dot = [-q_HMPC(6,i)*sin(q_HMPC(3,i)) -q_HMPC(6,i)*cos(q_HMPC(3,i)) 0;
             q_HMPC(6,i)*cos(q_HMPC(3,i)) -q_HMPC(6,i)*sin(q_HMPC(3,i)) 0;
             0 0 0];
    eta_HMPC_b = g_0(1:3,4:6) \ (eta_HMPC(:,i) - J_dot * q_HMPC(4:6,i));
    tao_HMPC(:,i) = M * eta_HMPC_b + Coriolis_HMPC * q_HMPC(4:6,i) + D * q_HMPC(4:6,i);
    %[tao_HMPC(1,i),tao_HMPC(3,i)] = sat(tao_HMPC(1,i),tao_HMPC(3,i));
    f_HMPC(1,i) = tao_HMPC(1,i)*0.5 + tao_HMPC(3,i)/B;
    f_HMPC(2,i) = tao_HMPC(1,i)*0.5 - tao_HMPC(3,i)/B;
    
    %% Simulate disturbance
    tao_wave(1,i) = ExwaveForce(PSD_wave(1,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(1,i));
    tao_wave(2,i) = ExwaveForce(PSD_wave(2,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(2,i));
    tao_wave(3,i) = ExwaveForce(PSD_wave(3,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_FSF(:,i),Drift_wave(3,i));
    tao_wind(:,i) = ExwindForce(beta_w(1,i),V_w(1,i),q_FSF(4:6,i),q_FSF(3,i),rho_air,L);
    %tao_HMPC(:,i) = tao_HMPC(:,i) - tao_wind(:,i) - tao_wave(:,i);
end
toc    
```  
  
  
  
  
