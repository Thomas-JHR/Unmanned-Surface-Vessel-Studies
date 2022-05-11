# Gaussian Process Regression

A Gaussian Process Regression approach, regarded as non-parametric Bayesian regression, is proposed to combine either the PD controller or the hybrid MPC to enhance tracking performance under environmental forces. The basic principle of the GPR approach is generate a prediction function by observed dataset training, and predict outputs as close as possible based on the new dataset and the function learning in statistical learning theory [Machine]. To understand some basics behind the learning tool, this section involves some backgrounds for the GPR algorithm to make prediction of the future inputs and uncertainty learning implementation. The proposed GPR solution was already applied in the hydraulic manipulator design for improving desired angle tracking performance in [Jack], and the Heron M300 vessel model takes this method to test if the vessel can overcome the disturbances.

Let ```X``` be the observed training dataset, ```X*``` be the testing dataset, and ```y``` defined as a label or an output vector. We introduce a kernel function or covariance function of the GPR, also known as squared exponential function between the training dataset and testing dataset. It can be described as

 <p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/GPR1.svg>
</p><p align="center">

```σ_f^2``` and ```l``` are defined as tunable signal-variance and length-scale hyperparameters vectors. These hypermeters are optimized during the model training and outputs prediction by maximizing the log of the marginal likelihood defined as,
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/GPR2.svg>
</p><p align="center">

```σ_n^2``` is the noise variance and n is the total number of observed training dataset points. The GPR method generates mean and covariance of the predicted output vectors. Their mathematical expression is presented as   
  
 <p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/GPR3.svg>
</p><p align="center">

```σ_y^2``` is the observed output noise.   
  
### GPR+PD+FBL
  
<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/GPR1.JPG width=70% height=70%>
</p><p align="center">
  
```matlab
  
  %% linear GPR+PD+FBL define
q_GPFSF = zeros(6,N); % x y psi u v r
tao_GPFSF = zeros(3,N);
f_GPFSF = zeros(2,N);
q_GPFSF(1:3,1) = [0.5 1 pi/6]';% Set up initial conditions
q_dot_GPFSF = zeros(6,N);
err_GPFSF = zeros(6,N);
new_a_n = zeros(3,N);
new_a_b = zeros(3,N);
new_a_b_new = zeros(3,N);

%% Estimator part initialization
q_hat_GPFSF = zeros(6,N);
aa = [0.3 -0.3 0.1 0.1 0.1 0.1]; 
P_hat_GPFSF = zeros(6,6,N); % Define covariance matrix
P_hat_GPFSF(:,:,1) = P_guess;
q_hat_GPFSF(:,1) = q_guess + transpose(aa); % set the initial guess for estimator
u_m_GPFSF = zeros(3,N);
z_m_GPFSF = zeros(6,N);

%% GPR variables define
a_GPFSF = zeros(9,N); % [q;qdot;a_b]; 
g_function_GPFSF = zeros(3,N);
sigma_noise = 0.002; % noise variance
dX_PD = size(a_GPFSF,1);
dY_PD = size(g_function_GPFSF,1);

meanfun1_GPFSF = [];
covfun1_GPFSF = @covSEiso;
likfun1_GPFSF = @likGauss;
meanfun2_GPFSF = [];
covfun2_GPFSF = @covSEiso;
likfun2_GPFSF = @likGauss;
meanfun3_GPFSF = [];
covfun3_GPFSF = @covSEiso;
likfun3_GPFSF = @likGauss;

hyp1_GPFSF = struct('mean',[],'cov',[1;1], 'lik', -1);
hyp2_GPFSF = struct('mean',[],'cov',[1;1], 'lik', -1);
hyp3_GPFSF = struct('mean',[],'cov',[1;1], 'lik', -1);
hyp1_min_GPFSF = hyp1_GPFSF;
hyp2_min_GPFSF = hyp2_GPFSF;
hyp3_min_GPFSF = hyp3_GPFSF;
  
MaxOber = 20;

%% Update data to the system
for i = 2:N
    [q_dot_GPFSF(1:3,i-1),q_dot_GPFSF(4:6,i-1)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_GPFSF(4:6,i-1),tao_GPFSF(:,i-1),q_GPFSF(1:3,i-1),t);
    [q_GPFSF(1:3,i),q_GPFSF(4:6,i)]= Vehicle(q_dot_GPFSF(1:3,i-1),q_dot_GPFSF(4:6,i-1),q_GPFSF(1:3,i-1),q_GPFSF(4:6,i-1),t); 
    
     %% Model propriocdeptive sensors with processing noise
    u_m_GPFSF(1,i-1) = tao_GPFSF(1,i-1) + sqrt(Q_EKF(1,1))*randn(1); 
    u_m_GPFSF(2,i-1) = tao_GPFSF(2,i-1) + sqrt(Q_EKF(2,2))*randn(1);
    u_m_GPFSF(3,i-1) = tao_GPFSF(3,i-1) + sqrt(Q_EKF(3,3))*randn(1);
        
    %% Model GPS,COMPASS and IMU with noise    
    z_m_GPFSF(1,i) = q_GPFSF(1,i) + sqrt(R_EKF(1,1))*randn(1);
    z_m_GPFSF(2,i) = q_GPFSF(2,i) + sqrt(R_EKF(2,2))*randn(1);
    z_m_GPFSF(3,i) = q_GPFSF(3,i) + sqrt(R_EKF(3,3))*randn(1);
    z_m_GPFSF(4,i) = q_GPFSF(4,i) + sqrt(R_EKF(4,4))*randn(1);
    z_m_GPFSF(5,i) = q_GPFSF(5,i) + sqrt(R_EKF(5,5))*randn(1);
    z_m_GPFSF(6,i) = q_GPFSF(6,i) + sqrt(R_EKF(6,6))*randn(1);
    
    %% Simulate EKF fusion of GPS 
    [q_hat_GPFSF(:,i),P_hat_GPFSF(:,:,i)] = Sensor_EKF(D,M,X_u_dot,Y_v_dot,t,q_hat_GPFSF(:,i-1),u_m_GPFSF(:,i-1),z_m_GPFSF(:,i),R_EKF,Q_EKF,P_hat_GPFSF(:,:,i-1));
    
    err_GPFSF(1,i) = q_hat_GPFSF(1,i-1) - q_d(1,i);
    err_GPFSF(2,i) = q_dot_GPFSF(1,i-1) - q_d_dot(1,i);
    err_GPFSF(3,i) = q_hat_GPFSF(2,i-1) - q_d(2,i);
    err_GPFSF(4,i) = q_dot_GPFSF(2,i-1) - q_d_dot(2,i);
    err_GPFSF(5,i) = q_hat_GPFSF(3,i-1) - q_d(3,i);
    err_GPFSF(6,i) = q_dot_GPFSF(3,i-1) - q_d_dot(3,i);
    new_a_n(:,i) = -K_FSF_MATRIX * err_GPFSF(:,i);
    g_0(1,4) = cos(q_hat_GPFSF(3,i)); 
    g_0(1,5) = -sin(q_hat_GPFSF(3,i)); 
    g_0(2,4) = sin(q_hat_GPFSF(3,i)); 
    g_0(2,5) = cos(q_hat_GPFSF(3,i)); 
    g_0(3,6) = 1;
    g_0(4,4) = -D(1,1)/M(1,1);
    g_0(4,5) = M(2,2)*q_hat_GPFSF(6,i)/M(1,1);
    g_0(5,5) = -D(2,2)/M(2,2);
    g_0(5,6) = M(1,1) * q_hat_GPFSF(4,i)/M(2,2);
    g_0(6,4) = (M(1,1)-M(2,2)) * q_hat_GPFSF(5,i) / M(3,3);
    g_0(6,6) = -D(3,3) / M(3,3);
    Coriolis_GPFSF = [0 0 -M(2,2)*q_hat_GPFSF(5,i);
                      0 0 M(1,1)*q_hat_GPFSF(4,i);
                      M(2,2)*q_hat_GPFSF(5,i) -M(1,1)*q_hat_GPFSF(4,i) 0];
    J_dot = [-q_hat_GPFSF(6,i)*sin(q_hat_GPFSF(3,i)) -q_hat_GPFSF(6,i)*cos(q_hat_GPFSF(3,i)) 0;
            q_hat_GPFSF(6,i)*cos(q_hat_GPFSF(3,i)) -q_hat_GPFSF(6,i)*sin(q_hat_GPFSF(3,i)) 0;
            0 0 0];     
   new_a_b(:,i) = g_0(1:3,4:6) \ (new_a_n(:,i) - J_dot * q_hat_GPFSF(4:6,i)); % PD control acceleration
    %% GP training   
Using GP_Training function to train the sample data and predict the output response, followed by updating the data to the system 
    g_function_GPFSF(:,i-1) = q_dot_GPFSF(4:6,i-1) - new_a_b_new(:,i-1) + sigma_noise*randn(3,1); % Ydataset
    a_GPFSF(:,i) = [err_GPFSF(:,i);new_a_b(:,i)]; % Xdataset(Train + Test) 
    if i > MaxOber
        n = i - MaxOber;
        [~, ~, g_function_GPFSF(1,i), ~] = gp(hyp1_min_GPFSF,@infGaussLik,meanfun1_GPFSF,covfun1_GPFSF,likfun1_GPFSF,a_GPFSF(:,n:i-1)',g_function_GPFSF(1,n:i-1)',a_GPFSF(:,i)');
        [~, ~, g_function_GPFSF(2,i), ~] = gp(hyp2_min_GPFSF,@infGaussLik,meanfun2_GPFSF,covfun2_GPFSF,likfun2_GPFSF,a_GPFSF(:,n:i-1)',g_function_GPFSF(2,n:i-1)',a_GPFSF(:,i)');
        [~, ~, g_function_GPFSF(3,i), ~] = gp(hyp3_min_GPFSF,@infGaussLik,meanfun3_GPFSF,covfun3_GPFSF,likfun3_GPFSF,a_GPFSF(:,n:i-1)',g_function_GPFSF(3,n:i-1)',a_GPFSF(:,i)');
    else
         n = min(i,MaxOber);
        [~, ~, g_function_GPFSF(1,n), ~] = gp(hyp1_min_GPFSF,@infGaussLik,meanfun1_GPFSF,covfun1_GPFSF,likfun1_GPFSF,a_GPFSF(:,1:n-1)',g_function_GPFSF(1,1:n-1)',a_GPFSF(:,n)');
        [~, ~, g_function_GPFSF(2,n), ~] = gp(hyp2_min_GPFSF,@infGaussLik,meanfun2_GPFSF,covfun2_GPFSF,likfun2_GPFSF,a_GPFSF(:,1:n-1)',g_function_GPFSF(2,1:n-1)',a_GPFSF(:,n)');
        [~, ~, g_function_GPFSF(3,n), ~] = gp(hyp3_min_GPFSF,@infGaussLik,meanfun3_GPFSF,covfun3_GPFSF,likfun3_GPFSF,a_GPFSF(:,1:n-1)',g_function_GPFSF(3,1:n-1)',a_GPFSF(:,n)');
    end
    
    if abs(g_function_GPFSF(1,i)) > 1 
        g_function_GPFSF(1,i) = sign(g_function_GPFSF(1,i)) * 1;
    end      
    if abs(g_function_GPFSF(2,i)) > 1
        g_function_GPFSF(2,i) = sign(g_function_GPFSF(2,i)) * 1;
    end
    if abs(g_function_GPFSF(3,i)) > 1
        g_function_GPFSF(3,i) = sign(g_function_GPFSF(3,i)) * 1;
    end
    
    if mod(i,MaxOber) == 0 && i > MaxOber
        hyp1_min_GPFSF = minimize(hyp1_GPFSF,@gp,-5,@infGaussLik,meanfun1_GPFSF,covfun1_GPFSF,likfun1_GPFSF,a_GPFSF(:,n:i-1)',g_function_GPFSF(1,n:i-1)');
        hyp2_min_GPFSF = minimize(hyp2_GPFSF,@gp,-5,@infGaussLik,meanfun2_GPFSF,covfun2_GPFSF,likfun2_GPFSF,a_GPFSF(:,n:i-1)',g_function_GPFSF(2,n:i-1)');
        hyp3_min_GPFSF = minimize(hyp3_GPFSF,@gp,-5,@infGaussLik,meanfun3_GPFSF,covfun3_GPFSF,likfun3_GPFSF,a_GPFSF(:,n:i-1)',g_function_GPFSF(3,n:i-1)');
    elseif mod(i,MaxOber) == 0 && i <= MaxOber
        hyp1_min_GPFSF = minimize(hyp1_GPFSF,@gp,-5,@infGaussLik,meanfun1_GPFSF,covfun1_GPFSF,likfun1_GPFSF,a_GPFSF(:,1:i-1)',g_function_GPFSF(1,1:i-1)');
        hyp2_min_GPFSF = minimize(hyp2_GPFSF,@gp,-5,@infGaussLik,meanfun2_GPFSF,covfun2_GPFSF,likfun2_GPFSF,a_GPFSF(:,1:i-1)',g_function_GPFSF(2,1:i-1)');
        hyp3_min_GPFSF = minimize(hyp3_GPFSF,@gp,-5,@infGaussLik,meanfun3_GPFSF,covfun3_GPFSF,likfun3_GPFSF,a_GPFSF(:,1:i-1)',g_function_GPFSF(3,1:i-1)');
    end
    
    %% GP training

    new_a_b_new(:,i) =  new_a_b(:,i) - g_function_GPFSF(:,i); 
    tao_GPFSF(:,i) = M * new_a_b_new(:,i) + Coriolis_GPFSF * q_GPFSF(4:6,i) + D * q_GPFSF(4:6,i);
    f_GPFSF(1,i) = tao_GPFSF(1,i)*0.5 + tao_GPFSF(3,i)/B;
    f_GPFSF(2,i) = tao_GPFSF(1,i)*0.5 - tao_GPFSF(3,i)/B;
    
    %% Saturation the force
    %[tao_GPFSF(1,i),tao_GPFSF(3,i)] = sat(tao_GPFSF(1,i),tao_GPFSF(3,i));
    
    %% Simulate disturbances    
    tao_wave(1,i) = ExwaveForce(PSD_wave(1,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_GPFSF(:,i), Drift_wave(1,i));
    tao_wave(2,i) = ExwaveForce(PSD_wave(2,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_GPFSF(:,i), Drift_wave(2,i));
    tao_wave(3,i) = ExwaveForce(PSD_wave(3,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_GPFSF(:,i), Drift_wave(3,i));
    tao_wind(:,i) = ExwindForce(beta_w(1,i),V_w(1,i),q_GPFSF(4:6,i),q_GPFSF(3,i),rho_air,L);
    %tao_GPFSF(:,i) = tao_GPFSF(:,i) -tao_wind(:,i) - tao_wave(:,i);    
end
toc  
```

### GPR+HMPC+FBL

<p align="center">
<img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/GPR2.JPG width=70% height=70%>
</p><p align="center">
  
```matlab
%% linear hybrid MPC+GPR define
q_GPMPC = zeros(6,N); % x y psi u v r
tao_GPMPC = zeros(3,N);
f_GPMPC = zeros(2,N);
q_GPMPC(1:3,1) = [0.5 1 pi/6]';% Set up initial conditions
q_dot_GPMPC = zeros(6,N+p);
z_GPMPC = zeros(nr*p,N);
eta_GPMPC = zeros(nc,N+p);
eta_GPMPC_b = zeros(nc,N+p);
tao_head_GPMPC = zeros(nc,N+p);

%% Estimator part initialization
q_hat_GPMPC = zeros(6,N);
P_guess = diag([5 -5 1 1 1 1].^2);  % Set the initial pose covariance estimate as a diagonal matrix
P_hat_GPMPC = zeros(6,6,N); % Define covariance matrix
P_hat_GPMPC(:,:,1) = P_guess;
q_hat_GPMPC(:,1) = q_guess + transpose(aa); % set the initial guess for estimator
u_m_GPMPC = zeros(3,N);
z_m_GPMPC = zeros(6,N);
delta_z_GPMPC = zeros(nr,N);
GP+MPC parameters initialization
%% GPR variables define
a_GPMPC = zeros(9,N+p); % [q;qdot;a_b]; 
g_function_GPMPC = zeros(6,N+p);
dX_GPMPC = size(a_GPMPC,1);
dY_GPMPC = size(g_function_GPMPC,1);

meanfun_GPMPC = [];    
covfun_GPMPC = @covSEiso;    
likfun_GPMPC = @likGauss;   
hyp_GPMPC = struct('mean',[],'cov',[1;1], 'lik', -1);
hyp_min_GPMPC = hyp_GPMPC;  

MaxGP = 200;
tic
for i = 2:N
    [q_dot_GPMPC(1:3,i-1),q_dot_GPMPC(4:6,i-1)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_GPMPC(4:6,i-1),tao_GPMPC(:,i-1),q_GPMPC(1:3,i-1),t);
    [q_GPMPC(1:3,i),q_GPMPC(4:6,i)]= Vehicle(q_dot_GPMPC(1:3,i-1),q_dot_GPMPC(4:6,i-1),q_GPMPC(1:3,i-1),q_GPMPC(4:6,i-1),t);
    
    %% Model propriocdeptive sensors with processing noise
    u_m_GPMPC(1,i-1) = tao_GPMPC(1,i-1) + sqrt(Q_EKF(1,1))*randn(1); 
    u_m_GPMPC(2,i-1) = tao_GPMPC(2,i-1) + sqrt(Q_EKF(2,2))*randn(1);
    u_m_GPMPC(3,i-1) = tao_GPMPC(3,i-1) + sqrt(Q_EKF(3,3))*randn(1);
        
    %% Model GPS,COMPASS and IMU with noise    
    z_m_GPMPC(1,i) = q_GPMPC(1,i) + sqrt(R_EKF(1,1))*randn(1);
    z_m_GPMPC(2,i) = q_GPMPC(2,i) + sqrt(R_EKF(2,2))*randn(1);
    z_m_GPMPC(3,i) = q_GPMPC(3,i) + sqrt(R_EKF(3,3))*randn(1);
    z_m_GPMPC(4,i) = q_GPMPC(4,i) + sqrt(R_EKF(4,4))*randn(1);
    z_m_GPMPC(5,i) = q_GPMPC(5,i) + sqrt(R_EKF(5,5))*randn(1);
    z_m_GPMPC(6,i) = q_GPMPC(6,i) + sqrt(R_EKF(6,6))*randn(1);
    
    %% Simulate EKF fusion of GPS 
    [q_hat_GPMPC(:,i),P_hat_GPMPC(:,:,i)] = Sensor_EKF(D,M,X_u_dot,Y_v_dot,t,q_hat_GPMPC(:,i-1),u_m_GPMPC(:,i-1),z_m_GPMPC(:,i),R_EKF,Q_EKF,P_hat_GPMPC(:,:,i-1));
   
    %% Prediction horizon loop
    for j = 1:p
         g_0(1,4) = cos(q_hat_GPMPC(3,i-1+j)); 
         g_0(1,5) = -sin(q_hat_GPMPC(3,i-1+j)); 
         g_0(2,4) = sin(q_hat_GPMPC(3,i-1+j)); 
         g_0(2,5) = cos(q_hat_GPMPC(3,i-1+j)); 
         g_0(3,6) = 1;
         g_0(4,4) = -D(1,1)/M(1,1);
         g_0(4,5) = M(2,2)*q_hat_GPMPC(6,i-1+j)/M(1,1);
         g_0(5,5) = -D(2,2)/M(2,2);
         g_0(5,6) = M(1,1) * q_hat_GPMPC(4,i-1+j)/M(2,2);
         g_0(6,4) = (M(1,1)-M(2,2)) * q_hat_GPMPC(5,i-1+j) / M(3,3);
         g_0(6,6) = -D(3,3) / M(3,3);
         Coriolis_GPMPC = [0 0 -M(2,2)*q_hat_GPMPC(5,i-1+j);
                        0 0 M(1,1)*q_hat_GPMPC(4,i-1+j);
                        M(2,2)*q_hat_GPMPC(5,i-1+j) -M(1,1)*q_hat_GPMPC(4,i-1+j) 0];
         J_dot = [-q_hat_GPMPC(6,i-1+j)*sin(q_hat_GPMPC(3,i-1+j)) -q_hat_GPMPC(6,i-1+j)*cos(q_hat_GPMPC(3,i-1+j)) 0;
                    q_hat_GPMPC(6,i-1+j)*cos(q_hat_GPMPC(3,i-1+j)) -q_hat_GPMPC(6,i-1+j)*sin(q_hat_GPMPC(3,i-1+j)) 0;
                    0 0 0];
         eta_GPMPC_b(:,i-1+j) = (g_0(1:3,4:6)^(-1))*(eta_GPMPC(:,i-1+j) - J_dot * q_hat_GPMPC(4:6,i-1+j));
         tao_head_GPMPC(:,i-1+j) = M * eta_GPMPC_b(:,i-1+j)   + Coriolis_GPMPC * q_hat_GPMPC(4:6,i-1+j) + D * q_hat_GPMPC(4:6,i-1+j);               
            
         [q_dot_GPMPC(1:3,i-1+j),q_dot_GPMPC(4:6,i-1+j)] = HERON(X_u_dot,Y_v_dot,N_r_dot,I_z,q_hat_GPMPC(4:6,i-1+j),tao_head_GPMPC(:,i-1+j),q_hat_GPMPC(1:3,i-1+j),t);
         [q_hat_GPMPC(1:3,i+j),q_hat_GPMPC(4:6,i+j)]= Vehicle(q_dot_GPMPC(1:3,i-1+j),q_dot_GPMPC(4:6,i-1+j),q_hat_GPMPC(1:3,i-1+j),q_hat_GPMPC(4:6,i-1+j),t); 
         eta_lookbehind(nc*j-nc+1:nc*j,1) = eta_GPMPC(:,i-1+j);            
         z_GPMPC(nr*j-nr+1,i) = q_hat_GPMPC(1,i+j) - q_d_new(1,i+j);
         z_GPMPC(nr*j-nr+2,i) = q_dot_GPMPC(1,i+j) - q_d_dot_new(1,i+j);
         z_GPMPC(nr*j-nr+3,i) = q_hat_GPMPC(2,i+j) - q_d_new(2,i+j);
         z_GPMPC(nr*j-nr+4,i) = q_dot_GPMPC(2,i+j) - q_d_dot_new(2,i+j);
         z_GPMPC(nr*j-nr+5,i) = q_hat_GPMPC(3,i+j) - q_d_new(3,i+j);
         z_GPMPC(nr*j-nr+6,i) = q_dot_GPMPC(3,i+j) - q_d_dot_new(3,i+j);
         %% GP Training Dataset Collection
         %a_GPMPC(:,i-1+j) = [z_GPMPC(1:6,i);eta_GPMPC_b(:,i-1+j)]; % Xtrain + Xtest   
         a_GPMPC(:,i-1+j) =[reshape(z_GPMPC(nr*j-nr+1:nr*j-nr+6,i),[6,1]);eta_GPMPC_b(:,i-1+j)];
    end    
    g_function_GPMPC(:,i-1) = z_GPMPC(1:6,i) + sigma_noise*randn(6,1); % Ytrain
    
    if i > MaxGP   
        %% GP model training and prediction for p horizons
        [~, ~, Testing, ~] = gp(hyp_min_GPMPC,@infGaussLik,meanfun_GPMPC,covfun_GPMPC,likfun_GPMPC,a_GPMPC(:,i-MaxGP:i-1)',g_function_GPMPC(:,i-MaxGP:i-1)',a_GPMPC(:,i:i+p-1)');  
        g_function_GPMPC(:,i:i+p-1) = transpose(Testing);
        for k = 1:6
            for x = i:i+p-1
                 if abs(g_function_GPMPC(k,x)) > 1 
                     g_function_GPMPC(k,x) = sign(g_function_GPMPC(k,x)) * 1;
                 end 
            end
        end
        hyp_min_GPMPC = minimize(hyp_GPMPC,@gp,-5,@infGaussLik,meanfun_GPMPC,covfun_GPMPC,likfun_GPMPC,a_GPMPC(:,i-MaxGP:i-1)',g_function_GPMPC(:,i-MaxGP:i-1)');
        z_GPMPC(:,i) = z_GPMPC(:,i) + reshape(g_function_GPMPC(:,i:i+p-1),[600,1]);
    end
    %% MPC Updates
    delta_z_GPMPC(:,i) = z_GPMPC(1:6,i) - z_GPMPC(1:6,i-1);
    new_eta = -K_HMPC \ (transpose(M_HMPC)*Q_HMPC*((z_GPMPC(:,i) + L_HMPC * delta_z_GPMPC(:,i))) + R_HMPC * eta_lookbehind); % delta eta
    eta = new_eta + eta_lookbehind;
    eta_GPMPC(:,i) = eta(1:3);
    g_0(1,4) = cos(q_GPMPC(3,i)); 
    g_0(1,5) = -sin(q_GPMPC(3,i)); 
    g_0(2,4) = sin(q_GPMPC(3,i)); 
    g_0(2,5) = cos(q_GPMPC(3,i)); 
    g_0(3,6) = 1;
    g_0(4,4) = -D(1,1)/M(1,1);
    g_0(4,5) = M(2,2)*q_GPMPC(6,i)/M(1,1);
    g_0(5,5) = -D(2,2)/M(2,2);
    g_0(5,6) = M(1,1) * q_GPMPC(4,i)/M(2,2);
    g_0(6,4) = (M(1,1)-M(2,2)) * q_GPMPC(5,i) / M(3,3);
    g_0(6,6) = -D(3,3) / M(3,3);
    Coriolis_MPC = [0 0 -M(2,2)*q_GPMPC(5,i);
                0 0 M(1,1)*q_GPMPC(4,i);
                M(2,2)*q_GPMPC(5,i) -M(1,1)*q_GPMPC(4,i) 0];
    J_dot = [-q_GPMPC(6,i)*sin(q_GPMPC(3,i)) -q_GPMPC(6,i)*cos(q_GPMPC(3,i)) 0;
             q_GPMPC(6,i)*cos(q_GPMPC(3,i)) -q_GPMPC(6,i)*sin(q_GPMPC(3,i)) 0;
             0 0 0];
    eta_GPMPC_b(:,i) = g_0(1:3,4:6) \ (eta_GPMPC(:,i) - J_dot * q_GPMPC(4:6,i));
    tao_GPMPC(:,i) = M * eta_GPMPC_b(:,i) + Coriolis_MPC * q_GPMPC(4:6,i) + D * q_GPMPC(4:6,i);    
    f_GPMPC(1,i) = tao_GPMPC(1,i)*0.5 + tao_GPMPC(3,i)/B;
    f_GPMPC(2,i) = tao_GPMPC(1,i)*0.5 - tao_GPMPC(3,i)/B;
    
    %% Simulate disturbance
    tao_wave(1,i) = ExwaveForce(PSD_wave(1,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_GPMPC(:,i), Drift_wave(1,i));
    tao_wave(2,i) = ExwaveForce(PSD_wave(2,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_GPMPC(:,i), Drift_wave(2,i));
    tao_wave(3,i) = ExwaveForce(PSD_wave(3,i),beta_w(1,i), omega_o, lambda_wave,omega(i), q_GPMPC(:,i), Drift_wave(3,i));
    tao_wind(:,i) = ExwindForce(beta_w(1,i),V_w(1,i),q_GPMPC(4:6,i),q_GPMPC(3,i),rho_air,L);
    %tao_GPMPC(:,i) = tao_GPMPC(:,i) -tao_wind(:,i) - tao_wave(:,i);    
end
toc  
```
  
