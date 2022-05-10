# Simple Random Walk

A random walk mechanism acting as the anemometers is firstly applied to generate the speed and angle data from wind for software simulation. ```p_value``` is assumed to be the probability that decides if the next state increases or decreases by ∆=0.009. The general algorithm is described as follows

</p><p align="center">
 <img src= https://user-images.githubusercontent.com/45107735/167704691-4ae1f260-399b-4dcd-b2a0-3a021a23dab0.png>
</p>

```matlab

function [new_P] = RandomWalk(P,Dect)
  p = rand;
  if p < 0.5
    S = -Dect;
  elseif p > 0.5
     S = Dect;
  end
  
  new_P = S + P; % Gives the next random walk from the new position P every time

end
```
# Wind

A wind feed-forward controller designed from is modeled to show how the wind resistance data affects the vessel trajectories. A wind disturbance forces vector  acting on the vessel in the body-fixed frame are considered to be the final output of the controller.```[V_w beta_w]```created by the random walk is assumed to be a true wind velocity and angle vector. The procedure to derive the wind forces is shown below.

</p><p align="center">
 <img src= https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/Wind1.svg>
</p>

```υ_w``` is a wind velocity vector in the body-fixed frame. ```υ_rw``` and ```γ```_rw are a relative wind velocity vector and angle of attack in the body-fixed coordinate in each direction, respectively. ```q_air``` is the dynamic pressure where ρ_air=1.27kg/m^3 is the density of the air. ```A_FW``` and ```A_LW``` are the frontal and lateral projected windage area. ```c_x```, ```c_y```, and ```c_z``` are wind coefficients in surge, sway, and yaw axes. The three wind coefficients are set within the recommended range ```0.5≤c_x≤0.9```,```0.7≤c_y≤0.95```,and ```0.05≤c_z≤0.2``` at the beginning, followed by being manually modified based on the environment and vehicle performance. The coefficients are set to be the mean of the range in the controller model for simplicity. 

```matlab

function tao_wind = ExwindForce(beta_w,V_w,V,psi,rho_air,Laa)

v_w(1) = V_w * cos(beta_w(1) - psi(1));
v_w(2) = V_w * sin(beta_w(1) - psi(1));

v_rw = zeros(2,1); % relative wind velocity vector
v_rw(1) = V(1) - v_w(1);
v_rw(2) = V(2) - v_w(2);

gamma_rw = zeros(1,1);
gamma_rw(1) = -atan(v_rw(2)/v_rw(1));

Afw = 0.202; % front windage area(m^2)
Alw = 0.3025; % lateral windage area (m^2)

cx = 0.7; % surge wind coefficient
cy = 0.825; % sway wind coefficient
cz = 0.125; % yaw wind coefficient

tao_wind(1) = -0.5* rho_air * (v_rw(1)^(2) + v_rw(2)^(2)) * cx * cos(gamma_rw(1)) * Afw;
tao_wind(2) = 0.5* rho_air * (v_rw(1)^(2) + v_rw(2)^(2)) * cy * sin(gamma_rw(1)) * Alw;
tao_wind(3) = 0.5* rho_air * (v_rw(1)^(2) + v_rw(2)^(2)) * cz * sin(2*gamma_rw(1)) * Alw * Laa;
end
```

# Wave

The wave disturbance is assumed to be created by the wind in this project so that the wind angle and direction can be directly used for generating wave spectrum energy. The wave model is designed based on the **Marine System Simulator** library originally created by Thor Fossen.

Wave-induced forces and moment is regarded as wave-frequency motion with zero mean oscillatory motionsThe significance height wave ```H_s```, the mean wave height of the one-third highest wave, can directly tell which sea states the vessel system is operating. This project assumes that a moderate wave disturbance is observed to exert to the Heron M300 vessel. A few wave parameters are defined as the wave model initialization. These include angular peak frequency ```ω_peak```, angular maximum frequency ```ω_max```, and wave sequence period ```ω```. Then, we use one of the wave spectrum approaches, Torsethaugen spectrum from curve-fitting experimental data from the North sea, to create empirical two peaked spectrum for swell (low frequency peak) and newly developed waves (high frequency peak) by ```torset_spec``` function in the library.

Afterwards, the second-order linear wave transfer function is approximated to estimate the wave forces and moment exerted on the Heron vessel in the body-fixed frame at the real time described as 


</p><p align="center">
 <img src= https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/wave1.svg>
</p>

```w_i (i=1,2,…,6)``` are Gaussian white noise processes, ```λ_i (i=1,2,…,6)``` is the damping factors and ```σ_wave^i (i=1,2,…,6)``` describe wave intensity as square root of the maximum power spectrum density. ```ω_e^i (i=1,2,…,6)``` are encounter frequencies, and presented as 

</p><p align="center">
 <img src= https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/Wind2.svg>
</p>

```g``` is the gravity constant. Although there are six degrees of freedom from expressions, only three degrees of freedom are used for modelling the wave force for simplicity. Finally, we need to compute their magnitudes from the wave frequency response by substituting ```s→jω``` and taking ```|τ_wave (jω)|``` as the wave disturbance forces and moment exerted to the vessel. 


```matlab

function tao_wave = ExwaveForce(Wave,beta, omega_peak, lambda_wave,omega_range, q, Drift)

  g = 9.81;
  sigma_wave = sqrt(max(Wave));
  omega_encounter = abs(omega_range - power(omega_range,2)*sqrt(power(q(4),2)+power(beta-q(5),2))/g);
  tao_wave = sqrt(abs( 4*(lambda_wave*omega_peak*sigma_wave)^2*omega_range / ( ((omega_encounter)^2-(omega_range)^2).^2 +4*(lambda_wave*omega_encounter*omega_range)^2)));
  %% Add Gaussian white noise process parameters to the PSD_wave
  tao_wave = awgn(tao_wave,10) + Drift;

end
```
