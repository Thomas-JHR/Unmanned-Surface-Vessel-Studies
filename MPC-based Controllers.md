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
  
  
  
  
  
  
  
