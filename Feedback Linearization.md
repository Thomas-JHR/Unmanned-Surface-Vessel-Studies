# Feedback Linearization Technique

Let ![image](https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/z.svg) be the FBL state of the controller and ![image](https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/eta.svg) be the linear FBL control input


<p align="center">
 <img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/12.svg
</p>

The next step is to construct the second order linear system as  
<p align="center">
 <img src=https://github.com/Thomas-JHR/Unmanned-Surface-Vessel-Studies/blob/main/Tex/equation_FBL.svg 
</p> 
  
Based on the equation, it can be easily depicted that the linear control input is defined as an acceleration vector state of the vessel in the global frame  
```
A_FSF = [0 1 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 0 0];
B_FSF = [0 0 0;
        1 0 0;
        0 0 0;
        0 1 0;
        0 0 0;
        0 0 1];
pole_FSF = [-1 -0.35 -1.2 -0.55 -1.3 -0.6];
K_FSF_MATRIX = place(A_FSF,B_FSF,pole_FSF);
```
