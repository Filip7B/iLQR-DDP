# Trajectory Optimization for 2 DoF arm
This project uses iLQR/DDP to generate a swing up motion of the arm. The main optimization method used in this code is based on [[1](#1)]. This optimization method can be classified as indirect, shooting method. Each iteration of the optimization consist of a backward and a a forward pass through the time steps. Backward pass runs from the final  k=N to the initial k=0 time step and minimizes cost-to-go function along the way. Forward pass runs in the other direction by forward stepping the dynamics. Each iteration generates a slightly better trajectory in terms of the cost function until convergence. 

The cost function consists of the running and the final cost. The final cost penalizes deviations from the final state, and the running cost penalizes large torques.

The main optimization method is based on [[2](#2)]. It adds, however, the possibility to use more accurate symbolic derivatives.

To generate and visualize the motion run iLQR_and_DDP_for_DoF_Arm.m. Optimization parameters, initial state as well as the time horizon can all be changed in this file as well. If you need help implementing this method in your own use case, feel free to contact me: bojadzievski.f@gmail.com




https://user-images.githubusercontent.com/103943635/164108765-287b14ca-baa8-4d10-aabb-bd6521dcb733.mp4




## References
<a id="1">[1]</a> 
Tassa, Yuval and Mansard, Nicolas and Todorov, Emo,
Control-limited differential dynamic programming,
2014 IEEE International Conference on Robotics and Automation (ICRA)

<a id="2">[2]</a>
Yuval (2022). iLQG/DDP trajectory optimization (https://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization), MATLAB Central File Exchange. Retrieved April 19, 2022.
