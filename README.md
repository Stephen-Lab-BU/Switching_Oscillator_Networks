# State_Space_Oscillator_Models

Code repository for all of the methods described and the simulation studies performed in the paper Switching Models of Oscillatory Networks Greatly Improve Inference of Dynamic Functional Connectivity.

The detailed explantion of the code is at the top of each Matlab script file. Below is a brief description of all scripts in this repo.

* core_functions
  * skf.m -- switching Kalman filter
  * smoother.m -- switching RTS smoother
  * em_B.m -- EM on B matrices (for the Common Oscillator Model)
  * em_projA.m -- EM on A matrices (for the Correlated Noise Model)
  * em_projQ.m -- EM on Q (Sigma) matrices (for the Directed Influence Model)

* helper_functions
  * funcs2.m -- a combination of helper functions for simulation settings and visualizations
    
* toy_examples
  * sim_COM_toy.m -- simulation file for the toy example of the Common Oscillator Model (COM)
  * sim_CNM_toy.m -- simulation file for the toy example of the Correlated Noise Model (CNM)
  * sim_DIM_toy.m -- simulation file for the toy example of the Directed Influence Model (DIM)

* model_evaluations
  * sim_COM_eval.m -- simulation file for model evaluations of the Common Oscillator Model (COM)
  * sim_CNM_eval.m -- simulation file for model evaluations of the Correlated Noise Model (CNM)
  * sim_DIM_eval.m -- simulation file for model evaluations of the Directed Influence Model (DIM)
  * eval_COM_fit.m -- model evaluations for the fits from the Common Oscillator Model
  * eval_CNM_fit.m -- model evaluations for the fits from the Correlated Noise Model
  * eval_DIM_fit.m -- model evaluations for the fits from the Directed Influence Model
  * testing_files -- including a set of model fits (results obtained from the EM) from the COM, CNM, and DIM. Files in this folder can be used to test if the eval_COM_fit.m, eval_CNM_fit.m and eval_DIM_fit.m are working. For example, "COM_sim_e10_mleQ.mat" corresponds to results that were simulated under the Common Oscillator Model (COM) and fit the data with the Correlated Noise Model (CNM).

