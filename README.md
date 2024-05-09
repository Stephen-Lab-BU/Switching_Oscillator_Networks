# Switching_Oscillator_Networks

Code repository for all of the methods described and the simulation studies performed in the paper *Switching Models of Oscillatory Networks Greatly Improve Inference of Dynamic Functional Connectivity*, [arXiv:2404.18854](https://arxiv.org/abs/2404.18854).

The detailed explanation of the code is at the top of each Matlab script file. Below is a brief description of all scripts in this repo.

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
  * sim_COM_eval.m -- simulation setting file for model evaluations of the Common Oscillator Model (COM)
  * sim_CNM_eval.m -- simulation setting file for model evaluations of the Correlated Noise Model (CNM)
  * sim_DIM_eval.m -- simulation setting file for model evaluations of the Directed Influence Model (DIM)
  * eval_COM_fit.m -- model evaluations for the fits from the Common Oscillator Model
  * eval_CNM_fit.m -- model evaluations for the fits from the Correlated Noise Model
  * eval_DIM_fit.m -- model evaluations for the fits from the Directed Influence Model
  * testing_files -- including a set of model fits (results obtained from the EM) from the COM, CNM, and DIM. Files in this folder can be used to test if the eval_COM_fit.m, eval_CNM_fit.m and eval_DIM_fit.m are working. For example, "COM_sim_e10_mleQ.mat" corresponds to results that were simulated under the Common Oscillator Model (COM) and fit the data with the Correlated Noise Model (CNM).


Reproduce figures in the paper:
* Section 3.1 Model Visualization:
  * Figure 4: run sim_COM_toy.m from the top to the end
  * Figure 5: run sim_CNM_toy.m from the top to the end
  * Figure 6: run sim_DIM_toy.m from the top to the end

* Section 3.2 Model Performance Evaluations
  * Figure 7: For the sub-figures in the first column, start by running sim_DIM_eval.m to generate simulated data from the Directed Influence Model. The last three chunks perform EM on A, B, and Q, respectively. Then, evaluate the EM results by running eval_DIM_fit.m for outputs from EM on A to obtain the norm, sensitivity, and switching accuracy. Similarly, run eval_COM_fit.m for outputs from EM on B, and eval_CNM_fit.m for outputs from EM on Q. The error bars in the figures are based on four realizations for each SNR level. The first realization uses a random seed rng(35) at the beginning on line 15 and rng(22) for the network structure on line 58. The second realization uses rng(22) for both. The third realization uses rng(50) at the beginning and rng(22) for the network structure. The last realization uses rng(52) at the beginning and rng(20) for the network structure. For the sub-figures in the middle column, start by running sim_CNM_eval.m to generate the simulated data from the Correlated Noise Model. Then, run eval_CNM_fit.m, eval_DIM_fit.m, and eval_COM_fit.m to evaluate EM results from the Correlated Noise Model (Q), Directed Influence Model (A), and Common Oscillator Model (B), respectively. The setup for the four realizations is the same as described above. Finally, for the sub-figures in the third column, start by running sim_COM_eval.m to generate simulated data from the Common Oscillator Model. Then, run eval_COM_fit.m, eval_CNM_fit.m, and eval_DIM_fit.m to evaluate EM results from the Common Oscillator Model (B), Correlated Noise Model (Q), and Directed Influence Model (A), respectively. The setup for the four realizations remains consistent as previously described.
