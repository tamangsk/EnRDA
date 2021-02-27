This repository contains the source code for implementing Ensemble Riemannian Data Assimilation. The code is referenced from the following: Tamang et. al. Ensemble Riemannian Data Assimilation over the Wasserstein Space.

There are three included matlab code files:

(1) Demo.m : This code provides an example implementation of the Ensemble Riemannian Data Assimilation on the chaotic Lorenz-63 system. Details of the experimental setting can be found in Tamang et. al. Ensemble Riemannian Data Assimilation over the Wasserstein Space.

(2) entrop_OMT.m: This function approximates the optimal transportation plan using entropic regularization.

(3) lorenz63.m: This function integrates the Lorenz model using fourth-order Runge-Kutta Scheme.
 
(4) Covariance.m: Sample covariance computation