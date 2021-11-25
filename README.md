This repository contains the source code for implementing Ensemble Riemannian Data Assimilation. The code is referenced from the following: Tamang et. al. Ensemble Riemannian Data Assimilation over the Wasserstein Space.

There are four included matlab code files:

(1) Demo.m : This code provides an example implementation of the Ensemble Riemannian Data Assimilation on the chaotic Lorenz-63 system. Details of the experimental setting can be found in Tamang et. al. Ensemble Riemannian Data Assimilation over the Wasserstein Space.

(2) entrop_OMT.m: This function approximates the optimal transportation plan using entropic regularization.

(3) lorenz63.m: This function integrates the Lorenz model using fourth-order Runge-Kutta Scheme.
 
(4) Covariance.m: Sample covariance computation


![EnRDA](https://user-images.githubusercontent.com/46984734/143469929-f0ce5514-5254-4a26-adea-33a18ad13e1e.PNG)
![Fig2](https://user-images.githubusercontent.com/46984734/143469962-962a42c1-0033-4409-ac83-55e4d0a7b8b3.png)
![Fig7](https://user-images.githubusercontent.com/46984734/143469969-9da52e68-c5b9-4ffb-89c5-ae43a2a7f80a.png)
