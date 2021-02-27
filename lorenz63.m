%  Created by Sagar K. Tamang (January 2021)
%
%  This program performs integration of chaotic Lorenz-63 system using
%  fourth-order Runge Kutta method 
%  
%  The algorithm used in this code is referenced from the following:
%  Tamang et. al. "Ensemble Riemannian Data Assimilation over the 
%  Wasserstein Space"

function xout = lorenz63(xin,dt,parm)
    sigma = parm(1);
    rho = parm(2);
    beta = parm(3);
    
    w1 = 1/6; w2 = 1/3; w3 = 1/3; w4 = 1/6;
    
    k1 = dt*lorenz63rhs(xin,sigma,rho,beta);
    k2 = dt *lorenz63rhs(xin+0.5*k1,sigma,rho,beta);
    k3 = dt *lorenz63rhs(xin+0.5*k2,sigma,rho,beta);
    k4 = dt *lorenz63rhs(xin+k3,sigma,rho,beta);
    
    xout = xin + w1*k1 + w2*k2 + w3*k3 + w4*k4;
    
    function rhs = lorenz63rhs(x,sigma,rho,beta)
        rhs  = [-sigma*(x(1)-x(2)),rho*x(1)-x(2)-x(1)*x(3),x(1)*x(2)-beta*x(3)];
    end

end
