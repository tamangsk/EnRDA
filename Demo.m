%  Created by Sagar K. Tamang (February 2021)
%
%  This program provides an example implementation of the Ensemble
%  Riemannian Data Assimilation on chaotic Lorenz-63 system
%  
%  The algorithm used in this code is referenced from the following:
%  Tamang et. al. "Ensemble Riemannian Data Assimilation over the 
%  Wasserstein Space"


clear; clc; close all;
tic
options = odeset('RelTol',1e-6,'AbsTol',1e-12*ones(3,1));

parm = [10,28,8/3]; 
parm_sys = [10.5,27,10/3]; 
dt = 0.01; 
t_f = 20;    
tspan = 0:dt:t_f; 
xtrue = [1.508870,-1.531271,25.46091];      

gamma = 10;                                 
sig_mod = sqrt(2*dt);                       
sig_init = sqrt(2);                         
sig_obs = sqrt(2);                          
N = 10;                                     
obs_int = 40*dt;                            
corr_mat = [1 0.5 0.25;                     
            0.5 1 0.5;
            0.25 0.5 1];
obs_num = 41:40:length(tspan);
t_obs = tspan(obs_num);
Obs_cov = sig_obs^2.*corr_mat;              
mod_cov = sig_mod^2.*eye(3);                

xa_EnRDA = nan(length(tspan),3,N);

for j = 1:N
    xa_EnRDA(1,:,j) = xtrue + normrnd(0,sig_init);    
end

count = 1;

for i = 2 : length(tspan)
    xtrue(i,:) = lorenz63(xtrue(i-1,:),dt,parm);
end
Obs = xtrue(41:40:end,:) + mvnrnd(zeros(3,1),Obs_cov,50); 

for i = 2 : length(tspan) 
    for j = 1:N
        xa_EnRDA(i,:,j) = lorenz63(xa_EnRDA(i-1,:,j),dt,parm_sys) + mvnrnd(zeros(3,1),mod_cov);
    end
    
    if i == obs_num(count)
        x = reshape(xa_EnRDA(i,:,:),[3,N])';
        px = ones(N,1)./N;
        
        for j = 1:N
            y(j,:) = Obs(count,:) + mvnrnd(zeros(3,1),Obs_cov);
        end
        
        B = covariance(x',x');
        
        eta = trace(Obs_cov)/trace(Obs_cov+B);
        
        py = ones(N,1)./N;
        
        U = entrop_OMT(x,y,px,py,gamma,300);
        
        [I,J,U1] = find(U);
        
        Xt = eta*x(I,:)+(1-eta)*y(J,:);
        
        for j = 1:N
            xa_EnRDA(i,:,j) = Xt(find(rand<cumsum(U1),1,'first'),:);
        end
        count = count + 1;
    end
end


subplot(1,3,1)
plot(tspan,xtrue(:,1),'color','k')   
hold on 
plot(tspan,mean(xa_EnRDA(:,1,:),3),'color','r')
scatter(t_obs,Obs(:,1),15,'+','g')

subplot(1,3,2)
plot(tspan,xtrue(:,2),'color','k')   
hold on 
plot(tspan,mean(xa_EnRDA(:,2,:),3),'color','r')
scatter(t_obs,Obs(:,2),15,'+','g')

subplot(1,3,3)
plot(tspan,xtrue(:,3),'color','k')   
hold on 
plot(tspan,mean(xa_EnRDA(:,3,:),3),'color','r')
scatter(t_obs,Obs(:,3),15,'+','g')
legend({'Truth','$\mathbf{x}_a$(EnRDA)','Obs.'},'interpreter','latex')
legend boxoff

xa_EnRDA = mean(xa_EnRDA,3);
mse_EnRDA = mse(xtrue(:),xa_EnRDA(:));
bias_EnRDA = abs(mean(xtrue(:)-xa_EnRDA(:)));
ubrmse_EnRDA = sqrt(mse_EnRDA-bias_EnRDA^2);