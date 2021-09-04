function [sres, xest, pest, resp, wp]=kmodel1pv(par,stimrep)
% input par : current parameter of the model, here only q/r
% input stimrep : stimulus stimrep(:,1) and response stimrep(:,2)
% output sres : residual of simulation and response
% simulation=resp-sres

% par(1) : ratio q/r for kalman filter

% two-stage model published in
% Petzschner FH, Glasauer S. Iterative bayesian estimation as an explanation 
% for range and regression effects: a study on human path integration. 
% J Neurosci. 2011 Nov 23;31(47):17220-9.
% 
% S.Glasauer 2011/2021

r=1;
q=par(1)*r;
p=r;

a=10;
off=1;

dist=stimrep(:,1);
pest=dist*0;
xest=dist*0;
resp=dist;
wp=dist*0;
for i=1:length(dist) % walk each distance in randomized order
    d=dist(i);
    
    % this is the core estimation routine
    % measurement
    z=log(a*d+off); % transform distance to log
    if i==1 % initial mean prior equals measurement d
        x=z;
    end
    % kalman filter for the mean of the prior
    km=(p+q)/(p+q+r); % kalman gain
    p=km*r; % new variance
    pest(i)=p;
    %x=(1-km)*x+km*z; % prior
    x=x+km*(z-x);
    xest(i)=(exp(x)-off)/a;

    % final fusion
    w1=1/r/(1/p+1/r); % weight of measurement
    xe=(1-w1)*x+w1*z; % estimate
    wp(i)=1-w1; % weight of prior
    % backtransform
    resp(i)=(exp(xe)-off)/a; % estimate backtransform    
end

sres=stimrep(:,2)-resp;
% analytical steady state solutions
% k_final=0.5*q/r*(sqrt(1+4*r/q)-1);
% final prior weight: w_prior=(1-k_final)/(1+k_final);

% remove NaNs ... this treatment of NaNs should lead to better estimation
% than just setting the residuals to zero, since here the missing values 
% are really treated as missing
sres=sres(isfinite(sres));

end
