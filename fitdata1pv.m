function [px,ci,resnorm,simres]=fitdata1pv(stimrep)
% input stimrep : stimulus and response
% output px : parameters of the model kmodel(p,stimrep)
% output ci : confidence interval of px
%
% fit wrapper for the two-stage model published in
% Petzschner FH, Glasauer S. Iterative bayesian estimation as an explanation for range and regression effects: a study on human path integration. J Neurosci. 2011 Nov 23;31(47):17220-9.
% 
% S.Glasauer 2011/2021
p0=1;lb=0;
[px,resnorm,residual,~,~,~,jacobian]=lsqnonlin(@(p)kmodel1pv(p,stimrep),p0,lb);
ci=nlparci(px,residual,'jacobian',jacobian);
[~,~,~,simres]=kmodel1pv(px,stimrep);
% simres=stimrep(:,2)-kmodel(px,stimrep);
end
