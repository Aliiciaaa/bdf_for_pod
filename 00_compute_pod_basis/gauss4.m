function [w,c]=gauss4
% Gaussian quadrature formula for 4 nodes in [0,1].
% nodes
alpha0=sqrt(30)/72;
beta1=0.5*sqrt((15+2*sqrt(30))/35);
beta2=0.5*sqrt((15-2*sqrt(30))/35);
c=[0.5-beta1, 0.5-beta2, 0.5+beta2, 0.5+beta1];
% weights
w=[0.25-alpha0, 0.25+alpha0, 0.25+alpha0, 0.25-alpha0]';
end
