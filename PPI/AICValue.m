function [AIC_value]=AICValue(X,theta,norm)
sigma=norm/length(X); 
AIC_value=log(sigma)+2*length(theta)/length(X);