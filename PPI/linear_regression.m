 function [X,phi,theta,resnorm]=linear_regression(i,nnn,X,phi_all,bind,phi_cat,pi,options)
% X=X_all(:,i);
% pi=phi_all(:,i);
phi = phi_all(:,sort(bind));   
% phi=cat(2,phi.*pi,pi,pi,ones(nnn-1,1)); % [Gi(t) Pi(t) 1]
phi=cat(2,phi.*pi,phi_cat); % [Gi(t) Pi(t) 1]
b=[0;1]; 
A= [zeros(1,size(phi,2)-3),-1,0,0;zeros(1,size(phi,2)-3),0,1,0];
% theta = phi\X;
% tmp = A*theta;
% resnorm = (X-phi*theta)'*(X-phi*theta);
% if tmp(1) > 0 || tmp(2) > 1
phi = full(phi);
[theta,resnorm]=lsqlin(phi,X,A,b,[],[],[],[],[],options); % resnorm = (X-phi*theta)'*(X-phi*theta)
% end
