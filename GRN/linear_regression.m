 function [X,phi,theta,resnorm]=linear_regression(cons,X,phi_all,bind,lb_tmp,ub_tmp,options)
phi = phi_all(:,[sort(bind),cons]);   
lb = lb_tmp([sort(bind),cons]);
ub = ub_tmp([sort(bind),cons]);
try
    phi = full(phi);
    [theta,resnorm]=lsqlin(phi,X,[],[],[],[],lb,ub,[],options); % resnorm = (X-phi*theta)'*(X-phi*theta)
catch
    theta = phi\X;
    tmp = A*theta;
    resnorm = (X-phi*theta)'*(X-phi*theta);
    if tmp(1) > 0 || tmp(2) > 1
        phi = full(phi);
        [theta,resnorm]=lsqlin(phi,X,[],[],[],[],lb,ub,[],options); % resnorm = (X-phi*theta)'*(X-phi*theta)
    end
end

