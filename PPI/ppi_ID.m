clear;clc;
%% loading data
load('PPI_EXP.mat');load('PPI.mat');load('PPI_name.mat');
interaction=PPI; gene_expre=PPI_EXP; Name = PPI_name;
nnn=size(PPI_EXP,2); % n times
Final_reg_ability = double(zeros(size(interaction),'like',interaction));
basal = zeros(size(interaction,1),1);
[BB,index]=sort(sum(interaction(:,1:end)~=0,2));
id=find(BB~=0,1,'first'); 
options=optimset('Algorithm','interior-point','Display','off');
%% ID
try
    load('PPI_ID.mat');
    j = j + 1;
catch
    j =1;
end
fprintf('start ID(from %d to %d)\n',size(gene_expre,1)-id,0)
X_all=gene_expre(:,2:nnn)';
phi_all=gene_expre(:,1:nnn-1)'; 
start = j; % start order
clear PPI PPI_EXP PPI_name

delta_complex1 = -2/(nnn-1);
delta_complex2 = 5/(nnn-1);
for j = start:18999
    tic % start time
    i = index(j); % index of the instant protein
    mmm = size(gene_expre,1)-j; % the amount of protein left without processing
    bind1=find(interaction(i,:)==1); % ppi information for 1 of the instant protein
    bindout = bind1;
    % initial AIC_value
    AIC_value_initial = 10000;
    min_AIC = 100;
    phi_cat = cat(2,phi_all(:,i),phi_all(:,i),ones(nnn-1,1));
    pi = phi_all(:,i);
    X = X_all(:,i);
    [X,phi,theta,resnorm]=linear_regression(i,nnn,X,phi_all,bind1,phi_cat,pi,options);
    AIC_value=AICValue(X,theta,resnorm);
    thetaout=theta;
    if length(bind1) > 1
        binding=bind1;
        Bind1 = cell(1,length(binding)+ nnn/3);
        Theta = cell(1,length(binding)+ nnn/3);
        AIC_value_pack = cell(1,length(binding)+ nnn/3);
        for m = 1:240% Foward
            for u = 1:length(binding)
                if m==1
                    bindtemp = binding(u);
                else
                    bindtemp = cat(2,Bind1{m-1},binding(u));
                end
                % linear regression
                [X,phi,theta_temp,resnorm_temp]=linear_regression(i,nnn,X,phi_all,bindtemp,phi_cat,pi,options);
                pre_AIC_value=AICValue(X,theta_temp,resnorm_temp);
                if u==1
                    if ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                        AIC_value_temp = pre_AIC_value;
                        bindinfo = bindtemp;
                    else
                        AIC_value_temp = AIC_value_initial;
                        bindinfo = bindtemp;
                    end
                else
                    if pre_AIC_value<=AIC_value_temp && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                        AIC_value_temp = pre_AIC_value;
                        bindinfo = bindtemp;
                    end
                end
            end
            if m ~= 1 && (AIC_value_pack{m-1}-AIC_value_temp) < delta_complex1 && AIC_value_temp - min_AIC > delta_complex2
                break
            end
            if AIC_value_temp < min_AIC
                min_AIC = AIC_value_temp;
            end
            Bind1{m}= bindinfo;
            Theta{m}= theta_temp;
            AIC_value_pack{m}=AIC_value_temp;
            binding = setdiff(bind1,Bind1{m});
        end
        [~,min_index] = min(cell2mat(AIC_value_pack));
        bindout= Bind1{min_index};
        if length(bindout) > 1
            for y=1:length(bindout)-1 % Backward
                for z = 1:length(bindout)
                    bindtemp = bindout;
                    bindtemp(z)= [];
                    % linear regression
                    [X,phi,theta_temp,resnorm_temp]=linear_regression(i,nnn,X,phi_all,bindtemp,phi_cat,pi,options);
                    pre_AIC_value=AICValue(X,theta_temp,resnorm_temp);
                    if z == 1
                        if ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                            AIC_value_temp = pre_AIC_value;
                            bindinfo = bindtemp;
                        else 
                            AIC_value_temp = AIC_value_initial;
                            bindinfo = bindtemp;
                        end
                    else
                        if pre_AIC_value<=AIC_value_temp && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                            AIC_value_temp = pre_AIC_value;
                            bindinfo = bindtemp;
                        end
                    end
                end
                bindout= bindinfo;
                Bind1{m+y}= bindinfo;
                Theta{m+y}= theta_temp;
                AIC_value_pack{m+y}=AIC_value_temp;
            end
        end
        [AIC_value,min_index] = min(cell2mat(AIC_value_pack));
        thetaout=Theta{min_index};
        bindout= Bind1{min_index};
    end
    basal(i,:) = thetaout(end); % basal level for epigenetic regulation
    null = find(thetaout(1:end-3)==0);
    bindout(null)= [];thetaout(null)= [];
    for n = 1:length(bindout)
        Final_reg_ability(i,bindout(n))=thetaout(n);
    end    
    if mod(mmm,1)==0
        save ('PPI_ID_TEMP_18750.mat','j','Final_reg_ability','basal');
        fprintf('Saved\n');
    end  
    fprintf('%d (%d nodes --> %d nodes)  ',mmm,length(bind1),length(bindout))
    toc % start time
end
load('PPI_ID.mat');
% load('PPI.mat');
load('PPI_name.mat');

Name = PPI_name;
interaction = PPI;
A = PPI;
fprintf('Done\n')
fprintf('Interaction:[%6d ------> %-6d]\n',size(find(interaction~=0),1),size(find(A~=0),1))
fprintf('       Node:[%6d ------> %-6d]\n',size(interaction,1),length(find(sum(A,2)~=0)))
% toc % elapsed time
%%
PPI_PNP = A;
PPI_Edge = (A+A')./2;
Basal_PPI = table(Name,basal);
save('PPI_PNP.mat','PPI_PNP', '-v7.3')
save('PPI_Edge.mat','PPI_Edge', '-v7.3')
writetable(Basal_PPI,'Basal_PPI.txt')