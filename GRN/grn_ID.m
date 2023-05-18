clear;clc;
%% loading data
load('GRN_EXP.mat');load('GRN.mat');load('GRN_name.mat');
interaction = GRN; gene_expre = GRN_EXP; Name = GRN_name;

% ######## you need to specify the columns of microRNA ########
micro = 19066:19616;
% #############################################################

%% parameters for ID
time_point = size(gene_expre,2);  % n times
gene_num = size(gene_expre,1);  %number of gene
Final_reg_ability = double(zeros(size(interaction),'like',interaction));
basal = zeros(size(interaction,1),1);
[BB,index] = sort(sum(interaction,2));
id = find(BB~=0,1,'first');   %you can use for the start point

fprintf('start GRN ID:\n');
acc = 120;
accelerate = 1:acc;  %used when forward model selection
All = 1:size(interaction,2);  % all binds (used for diff to get the genes that don't bind with ith gene)
Ones = ones(time_point-1,1);  %used for phi in linear regression
cons = [gene_num+1,gene_num+2];

%----parameters for lsqlin (algorithm : 'interior-point'(default) )----
options = optimset('Display','off');
lb_tmp = -1*Inf(gene_num+2,1);    %lower bound
ub_tmp = Inf(gene_num+2,1);    %upper bound
ub_tmp(micro) = zeros(length(micro),1);  ub_tmp(end-1) = 1;

% default : 'MaxIter',200,'TolCon',1e-8,'TolFun',1e-8,'TolX',1e-12


%% Start (or continue. To continue, you need to load the saved variables and enter the lastest j.)
File1 = fopen('grn_early_record.txt','a');
try
    load('GRN_ID_TEMP.mat');
    j = j+1;
catch
    j = id;
end
start = j; % start order  (or enter the last j to continue)
total_time = tic;
X_all = gene_expre(:,1:time_point-1)';
phi_temp_all = gene_expre(:,1:time_point-1)';
phi_temp_all(:,micro) = phi_temp_all(:,micro).*pi;

delta_complex1 = -2/(time_point-1);
delta_complex2 = 5/(time_point-1);
for j = start:gene_num
    tic
    i = index(j); % index of the instant protein
    remaining = gene_num-j;  % the amount of protein left without processing (only be used to print information)
    bind1 = find(interaction(i,:));  % ppi = 1 with the instant (ith) protein
    % initial AIC_value
    AIC_value_initial = 90000;
	min_AIC = 100;
    %----linear regression (with all known genes interacting with ith gene)----
    X = gene_expre(i,2:time_point)';
    if ~any(X) || isempty(bind1)
        continue
    end
    pi = X_all(:,i);
	phi_temp = cat(2,phi_temp_all,pi,Ones);  % [PPI Pi(t) Pi(t) 1], Pi(t) are use to multiply the translation and degration rate
    % above variable would skip in the following linear regression
    
    
    % Foward
    Bind1 = {}; Theta = {}; AIC_value_pack = {};
    binding = bind1;
    if length(bind1) < acc  %to accelerate, here skip some combinations of parameter
        tmp = 1:length(bind1);
    else
        tmp = accelerate;
    end
    for m = tmp  
        AIC_value_temp = AIC_value_initial;
        for u = 1:length(binding)
            if m == 1
                bindtemp = binding(u);
            else
                bindtemp = cat(2,Bind1{m-1},binding(u));
            end
            
            %----linear regression forward----
            Sort = [sort(bindtemp),cons];
            phi = phi_temp(:,Sort);
            lb = lb_tmp(Sort);
            ub = ub_tmp(Sort);
            phi_len = size(phi,2);
            
            [theta_temp,resnorm_temp] = lsqlin(phi,X,[],[],[],[],lb,ub,[],options); % return theata that minimize the function and get the value(resnorm)
            %----linear regression end and calculate AIC ----
            pre_AIC_value = log( resnorm_temp/(time_point-1) ) + 2*phi_len/(time_point-1);  %AIC
            
            if u == 1 && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                AIC_value_temp = pre_AIC_value;
                bindinfo = bindtemp;
            elseif pre_AIC_value<=AIC_value_temp && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                AIC_value_temp = pre_AIC_value;
                bindinfo = bindtemp;
            end
        end
        disp(m)
        if m ~= 1 && (AIC_value_pack{m-1}-AIC_value_temp) < delta_complex1 && AIC_value_temp - min_AIC > delta_complex2
            break
        end
        if AIC_value_temp < min_AIC
            min_AIC = AIC_value_temp;
        end
        Bind1{m} = bindinfo; Theta{m} = theta_temp; AIC_value_pack{m} = AIC_value_temp;  %record variable
        binding = setdiff(bind1,Bind1{m});
    end
    [~,min_index] = min(cell2mat(AIC_value_pack));    
    bindout = Bind1{min_index};
    index1 = min_index;
    Bind1 = {bindout}; Theta = {Theta{min_index}}; AIC_value_pack = {AIC_value_pack{min_index}};
    % Backward
    if length(bindout) > 1
        for y = 1:length(bindout)-2  
            AIC_value_temp = AIC_value_initial;
            disp(y)
            for z = 1:length(bindout)
                bindtemp = bindout;
                bindtemp(z) = [];
                
                %----linear regression backward----
                Sort = [sort(bindtemp),cons];
                phi = phi_temp(:,Sort);
                lb = lb_tmp(Sort);
                ub = ub_tmp(Sort);
                phi_len = size(phi,2);
                
                [theta_temp,resnorm_temp] = lsqlin(phi,X,[],[],[],[],lb,ub,[],options);  % return theata that minimize the function and get the value(resnorm)
                %----linear regression end  and calculate AIC----
                pre_AIC_value = log( resnorm_temp/(time_point-1) ) + 2*phi_len/(time_point-1);  %AIC
                if y ~= 1 && (AIC_value_pack{y-1}-AIC_value_temp) < delta_complex1 && AIC_value_temp - min_AIC > delta_complex2
                    break
                end
                if AIC_value_temp < min_AIC
                    min_AIC = AIC_value_temp;
                end
                if z == 1 && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                    AIC_value_temp = pre_AIC_value;
                    bindinfo = bindtemp;
                elseif pre_AIC_value<=AIC_value_temp && ~isnan(pre_AIC_value) && ~isinf(pre_AIC_value)
                    AIC_value_temp = pre_AIC_value;
                    bindinfo = bindtemp;
                end
            end
            Bind1{1+y} = bindinfo; Theta{1+y} = theta_temp; AIC_value_pack{1+y} = AIC_value_temp;
            bindout = bindinfo;
        end
    end
    
    [~,min_index] = min(cell2mat(AIC_value_pack));
    thetaout = Theta{min_index};
    bindout = Bind1{min_index};
	
    
    basal(i) = thetaout(end); % basal level for epigenetic regulation
    Final_reg_ability(i,bindout) = thetaout(1:end-2);
    
    save ('GRN_ID_TEMP.mat','j','Final_reg_ability','basal');
    fprintf('%d (%d nodes --> %d nodes)  ',j,length(bind1),length(bindout))
    toc % start time
end
fclose(File1);
fprintf('Done\n');
A = Final_reg_ability;
fprintf('Interaction:[%6d ------> %-6d]\n',size(find(interaction~=0),1),size(find(A~=0),1))
fprintf('       Node:[%6d ------> %-6d]\n',size(interaction,1),length(find(sum(A,2)~=0)))

%
GRN_PNP = A;
GRN_Edge = A;
Basal_GRN = table(Name,basal);
save('GRN_PNP.mat','GRN_PNP', '-v7.3')
save('GRN_Edge.mat','GRN_Edge', '-v7.3')
writetable(Basal_GRN,'Basal_GRN.txt')