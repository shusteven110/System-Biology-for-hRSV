clc;clear; 
GRN_name = textread('gene_symbol.txt','%s'); 
[~,~,list]=xlsread('gene symbol expression.xlsx');
PPI_list = {};
GRN_list = {};
load('gene_dir.mat');
Inlnc = cellstr(Inlnc);InMir = cellstr(InMir);InProtein = cellstr(InProtein);InTF = cellstr(InTF);InRcp = cellstr(InRcp);
InVir = {'NS1';'NS2';'N';'P';'M';'SH';'G';'F';'M2';'L'};
Protein = cat(1,InProtein,InRcp,InTF);
Name = cat(1,Protein,InMir,Inlnc,InVir);
Protein = cat(1,Protein,InVir);
%% construct PPI and GRN_list
% for i = 1:size(GRN_name,1)
%     if ismember(GRN_name(i),Protein)
%         PPI_list = cat(1,PPI_list,list(i,:));
%     end
%     if ismember(GRN_name(i),Name)
%         GRN_list = cat(1,GRN_list,list(i,:));
%     end
% end
for i = 1:size(Protein,1)
    index = find(strcmp(GRN_name,Protein{i}));
    PPI_list = cat(1,PPI_list,list(index,:));
end
for i = 1:size(Name,1)
    index = find(strcmp(GRN_name,Name{i}));
    GRN_list = cat(1,GRN_list,list(index,:));
end
save('assignment_temp.mat');
% %% Arrangement
% fprintf('PPI Construction');
% % [PPI]=interMtx_20151228(Protein,Protein,1);
% [PPI]=candidate_network_construction(Protein,Protein,"PPI");
% fprintf('GRN Construction');
% % [GRN]=interMtx2_20160324(Name,Name,2);
% [GRN]=candidate_network_construction(Name,Name,"GRN");
% save('assignment_temp.mat');
%% add virus candidate
PPI_size = size(PPI_list,1);
GRN_size = size(GRN_list,1);
VIR_size = size(InVir,1);
for i = 1:VIR_size
    PPI(PPI_size-i+1,1:end)=1;
    PPI(1:end,PPI_size-i+1)=1;
    PPI(PPI_size-i+1,PPI_size-i+1)=0;
    GRN(GRN_size-i+1,1:end)=1;
    GRN(GRN_size-i+1,GRN_size-i+1)=0;
end
PPI_name = Protein;
GRN_name = Name;
save('GRN_name.mat','GRN_name');save('PPI_name.mat','PPI_name');
save('PPI.mat','PPI');save('GRN.mat','GRN');
%% construct GRN_EXP 
fprintf('GRN_EXP Construction\n');
V=cell2mat(GRN_list(1:end,3:10));
%linspace(start,end,linspace)
yy_EBV=spline([0 60 120 240 360 480 600 720],V,linspace(0,720,720)); 

yyV_NAN=find(isnan(yy_EBV));%check
yyV_INF=find(isinf(yy_EBV));%check

V_max=max(V,[],2);
yy_EBV(find(yy_EBV<0))=0;

for j=1:size(V,1);
     yy_EBV(j,find(yy_EBV(j,:)>max(V(j,:))))=max(V(j,:));
 end
 yyV_max=max(yy_EBV,[],2);
 GRN_EXP=sparse(yy_EBV);
 save('GRN_EXP.mat','GRN_EXP')

%% construct PPI_EXP
fprintf('PPI_EXP Construction\n');
V=cell2mat(PPI_list(1:end,3:10));
%linspace(start,end,linspace)
yy_EBV=spline([0 60 120 240 360 480 600 720],V,linspace(0,720,720)); 

yyV_NAN=find(isnan(yy_EBV));%check
yyV_INF=find(isinf(yy_EBV));%check

V_max=max(V,[],2);
yy_EBV(find(yy_EBV<0))=0;

for j=1:size(V,1);
     yy_EBV(j,find(yy_EBV(j,:)>max(V(j,:))))=max(V(j,:));
 end
 yyV_max=max(yy_EBV,[],2);
 PPI_EXP=sparse(yy_EBV);
 save('PPI_EXP.mat','PPI_EXP')

% for i=1:6
%     for j=1:6
%         if i==1 a=1;b=14389; fprintf('protein ');end
%         if i==2 a=14390;b=17471; fprintf('rcp ');end
%         if i==3 a=17472;b=19065; fprintf('tf ');end
%         if i==4 a=19066;b=19616; fprintf('mi ');end
%         if i==5 a=19617;b=22861; fprintf('lnc ');end
%         if i==6 a=22862;b=22871; fprintf('vir ');end
%         if j==1 c=1;d=14389; fprintf('protein ');end
%         if j==2 c=14390;d=17471; fprintf('rcp ');end
%         if j==3 c=17472;d=19065; fprintf('tf ');end
%         if j==4 c=19066;d=19616; fprintf('mi ');end
%         if j==5 c=19617;d=22861; fprintf('lnc ');end
%         if j==6 c=22862;d=22871; fprintf('vir ');end
%         try
%             num = sum(sum(GRN(a:b,c:d)));
%             fprintf("%d\n",full(num));
%         catch
%             fprintf('0\n');
%         end
%     end
% end
 