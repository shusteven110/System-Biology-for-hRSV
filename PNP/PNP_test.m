clear all, clc;
% load data: PPI and GRN are regulation ability from PPI and GRN
load('PPI_name');  load('GRN_name');
load('PPI_PNP');  load('GRN_PNP');

ppi = PPI_PNP;
grn = GRN_PNP;

gene = 1:19065;
miRNA = 19066:19616;
lncRNA = 19617:22861;
gene_vir = 22862:22871;

PNP_up_name = [PPI_name;GRN_name(miRNA);GRN_name(lncRNA)];
PNP_down_name = [PPI_name;GRN_name(gene);GRN_name(miRNA);GRN_name(lncRNA);GRN_name(gene_vir)];
save('PNP_name.mat','PNP_up_name', 'PNP_down_name');

PPI = sparse(size(ppi,1),size(ppi,2)+length(miRNA)+length(lncRNA));
PPI(1:size(ppi,1),1:size(ppi,2)) = ppi;
HG = sparse(cat(2, grn(gene,gene), zeros(length(gene),length(gene_vir)), grn(gene,miRNA), grn(gene,lncRNA)));
HM = sparse(cat(2, grn(miRNA,gene), zeros(length(miRNA),length(gene_vir)), grn(miRNA,miRNA), grn(miRNA,lncRNA)));
HL = sparse(cat(2, grn(lncRNA,gene), zeros(length(lncRNA),length(gene_vir)), grn(lncRNA,miRNA), grn(lncRNA,lncRNA)));
PG = sparse(cat(2, grn(gene_vir,gene), grn(gene_vir,gene_vir), grn(gene_vir,miRNA), grn(gene_vir,lncRNA)));

M = cat(1,PPI, HG, HM, HL, PG);
lower=-5*10^(2); upper=5*10^(2);
M(M<=lower) = lower;  M(M>=upper) = upper;
M = full(M);
clearvars -except M

% 
tic
[U, S, V] = svd(M); % U1=mxm S1=mxn V1=nxn 
toc    % for early: 3639.9 sec, mid: 3485.4 sec,  late: 3467.6sec
Sdiag = diag(S);
Sdiag = Sdiag.^2;

%
e = 0;
index_85 = 0;
e85 = 0.85*sum(Sdiag);
while(e < e85)
    index_85 = index_85+1;
    e = e+Sdiag(index_85);
end

%  Rank
downstream = M*V(:,1:index_85);
upstream = M'*U(:,1:index_85);      
downstream(isnan(downstream))=0;
down = sqrt(sum(downstream.^2,2));
upstream(isnan(upstream))=0;
up = sqrt(sum(upstream.^2,2));

save('result_85.mat','down', 'up','-v7.3');

%%-------- rank--------
clear all, clc;
load('PNP_name.mat');  load('result_85.mat', 'down', 'up');
c = containers.Map('KeyType','char','ValueType','double');
% Protein and gene name are identical
n = zeros(1,4);
for i = 1:length(down)
    if c.isKey(char(PNP_down_name(i)))
        c(char(PNP_down_name(i))) = (c(char(PNP_down_name(i))) + down(i)) / 2;
        n(1) = n(1)+1;
    else
        c(char(PNP_down_name(i))) = down(i);
        n(2) = n(2)+1;
    end
end

for i = 1:length(up)
    if c.isKey(char(PNP_up_name(i)))
        c(char(PNP_up_name(i))) = (c(char(PNP_up_name(i))) + up(i)) / 2;
        n(3) = n(3)+1;
    else
        c(char(PNP_up_name(i))) = up(i);
        n(4) = n(4)+1;
    end
end


keys = c.keys;
values = cell2mat(c.values);
[rank_value, sortId] = sort(values,'descend');
rank_value = rank_value';
rank_gene = cell(length(sortId),1);
for i = 1:length(sortId)
    rank_gene{i} = keys{sortId(i)};
end
rank_gene = string(rank_gene);
save('Rank_85.mat','rank_value', 'rank_gene');


File = fopen('rank_85.txt','w');
File2 = fopen('rank_85_number_for_cytoscape.txt','w');
fprintf(File2,'name\torder\n');
for i = 1:length(rank_gene)
    if i == length(rank_gene)
        fprintf(File,'%s',rank_gene(i));
        fprintf(File2,'%s\t%d', rank_gene(i), i);
    else
        fprintf(File,'%s\n',rank_gene(i));
        fprintf(File2,'%s\t%d\n', rank_gene(i), i);
    end
end
fclose(File);
fclose(File2);
