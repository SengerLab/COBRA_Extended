function extract_COBRA(r)

%extract_SBML_20120424.m

%written by the Senger Research Group, 4/24/2012, last updated 10/3/2012

%this program extracts useful information from an SBML COBRA .mat
%file so it can be stored as a spreadsheet.

%the gene names and EC numbers are added only if they exist.

fprintf(['Extracting model to spreadsheets...']);

[file1,path1]=uiputfile('reactions.txt','Select the folder to put the reactions.txt file');
[file2,path2]=uiputfile('metabolites.txt','Select the folder to put the metabolites.txt file');

fid1=fopen(fullfile(path1,file1),'w');
fid2=fopen(fullfile(path2,file2),'w');

%write headers
fprintf(fid1,['ID' 9 'Rxn ID' 9 'Rxn Name' 9 'Sub-System' 9 'Gene' 9 'EC' 9 'LB' 9 'UB' 9 'Obj' 9 'Reaction' 9 'Reaction' '\n']);
fprintf(fid2,['ID' 9 'Met Name' 9 'Metabolite' '\n']);

%extract the gene matrix
r.gm=full(r.rxnGeneMat);

for i=1:size(r.rxns,1)
    
    rxnName=r.rxnNames{i,1};
    rxnID=r.rxns{i,1};
    sub=r.subSystems{i,1};
    if size(sub,2)==0
        sub='N/A';
    end
   
    lb=r.lb(i,1);
    ub=r.ub(i,1);
    obj=r.c(i,1);
    
    EC=r.rxnECNumbers{i,1};
    if size(EC,2)==0
        EC='Unassigned';
    end
    
    %reading the gene matrix
    gene_temp=find(r.gm(i,:)==1);
    if size(gene_temp,2)>0
        
        gene_temp2=[];
        for j=1:size(gene_temp,2)
            
            if j==1
                gene_temp2=r.genes{gene_temp(1,j),1};
            else
                gene_temp2=[gene_temp2 '|' r.genes{gene_temp(1,j),1}];
            end
        end
    else
        gene_temp2='Unassigned';
    end
    gene=gene_temp2;
    clear r.gm    
    
    %building the reaction elements
    s_temp=[];
    r_temp=r.S(:,i);
    c_temp=find(r_temp~=0);
    for j=1:size(c_temp,1)
        s_temp(j,1)=r.S(c_temp(j,1),i);
        cmpd{j,1}=r.mets{c_temp(j,1),1};
        cmpd2{j,1}=r.metNames{c_temp(j,1),1};
    end
    s_temp=full(s_temp);
    
    %bulding reactants and products
    reactants=[];
    reactants2=[];
    products=[];
    products2=[];
    
    for j=1:size(s_temp,1)
        if s_temp(j,1)<0
            if size(reactants,2)==0
                
                if s_temp(j,1)==-1
                    reactants=[cmpd{j,1}];
                    reactants2=[cmpd2{j,1}];
                else
                    reactants=[num2str(-1.*s_temp(j,1)) ' ' cmpd{j,1}];
                    reactants2=[num2str(-1.*s_temp(j,1)) ' ' cmpd2{j,1}];
                end
            else
                
                if s_temp(j,1)==-1
                    reactants=[reactants ' + ' cmpd{j,1}];
                    reactants2=[reactants2 ' + ' cmpd2{j,1}];
                else
                    reactants=[reactants ' + ' num2str(-1.*s_temp(j,1)) ' ' cmpd{j,1}];
                    reactants2=[reactants2 ' + ' num2str(-1.*s_temp(j,1)) ' ' cmpd2{j,1}];
                end
            end
        else
            if size(products,2)==0
                
                if s_temp(j,1)==1
                    products=[cmpd{j,1}];
                    products2=[cmpd2{j,1}];
                else
                    products=[num2str(s_temp(j,1)) ' ' cmpd{j,1}];
                    products2=[num2str(s_temp(j,1)) ' ' cmpd2{j,1}];
                end
                
            else
                
                if s_temp(j,1)==1
                    products=[products ' + ' cmpd{j,1}];
                    products2=[products2 ' + ' cmpd2{j,1}];
                else
                    products=[products ' + ' num2str(s_temp(j,1)) ' ' cmpd{j,1}];
                    products2=[products2 ' + ' num2str(s_temp(j,1)) ' ' cmpd2{j,1}];
                end
            end
        end
    end
    
    %bulding the reaction string
    marker=' <=> ';
    
    if lb==0
        marker=' => ';
    end
    
    if ub==0
        marker=' <= ';
    end
    
    if lb==0 && ub==0
        marker=' <=> ';
    end
    
%     if size(reactants,2)==0
%         reaction=[products ' <=> '];
%     else
        reaction=[reactants marker products];
        reaction2=[reactants2 marker products2];
%     end
        
    %print to file   
    
    fprintf(fid1,[num2str(i) 9 rxnID 9 rxnName 9 sub 9 gene 9 EC 9 num2str(lb) 9 num2str(ub) 9 num2str(obj) 9 reaction 9 reaction2 '\n']);
    fprintf(['working on reaction ' num2str(i) ' of ' num2str(size(r.rxns,1)) '\n']);
    
end

fid1=fclose(fid1);
   
for i=1:size(r.mets,1)
    
    met=r.mets{i,1};
    met_name=r.metNames{i,1};
    
    fprintf(fid2,[num2str(i) 9 met 9 met_name '\n']);
    fprintf(['working on metabolite ' num2str(i) ' of ' num2str(size(r.mets,1)) '\n']);
    
end

fid2=fclose(fid2);

fprintf(['done!' '\n']);
    