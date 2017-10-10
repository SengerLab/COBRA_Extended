function r=import_rxn_notation(r,import_data,rebuild)

%written by the Senger Research Group, 10/7/2012

%this simple function allows one to update reaction notation data in excel
%and directly import to a cobra model in matalb.

%the following columns are copied/pasted from Excel using the "paste to workspace" 
%command and are stored as a cell vector normally called "A_pastespecial".
%Line 1 contains column headers.

%When calling this function, use the name of this dataset for "import_data"

%col1: ID
%col2: Rxn ID --- the field rxns
%col3: Rxn Name --- the field rxnNames
%col4: Sub-system --- the field subSystems
%col5: Gene --- the field genes
%col6: EC --- the field rxnECNumbers

%Note: for the genes field, in the Excel spreadsheet, multiple genes for a
%single reaction must be separated by a "|" delimiter (no spaces)

%the function also calls for a "rebuild" option. This is a binary variable.
%(1) - the genes are updated and the gene reaction matrix and rules lists
%are re-built from scratch.
%(0) - the genes, gene reaction matrix, and rules field are un-altered

%% Parsing "import data" and updating notation

fprintf(['Parsing the import data and updating notation...']);

%checking line 1 (contains headers)
import_temp=import_data{1,1};
d1=findstr(9,import_temp);

if size(d1,2)~=5
    
    fprintf(['Error. The import data should contain: ID, Rxn ID, Rxn Name, Sub-System, Gene, EC fields sparated by tabs.' '\n']);
    return
end

import_data(1,:)=[];

for i=1:size(import_data,1)
    
    import_temp=import_data{i,1};
    
    d1=findstr(9,import_temp);
    
    %updating notation
    r.rxns{i,1}=import_temp(1,d1(1,1)+1:d1(1,2)-1);
    r.rxnNames{i,1}=import_temp(1,d1(1,2)+1:d1(1,3)-1);
    r.subSystems{i,1}=import_temp(1,d1(1,3)+1:d1(1,4)-1);
    r.genes2{i,1}=import_temp(1,d1(1,4)+1:d1(1,5)-1);
    r.rxnECNumbers{i,1}=import_temp(1,d1(1,5)+1:size(import_temp,2));
    
end
 
fprintf(['Done!' '\n']);

%% Re-building the gene list, gene ration matrix, and rules fields
if rebuild==1
    
    fprintf(['Rebuilding the genes field...'])
    %rebuilding the genes list
    r.genes={};
    
    for i=1:size(r.genes2,1)
        
        %look for multiple genes
        d1=findstr('|',r.genes2{i,1});
        
        if size(d1,2)==0 %only one gene
            
            gene_temp=r.genes2{i,1};
            
            if size(strmatch(gene_temp,r.genes,'exact'),1)==0
                
                %add to the new list
                r.genes{size(r.genes,1)+1,1}=gene_temp;
                
            end
            
        end
        
        if size(d1,2)==1 %two genes present
            
            gene_temp_line=r.genes2{i,1};
            
            gene_temp1=gene_temp_line(1,1:d1(1,1)-1);
            
            gene_temp2=gene_temp_line(1,d1(1,1)+1:size(gene_temp_line,2));
            
            %add to the list if not already present
            if size(strmatch(gene_temp1,r.genes,'exact'),1)==0
                r.genes{size(r.genes,1)+1,1}=gene_temp1;
            end
            
            %add to the list if not already present
            if size(strmatch(gene_temp2,r.genes,'exact'),1)==0
                r.genes{size(r.genes,1)+1,1}=gene_temp2;
            end
        
        end
        
        if size(d1,2)>1 %more than two genes present
            
            gene_temp_line=r.genes2{i,1};
            
            for j=1:size(d1,2)
        
                if j==1
                    
                    %extract the first gene
                    gene_temp=gene_temp_line(1,1:d1(1,j)-1);
                    
                    %add to the list if not already present
                    if size(strmatch(gene_temp,r.genes,'exact'),1)==0
                        r.genes{size(r.genes,1)+1,1}=gene_temp;
                    end
                    
                end
                
                if j==size(d1,2)
                    
                    %extract the second to last gene
                    gene_temp1=gene_temp_line(1,d1(1,size(d1,2)-1)+1:d1(1,size(d1,2))-1);
                    
                    %extract the last gene
                    gene_temp2=gene_temp_line(1,d1(1,size(d1,2))+1:size(gene_temp_line,2));
                    
                    %add to the list if not already present
                    if size(strmatch(gene_temp1,r.genes,'exact'),1)==0
                        r.genes{size(r.genes,1)+1,1}=gene_temp1;
                    end
                    
                    %add to the list if not already present
                    if size(strmatch(gene_temp2,r.genes,'exact'),1)==0
                        r.genes{size(r.genes,1)+1,1}=gene_temp2;
                    end
                    
                end
        
                if j~=1 && j~=size(d1,2)
                    
                    %extract the middle genes
                    gene_temp=gene_temp_line(1,d1(1,j-1)+1:d1(1,j)-1);
                    
                    %add to the list if not already present
                    if size(strmatch(gene_temp,r.genes,'exact'),1)==0
                        r.genes{size(r.genes,1)+1,1}=gene_temp;
                    end
                    
                end
                
            end
            
        end
    end
    
    %sort the new genes list (ascending)
    r.genes=sortrows(r.genes,1);
    fprintf(['Done!' '\n']);
    
    %% Re-building the gene reaction matrix and the rules fields
    fprintf(['Rebuilding the gene reaction matrix and rules...']);
    
    %reset the matrix and rules vector
    r.rxnGeneMat=full(zeros(size(r.rxns,1),size(r.genes,1)));
    r.grRules={};
    r.rules={};
    
    for i=1:size(r.genes2,1)
        
        %look for multiple genes
        d1=findstr('|',r.genes2{i,1});
        
        if size(d1,2)==0 %only one gene
            
            gene_temp=r.genes2{i,1};
            
            %find in r.genes
            a=strmatch(gene_temp,r.genes,'exact');
            
            %create rules and update the gene reaction matrix
            r.rules{i,1}=['x(' num2str(a) ')'];
            r.grRules{i,1}=gene_temp;
            r.rxnGeneMat(i,a)=1;
            
        else  %multiple genes present
            
            gene_temp_line=r.genes2{i,1};
            rules_temp=[];
            grRules_temp=[];
            
            for j=1:size(d1,2)+1    %for all genes present
                
                if j==1     %gene1|...
                    
                    %extract the first gene
                    gene_temp=gene_temp_line(1,1:d1(1,j)-1);
                    
                    %find in r.genes
                    a=strmatch(gene_temp,r.genes,'exact');
                    
                    %create rules and update the gene reaction matrix
                    rules_temp=['x(' num2str(a) ')|'];
                    grRules_temp=[gene_temp '|'];
                    r.rxnGeneMat(i,a)=1;
                    
                end
                
                if j==size(d1,2)+1     %...|geneN
                    
                    %extract the last gene
                    gene_temp=gene_temp_line(1,d1(1,size(d1,2))+1:size(gene_temp_line,2));
                    
                    %find in r.genes
                    a=strmatch(gene_temp,r.genes,'exact');
                    
                    %create rules and update the gene reaction matrix
                    rules_temp=[rules_temp 'x(' num2str(a) ')'];
                    grRules_temp=[grRules_temp gene_temp];
                    r.rxnGeneMat(i,a)=1;
                    
                end
                
                if j~=1 && j~=size(d1,2)+1
                    
                    %extract the middle genes
                    gene_temp=gene_temp_line(1,d1(1,j-1)+1:d1(1,j)-1);
                    
                    %find in r.genes
                    a=strmatch(gene_temp,r.genes,'exact');
                    
                    %create rules and update the gene reaction matrix
                    rules_temp=[rules_temp 'x(' num2str(a) ')|'];
                    grRules_temp=[grRules_temp gene_temp '|'];
                    r.rxnGeneMat(i,a)=1;
                    
                end
                
                r.rules{i,1}=rules_temp;
                r.grRules{i,1}=grRules_temp;
                
            end
        end
    end
    
    r.rxnGeneMat=sparse(r.rxnGeneMat);
    r=rmfield(r,'rxnGeneMat2');
    r=rmfield(r,'genes2');
    
else
    fprintf(['Rebuilding gene reaction matrix and rules bypassed...']);
end

fprintf(['Done!' '\n']);
fprintf(['Manually save the workspace.' '\n']);
end












