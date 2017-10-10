function extract_COBRAv2(model)

%written by the Senger Research Group, 4/24/2012, last updated 10/4/2013

%this program extracts useful information from an SBML COBRA .mat
%file so it can be stored as a spreadsheet.

%the update of this program produces a format that can be re-uploaded using
%the "xls2modelv2" function described in "Protocol for translating a
%genome-scale model in text format into a COBRA structure.docx"

%the gene names and EC numbers are added only if they exist.

% Last edited, 10-07-2013

fprintf(['Extracting model to spreadsheets...' '\n']);

[file1,path1]=uiputfile('reactions.txt','Select the folder to put the reactions.txt file');
[file2,path1]=uiputfile('metabolites.txt','Select the folder to put the metabolites.txt file');

fid1=fopen(fullfile(path1,file1),'w');
fid2=fopen(fullfile(path1,file2),'w');

%write headers
fprintf(fid1,['Abbreviation' 9 'Name' 9 'Reaction' 9 'GPR' 9 'Genes' 9 'Proteins' 9 'Subsystem' 9 'Reversible' 9 'Lower bound' 9 'Upper bound' 9 'Objective' 9 'Confidence Score' 9 'EC. Number' 9 'Notes' 9 'References' '\n']);
fprintf(fid2,['Abbreviation' 9 'Name' 9 'Formula (neutral)' 9 'Formula (charged)' 9 'Charge' 9 'Compartment' 9 'KEGG ID' 9 'PubChem ID' 9 'ChEBI ID' 9 'InChI string' 9 'Smiles' '\n']);

%load model

%extract the gene matrix
gm = full(model.rxnGeneMat);

for i=1:size(model.rxns,1)
    
    rxid = model.rxns{i,1};
        
    if (isfield(model, 'rxnNames')==1)
        if ~isempty(model.rxnNames{i,1})
            rxnm = model.rxnNames{i,1};
        else
            rxnm = 'Unassigned';
        end
    else
        rxnm = 'N/A';
    end
    
    %building the reaction elements
    s_temp=[];
    r_temp=model.S(:,i);
    c_temp=find(r_temp~=0);
    for j=1:size(c_temp,1)
        s_temp(j,1)=model.S(c_temp(j,1),i);
        cmpd{j,1}=model.mets{c_temp(j,1),1};
        cmpd2{j,1}=model.metNames{c_temp(j,1),1};
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
    
    lb = model.lb(i,1);
    ub = model.ub(i,1);
    obj = model.c(i,1);
    
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
    
    rxeq=[reactants marker products];
    reaction2=[reactants2 marker products2];
    
    if (isfield(model, 'grRules')==1)
        if ~isempty(model.grRules{i,1})
            grul = model.grRules{i,1};
        else
            grul = 'Unassigned';
        end
    else
        grul = 'N/A';
    end
    
    %reading the gene matrix
    gene_temp=find(gm(i,:)==1);
    if size(gene_temp,2)>0
        
        gene_temp2=[];
        for j=1:size(gene_temp,2)
            
            if j==1
                gene_temp2=model.genes{gene_temp(1,j),1};
            else
                gene_temp2=[gene_temp2 '|' model.genes{gene_temp(1,j),1}];
            end
        end
    else
        gene_temp2='Unassigned';
    end
    gene=gene_temp2;    

    if (isfield(model, 'proteins')==1)
        if ~isempty(model.proteins{i,1})
            prot = model.proteins{i,1};
        else
            prot = 'Unassigned';
        end
    else
        prot = 'N/A';
    end
    
    if (isfield(model, 'subSystems')==1)
        if ~isempty(model.subSystems{i,1})
            subs = model.subSystems{i,1};
        else
            subs = 'Unassigned';
        end
    else
        subs = 'N/A';
    end
   
    if (isfield(model, 'rev')==1)
        if ~isempty(model.rev(i,1))
            revr = model.rev(i,1);
        else
            revr = 1;
        end
    else
        revr = 1;
    end
    
    if (isfield(model, 'confidenceScores')==1)
        if ~isempty(model.confidenceScores{i,1})
            conf = model.confidenceScores{i,1};
        else
            conf = 'Unassigned';
        end
    else
        conf = 'N/A';
    end
    
    if (isfield(model, 'rxnECNumbers')==1)
        if ~isempty(model.rxnECNumbers{i,1})
            ecnb = model.rxnECNumbers{i,1};
        else
            ecnb = 'Unassigned';
        end
    else
        ecnb = 'N/A';
    end
    
    if (isfield(model, 'rxnNotes')==1)
        if ~isempty(model.rxnNotes{i,1})
            note = model.rxnNotes{i,1};
        else
            note = 'Unassigned';
        end
    else
        note = 'N/A';
    end
    
    if (isfield(model, 'rxnReferences')==1)
        if ~isempty(model.rxnReferences{i,1})
            refr = model.rxnReferences{i,1};
        else
            refr = 'Unassigned';
        end
    else
        refr = 'N/A';
    end
      
    %print to file   
    fprintf(fid1,[rxid 9 rxnm 9 rxeq 9 grul 9 gene 9 prot 9 subs 9 num2str(revr) 9 num2str(lb) 9 num2str(ub) 9 num2str(obj) 9 conf 9 ecnb 9 note 9 refr '\n']);
    fprintf(['working on reaction ' num2str(i) ' of ' num2str(size(model.rxns,1)) '\n']);
    
end

fid1=fclose(fid1);

for i=1:size(model.mets,1)
    
    mtid = model.mets{i,1};
    mtnm = model.metNames{i,1};
        
    if (isfield(model, 'metFormulasNeutral')==1)
        if ~isempty(model.metFormulasNeutral{i,1})
            frmn = model.metFormulasNeutral{i,1};
        else
            frmn = 'Unassigned';
        end
    else
        frmn = 'N/A';
    end
    
    if (isfield(model, 'metFormulas')==1)
        if ~isempty(model.metFormulas{i,1})
            frml = model.metFormulas{i,1};
        else
            frml = 'Unassigned';
        end
    else
        frml = 'N/A';
    end
    
    if (isfield(model, 'metCompartment')==1)
        if ~isempty(model.metCompartment{i,1})
            cmpt = model.metCompartment{i,1};
        else
            cmpt = 'Unassigned';
        end
    else
        cmpt = 'N/A';
    end
    
    if (isfield(model, 'metCharge')==1)
        if ~isempty(model.metCharge(i,1))
            chrg = model.metCharge(i,1);
        else
            chrg = 0;
        end
    else
        chrg = 0;
    end
    
    if (isfield(model, 'metKEGGID')==1)
        if ~isempty(model.metKEGGID{i,1})
            kegg = model.metKEGGID{i,1};
        else
            kegg = 'Unassigned';
        end
    else
        kegg = 'N/A';
    end
    
    if (isfield(model, 'metPubChemID')==1)
        if ~isempty(model.metPubChemID{i,1})
            pubc = model.metPubChemID{i,1};
        else
            pubc = 'Unassigned';
        end
    else
        pubc = 'N/A';
    end
    
    if (isfield(model, 'metChEBIID')==1)
        if ~isempty(model.metChEBIID{i,1})
            cheb = model.metChEBIID{i,1};
        else
            cheb = 'Unassigned';
        end
    else
        cheb = 'N/A';
    end
    
    if (isfield(model, 'metInChIString')==1)
        if ~isempty(model.metInChIString{i,1})
            inci = model.metInChIString{i,1};
        else
            inci = 'Unassigned';
        end
    else
        inci = 'N/A';
    end
    
    if (isfield(model, 'metSmiles')==1)
        if ~isempty(model.metSmiles{i,1})
            smle = model.metSmiles{i,1};
        else
            smle = 'Unassigned';
        end
    else
        smle = 'N/A';
    end
    
    fprintf(fid2,[mtid 9 mtnm 9 frmn 9 frml 9 num2str(chrg) 9 cmpt 9 kegg 9 pubc 9 cheb 9 inci 9 smle '\n']);
    fprintf(['working on metabolite ' num2str(i) ' of ' num2str(size(model.mets,1)) '\n']);
    
end

fid2=fclose(fid2);

fprintf(['done!' '\n']);
    