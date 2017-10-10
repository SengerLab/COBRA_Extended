function newModel=spf(model, spf_value)
%Speficies the specific proton flux of a COBRA formatted genome-scale model
%
%newModel=spf(model, reactions, protons, directions, spf_value)
%
%written by the Senger Research Group, 8/1/2009
%last updated by Ryan Senger, 2/15/2014
%
%INPUTS
%1. model           COBRA formatted genome-scale model
%2. reactions       row vector of reaction indices that transport protons
%3. protons         row vector of number of protons in each reaction
%4. directions      row vector of directions (1)=import, (-1)=export
%5. spf_value       numerical value (<0)=proton export
%
%OUPUTS
%1. newModel        COBRA formatted model with constrained spf

%Identify membrane transport reactions
[transRxns, nonTransRxns, membTransRxns, periTransRxns]=findTransRxnsModif(model);

%Ammend the stoichiometric matrix
newModel=model;
newModel.S(:,end+1)=0;
newModel.S(end+1,:)=0;
newModel.S(end,end)=-spf_value;     %insert the spf value

%Ammend other important fields (may be edited as needed)
newModel.rxns{end+1,1}='spf';
newModel.mets{end+1,1}='spf';
newModel.rev(end+1,1)=1;
newModel.lb(end+1,1)=1; %these constraints are required for the spf
newModel.ub(end+1,1)=1;
newModel.c(end+1,1)=0;
newModel.metCharge(end+1,1)=0;
newModel.rules{end+1,1}={};
newModel.genes{end+1,1}='spf';
newModel.rxnGeneMat(:,end+1)=0;
newModel.rxnGeneMat(end+1,:)=0;
newModel.grRules{end+1,1}='spf';
newModel.subSystems{end+1,1}='spf';
newModel.confidenceScores{end+1,1}={};
newModel.rxnReferences{end+1,1}={};
newModel.rxnECNumbers{end+1,1}='spf';
newModel.rxnNotes{end+1,1}='spf';
newModel.rxnNames{end+1,1}='spf';
newModel.metFormulas{end+1,1}={};
newModel.metChEBIID{end+1,1}={};
newModel.metKEGGID{end+1,1}={};
newModel.metPubChemID{end+1,1}={};
newModel.metInChIString{end+1,1}={};
newModel.b(end+1,1)=0;

%Locating indices of H[e], H[c], H[p] - also finds h[], H+[], h+[]
a1=regexpi(newModel.mets,'^h(+)?\[e\]');
b1=regexpi(newModel.mets,'^h(+)?\[p\]');
c1=regexpi(newModel.mets,'^h(+)?\[c\]');

a2=cellfun(@isempty,a1);
b2=cellfun(@isempty,b1);
c2=cellfun(@isempty,c1);

he=find(a2==0);
hp=find(b2==0);
hc=find(c2==0);

if size(he,2)>1 || size(hp,2)>1 || size(hc,2)>1 || size(he,2)==0 || size(hc,2)==0
    
    fprintf(['Error in identifying H[e], H[p], or H[c]. Abort.' '\n']);
    return
end

no_hp=0;    %dealing with the possible periplasm compartment

if size(hp,1)==0
    
    no_hp=1;
    
end

for i=1:size(membTransRxns,1)
    
    if membTransRxns(i,1)==1    %it is a transport reaction
        
        mets=find(newModel.S(:,i));     %create the vector of met indices
        
        stoich=[];
        for j=1:size(mets,1)
            
            stoich(j,1)=newModel.S(mets(j,1),i);    %create the vector of stoichiometries
            
        end
        
        if no_hp==1     %dealing with H[e] and H[c]
            
            he_stoich=0;
            hc_stoich=0;
            
            he_find=find(mets==he);
            
            if size(he_find,1)>0
                
                he_stoich=stoich(he_find,1);
                
            end
            
            hc_find=find(mets==hc);
            
            if size(hc_find,1)>0
                
                hc_stoich=stoich(hc_find,1);
                
            end
            
            if abs(he_stoich) > abs(hc_stoich)
                
                newModel.S(end,i)=-he_stoich;
                
            else
                
                newModel.S(end,i)=hc_stoich;
                
            end
            
        end
        
        if no_hp==0     %dealing with H[e] and H[p]
            
            he_stoich=0;
            hp_stoich=0;
            
            he_find=find(mets==he);
            
            if size(he_find,1)>0
                
                he_stoich=stoich(he_find,1);
                
            end
            
            hp_find=find(mets==hp);
            
            if size(hp_find,1)>0
                
                hp_stoich=stoich(hp_find,1);
                
            end
            
            if abs(he_stoich) > abs(hp_stoich)
                
                newModel.S(end,i)=-he_stoich;
                
            else
                
                newModel.S(end,i)=hp_stoich;
                
            end
            
        end
        
    end
    
end
            









