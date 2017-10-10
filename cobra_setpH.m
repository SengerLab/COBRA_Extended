function model=cobra_setpH(modeli,pH,skipFlag,stoichFlag)

%written by the Senger Research Group, 2/2/2011

%cobra_setpH modifies the stoichiometric matrix of an input model.
%The metabolic charges (in metCharge) are changed according to the input
%pH. The stoichiometric matrix is then adjusted accordingly by correcting
%the rows corresponding to H+ for all the reactions through charge balancing.
%
%model=cobra_setpH(modeli,pH,stoichFlag)
%
%The input model (modeli) must have the following additional fields:
%  apkas          Array containing acidic pka's of each metabolite
%  bpkas          Array containing basic pka's of each metabolite
%  formalCharge   Vector containing the formal charge of each metabolite
%  metFormulas    Array containing the metabolic formulas (needed for
%                 stoichCheck)
%pH               pH vector (one value for each compartment in alphabetical
%                 order) i.e. pH = [pH for 'c' compartment ; pH for 'e'
%                 compartment]
%skipFlag         input true if skipCheck (indicates which reactions have
%                 been not been proton-balanced) is needed in output
%                 model, false otherwise. Default true.
%stoichFlag       input true if stoichCheck (stores stoichiometric error
%                 values in terms of H+ for all reactions) is needed in 
%                 output model, false otherwise. Default true.
%
%The output model will contain the following additional fields:
%  metCharge      Vector containing the average charge of each metabolite
%  pH             Vector indicating what pH each compartment is set to
%skipped          Array containing rxn ID's of all reactions that could not
%                 be modified by charge balancing
%stoichCheck      Array containing stoichiometric error values in terms of
%                 H+ for all reactions
%
%This function relies on the calc_avgCharge function


%specify how many decimal points H+ stoichiometries must round to
H_round=2;

if nargin<3
    skipFlag=true;
end
if nargin<4
    stoichFlag=true;
end


model=modeli;
model.pH=pH;


%create the compartments (cmpt) array
for row=1:size(model.mets,1)
    met_temp=model.mets{row,1};
    bracket_finder=strfind(met_temp,'[');
    cmpt{row,1}=met_temp(1,bracket_finder+1);
end
model_cmpt=unique(cmpt);    %holds which compartments are in the model

%find where H+ is located
H_row=strmatch('H',model.metFormulas,'exact');
H_cmpt=cmpt(H_row);

%check if the model has the correct number of H+ entries as a metabolite
H_check=size(model_cmpt,1)-size(H_row,1);
if H_check>0
    
    %add H+ as a metabolite in compartments it is missing in
    for a=1:size(model_cmpt,1)
        cmpt_check=strmatch(model_cmpt{a,1},H_cmpt);
        if size(cmpt_check,1)==0
            fprintf(['H+ not found as a compound in compartment ' model_cmpt{a,1} '. New compound data for H+ in this compartment created in the model.\n'])
            H_row_temp=size(model.mets,1)+1;
            row_temp=size(H_row,1)+1;
            H_row(row_temp,1)=H_row_temp;
            H_cmpt{row_temp}=model_cmpt{a,1};
            model.mets{H_row_temp,1}='cid0000068';
            model.metFormula(H_row_temp,1)='H';
            model.S(H_row_temp,:)=zeros(1,size(model.S,2));
            model.metCharge(H_row_temp,1)=1;
            model.formalCharge(H_row_temp,1)=1;
            model.apkas(H_row_temp,1)=NaN;
            model.bpkas(H_row_temp,1)=NaN;
        end
    end
    
    %double check if the additions make sense
    H_check=size(model_cmpt,1)-size(H_row,1);
    if H_check~=0
        error('Error: H+ entries do not make sense!')
    end
    
%if there are too many entries for H+, display an error message
elseif H_check<0
    error('Error: Too many compound entries for H+ detected!')
end

if size(pH,1)~=size(model_cmpt,1)
    error('Error: size of pH vector does not match amount of compartments in model!')
end



%%%creating the metCharge vector

%prepare to make usable vectors without NaN entries
nan_mets=isnan(model.formalCharge);

bar=waitbar(0,'Calculating the average charge for each compound...');
for row=1:size(model.mets,1)
    
    %retrieve pka values from apkas and bpkas arrays
    a_pkas_temp=model.apkas{row,1};
    b_pkas_temp=model.bpkas{row,1};
    
    %retrieve formal charge data
    formalcharge_temp=model.formalCharge(row,1);
    
    %find the [H+] in the compartment the metabolite is in
    H_row_temp=strmatch(cmpt{row,1},H_cmpt);
    pH_temp=pH(H_row_temp);
    
    %calculate the metabolite's average charge
    average_charge_temp=calc_avgCharge(pH_temp,formalcharge_temp,a_pkas_temp,b_pkas_temp);
    model.metCharge(row,1)=average_charge_temp;
    
    %prepare the H_stoich vector to be used in stoichiometric checking
    if stoichFlag==1
        formula_temp=model.metFormulas{row,1};
        formula_temp=[formula_temp ' '];
        H_finder=regexp(formula_temp,'H[^a-z]');
        stoich_finder=regexp(formula_temp,'H[0-9]+','match');
        if size(H_finder,2)==1 && size(stoich_finder,2)==0
            H_stoich_temp=1;
        elseif size(H_finder,2)==1 && size(stoich_finder{1,1},2)>1
            H_stoich_temp=str2num(stoich_finder{1,1}(1,2:end));
        else
            H_stoich_temp=0;
        end
        H_stoich(row,1)=H_stoich_temp;
        
        %prepare to ignore metabolites whose average charges can't be
        %calculated
        if isnan(average_charge_temp);  %if size(formula_temp,2)>3 && size(a_pkas_temp,2)==0 && size(b_pkas_temp,2)==0   7/10/2011 (see calc_avgCharge)
            nan_mets(row)=1;
        end
        
    end
    
    %printing the program's progress
    waitbar(row/size(model.mets,1),bar);
    
    %renewing variables
    components=[];
    Dcharges=[];
    
end

close(bar)

%Make usable vectors without NaN entries
metCharge_nanfree=model.metCharge;
metCharge_nanfree(nan_mets)=0;
Dcharge=model.metCharge-model.formalCharge;
Dcharge_nanfree=Dcharge;
Dcharge_nanfree(nan_mets)=0;


%%%Modifying the H+ rows in the stoichiometric matrix
bar=waitbar(0,'Modifying H+ rows in the stoichiometric matrix...');
skipCheck={};
column=1;
S_end=size(model.S,2);
while column<=S_end                    %for column=1:size(model.S,2)
    
    %extracting relevant data for current rxn
    stoichcolumn_temp=full(model.S(:,column));
    involved_rows=find(model.S(:,column)~=0);
    involved_cmpts=cmpt(involved_rows);
    unique_cmpts=unique(involved_cmpts);
    nonreal_check=abs(stoichcolumn_temp).*nan_mets;
    nonreal_check=sum(nonreal_check);
    
    %prepare for specificity (if available)
    rxn_name=model.rxns{column,1};
    specific=strfind(rxn_name,'_S');
    
    if nonreal_check==0     %if any of the compounds involved in the reaction don't have a charge or formula listed, then the reaction must be skipped
        
        %exchange reactions do not need to be proton balanced
        if size(involved_rows,1)==1
            
        %if only one compartment is involved (only one H+ row needs to be
        %modified)
        elseif size(unique_cmpts,1)==1
            H_row_temp=strmatch(cmpt{involved_rows(1,1),1},H_cmpt);
            H_row_temp=H_row(H_row_temp);
            H_add=stoichcolumn_temp.*metCharge_nanfree;    %Dcharge_nanfree;
            H_add=sum(H_add);
            model.S(H_row_temp,column)=model.S(H_row_temp,column)-H_add;
        
            
        %%%%%if multiple compartments are involved
        
        %proton-balancing according to specificity (if available)
            
        elseif size(specific,2)==1
                stoichcolumn_hfree=stoichcolumn_temp;
                stoichcolumn_hfree(H_row,1)=0;
                metrows_temp=find(stoichcolumn_hfree~=0);
                metcmpts_temp=cmpt(metrows_temp);
                umetcmpts_temp=unique(metcmpts_temp);
                species_charge=str2double(rxn_name(1,specific(1,1)+2:end));  %find the charge that the reaction is specific to
                Dspeciesavg=metCharge_nanfree-species_charge;        %find the difference between average charges and species charge
                for a=1:size(umetcmpts_temp,1)
                    cmptstoich_temp=stoichcolumn_hfree;
                    in_cmpt=strcmp(cmpt,umetcmpts_temp{a,1});
                    cmptstoich_temp(in_cmpt==0)=0;                  %exclude all metabolites not in the current compartment
                    H_row_temp=strmatch(umetcmpts_temp{a,1},H_cmpt);
                    H_row_temp=H_row(H_row_temp);
                    H_add=cmptstoich_temp.*(Dspeciesavg);
                    H_add=sum(H_add);
                    model.S(H_row_temp,column)=-1*H_add;
                end
                
        else
            
%             H_involved=stoichcolumn_temp(H_row);
%             if nnz(H_involved)>2
%                 error(['Error: Protons involved in more than 2 compartments in reaction ' model.rxn{column,1} '! \n'])
%             end
%             
%             %calculating charge_transfer
%             charge_transfer=[];
%             for a=1:size(unique_cmpts,1)
%                 cmpt_temp=unique_cmpts{a,1};
%                 column_check=strcmp(cmpt_temp,cmpt);
%                 column_temp=stoichcolumn_temp.*column_check;  %exclude metabolites from other compartments
%                 transfer_temp=column_temp.*model.modelCharge;
%                 charge_transfer(a,1)=sum(transfer_temp);
%             end
%             if sum(charge_transfer)~=0
%                 error(['Original model not properly charge balanced at reaction ' model.rxns{column,1} ' (column ' num2str(column) '). Error of ' num2str(sum(charge_transfer)) ' protons encountered.'])
%             end
            


            %%%Method of no proton specificity
            
            %create a new rxn for reverse rxn if it occurs
%             model.S(H_row,column)=0;    %no need to delete previous H+ stoichiometry
%             stoichcolumn_temp=full(model.S(:,column));
%             pka_check=stoichcolumn_temp.*Dcharge_nanfree;
%             pka_check=find(pka_check~=0);
%             if model.rev(column,1)==1 && size(pka_check,1)~=0 && size(specific,2)==0
%                 column_new=size(model.S,2)+1;
%                 model.S(:,column_new)=-1.*model.S(:,column);
%                 model.lb(column_new,1)=0;
%                 model.ub(column_new,1)=-1*model.lb(column,1);
%                 model.rev(column_new,1)=0;
%                 model.c(column_new,1)=0;
%                 model.rxns{column_new,1}=['rev_' model.rxns{column,1}];
%                 model.rxnNames{column_new,1}=['rev_' model.rxnNames{column,1}];
%                 model.subSystems{column_new,1}=model.subSystems{column,1};
%                 model.grRules{column_new,1}=model.grRules{column,1};
%                 model.rxnGeneMat(column_new,1)=model.rxnGeneMat(column,1);
%                 model.lb(column,1)=0;
%                 model.rev(column,1)=0;
%                 S_end=S_end+1;
%             end
            
            %reverse the reaction to the forward direction if necessary
            if model.rev(column,1)==0 && model.lb(column,1)<0
                stoichcolumn_temp=-1.*stoichcolumn_temp;
                model.S(:,column)=stoichcolumn_temp;
                lb_temp=model.lb(column,1);
                ub_temp=model.ub(column,1);
                model.lb(column,1)=-1*ub_temp;
                model.ub(column,1)=-1*lb_temp;
            end
            
            %proton balance the transport reaction
            substraterows_temp=find(stoichcolumn_temp<0);
            productrows_temp=find(stoichcolumn_temp>0);
            productcmpts_temp=cmpt(productrows_temp);
            uproductcmpts_temp=unique(productcmpts_temp);
            substratecmpts_temp=cmpt(substraterows_temp);
            usubstrate_cmpts_temp=unique(substratecmpts_temp);
            product_no=size(productrows_temp,1);
            while (size(uproductcmpts_temp,1)*product_no)~=0
                if size(uproductcmpts_temp,1)==1  %Condition: all products are in one compartment
                    H_row_temp=strmatch(uproductcmpts_temp{1,1},H_cmpt);
                    H_row_temp=H_row(H_row_temp);
                    H_add=stoichcolumn_temp.*metCharge_nanfree;    %Dcharge_nanfree;
                    H_add=sum(H_add);
                    model.S(H_row_temp,column)=model.S(H_row_temp,column)-H_add;
                    uproductcmpts_temp={};
                else
                    cpdrow_match=strmatch(model.metNames{productrows_temp(product_no,1),1},model.metNames(substraterows_temp));
                    productcmpt_temp=cmpt{productrows_temp(product_no,1),1};
                    substratecmpt_check=strmatch(productcmpt_temp,substratecmpts_temp);
                    if size(cpdrow_match,1)==1   %Sub-rxn: one compound traveling across the membrane
                        H_row_temp=strmatch(productcmpts_temp{product_no,1},H_cmpt);
                        H_row_temp=H_row(H_row_temp);
                        H_match=(productrows_temp(product_no,1)==H_row);
                        if (stoichcolumn_temp(productrows_temp(product_no,1),1) ~= abs(stoichcolumn_temp(substraterows_temp(cpdrow_match,1),1))) && sum(H_match)~=0
                            stoich_new=[stoichcolumn_temp(productrows_temp(product_no,1),1); abs(stoichcolumn_temp(substraterows_temp(cpdrow_match,1),1))];
                            stoich_new=sort(stoich_new,1,'descend');
                            stoich_new=stoich_new(1,1);
                            productstoich_temp=stoich_new;
                            substratestoich_temp=-1*stoich_new;
                            model.S(productrows_temp(product_no,1),column)=productstoich_temp;
                            model.S(substraterows_temp(cpdrow_match,1),column)=substratestoich_temp;
                        else
                            productstoich_temp=stoichcolumn_temp(productrows_temp(product_no,1),1);
                            substratestoich_temp=stoichcolumn_temp(substraterows_temp(cpdrow_match,1),1);
                        end
                        H_add= productstoich_temp*metCharge_nanfree(productrows_temp(product_no,1)) + substratestoich_temp*metCharge_nanfree(substraterows_temp(cpdrow_match,1),1);    %Dcharge_nanfree;
                        H_add=sum(H_add);
                        model.S(H_row_temp,column)=model.S(H_row_temp,column)-H_add;
                        stoichcolumn_temp(productrows_temp(product_no,1),1)=0;
                        stoichcolumn_temp(substraterows_temp(cpdrow_match,1),1)=0;
                        substraterows_temp=find(stoichcolumn_temp<0);
                        productrows_temp=find(stoichcolumn_temp>0);
                        productcmpts_temp=cmpt(productrows_temp);
                    elseif size(cpdrow_match,1)==0 && size(substratecmpt_check,1)==0 && size(usubstrate_cmpts_temp,1)==1   %Sub-rxn: one substrate is modified before being transported
                        cpdrow_temp=productrows_temp(product_no,1);
                        othercmpt=usubstrate_cmpts_temp{1,1};
                        H_row_temp=strmatch(othercmpt,H_cmpt);
                        pH_temp=pH(H_row_temp);
                        formalCharge_temp=model.formalCharge(cpdrow_temp);
                        a_pkas_temp=model.apkas{cpdrow_temp,1};
                        b_pkas_temp=model.bpkas{cpdrow_temp,1};
                        othercmpt_charge=calc_avgCharge(pH_temp,formalCharge_temp,a_pkas_temp,b_pkas_temp);
                        H_add=stoichcolumn_temp(productrows_temp(product_no,1)) * (metCharge_nanfree(productrows_temp(product_no,1)) - othercmpt_charge);
                        model.S(H_row_temp,column)=model.S(H_row_temp,column)-H_add;
                        stoichcolumn_temp(productrows_temp(product_no,1),1)=stoichcolumn_temp(productrows_temp(product_no,1),1) * othercmpt_charge/metCharge_nanfree(productrows_temp(product_no,1));
                        productcmpts_temp{product_no}=othercmpt;
                    end
                    uproductcmpts_temp=unique(productcmpts_temp);
                    product_no=product_no-1;
                end
            end
            if size(uproductcmpts_temp,1)~=0  %Sub-rxn: compartment balance remaining mets
                involved_rows=find(stoichcolumn_temp~=0);
                involved_cmpts=cmpt(involved_rows);
                unique_cmpts=unique(involved_cmpts);
                for a=1:size(unique_cmpts,1)
                    cmpt_temp=unique_cmpts{a,1};
                    column_check=strcmp(cmpt_temp,cmpt);
                    column_temp=stoichcolumn_temp.*column_check;  %exclude metabolites from other compartments
                    H_row_temp=strmatch(cmpt_temp,H_cmpt);
                    H_row_temp=H_row(H_row_temp);
                    H_add=column_temp.*metCharge_nanfree;
                    H_add=sum(H_add);
                    model.S(H_row_temp,column)=model.S(H_row_temp,column)-H_add;
                end
            end
            
            
%             %method of copying the original model's charge transfer between
%             %compartments (charge_transfer needed)
%             for a=1:size(unique_cmpts,1)
%                 cmpt_temp=unique_cmpts{a,1};
%                 column_check=strcmp(cmpt_temp,cmpt);
%                 column_temp=stoichcolumn_temp.*column_check;  %exclude metabolites from other compartments
%                 H_row_temp=strmatch(cmpt_temp,H_cmpt);
%                 H_row_temp=H_row(H_row_temp);
%                 charge_temp=column_temp.*metCharge_nanfree;
%                 charge_temp=sum(charge_temp);
%                 model.S(H_row_temp,column)=model.S(H_row_temp,column)+(charge_transfer(a,1)-charge_temp);
%             end
%             
%             
%             %method of adding/subtracting from H row by looking at each
%             metabolite one by one
%             for a=1:size(involved_rows,1) %additions are made metabolite by metabolite
%                 met_row_temp=involved_rows(a,1);
%                 cmpt_temp=cmpt(met_row_temp,1);
%                 H_row_temp=strmatch(cmpt_temp,H_cmpt);
%                 H_row_temp=H_row(H_row_temp);
%                 H_add=stoichcolumn_temp(met_row_temp,1)*Dcharge_nanfree(met_row_temp,1);
%                 model.S(H_row_temp,column)=model.S(H_row_temp,column)-H_add;
%             end
%             
%             %method of balancing reactant and product sides for each
%             %compartment (proton transport is messy)
%             for a=1:size(unique_cmpts,1)
%                 cmpt_temp=unique_cmpts{a,1};
%                 column_check=strcmp(cmpt_temp,cmpt);
%                 column_temp=stoichcolumn_temp.*column_check;  %exclude metabolites from other compartments
%                 H_row_temp=strmatch(cmpt_temp,H_cmpt);
%                 H_row_temp=H_row(H_row_temp);
%                 H_add=column_temp.*Dcharge_nanfree;
%                 H_add=sum(H_add);
%                 model.S(H_row_temp,column)=model.S(H_row_temp,column)-H_add;
%             end
        end
        
        %round the results to the nearest 2 decimal places
        H_add=model.S(H_row,column);
        H_add=H_add*(10^H_round);
        H_add=round(H_add);
        H_add=H_add/(10^H_round);
        model.S(H_row,column)=H_add;
        
        %checking the stoichiometric balance of H in the current column of S
        if stoichFlag
            HTemp=full(model.S(:,column)).*(H_stoich+Dcharge_nanfree);
            stoichCheck(column,1)=sum(HTemp);
        end
                
        
    else
        skiprow=size(skipCheck,1)+1;
        skipCheck{skiprow,1}=model.rxns{column,1};
        skipCheck{skiprow,2}=num2str(column);
        
    end
    
    %printing the program's progress
    waitbar(column/S_end,bar);
    column=column+1;
end

close(bar)

%add optional fields to the model
if skipFlag
    model.skipCheck=skipCheck;
end
if stoichFlag
    model.stoichCheck=stoichCheck;
end


end