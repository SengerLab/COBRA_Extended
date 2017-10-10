function model=flux_ratio(model,top_rx,bottom_rx,top_dir,bottom_dir,ratio)

%written by the Senger Research Group, 7/21/2012, last updated 9/26/2012

%this function constrains flux ratios based on the concept that multiple
%enzymes are competing for the same limited metabolite.  

%the flux ratio is constrained for a specified (top) reaction over all
%competing reactions.

%example: A -> B (rxn1)
%         A -> C (rxn2)

%define the ratio based on rxn1 (e.g., ensure half of the total flux goes through rxn1):

%flux(rxn1)/[flux(rxn1 + flux(rxn2)] = 0.5

%note the reaction in the numerator of the ratio must also appear in the
%denominator

%Inputs:
%model: a cobra-style model
%top_rxn: the numerical index in the "rxns" variable of the reaction in the
%numerator of the ratio -- example: 12
%bottom_rx: a row vector of all indices of reactions in the denominator --
%example: [10 54 67]
%top_dir: 1 = top_rxn is written in the forward direction; -1 = written in
%the reverse direction -- example: 1
%bottom_dir: row vector of 1 or -1 values specifying the direction for each
%reaction in the bottom_rx vector -- example: [1 -1 -1]
%ratio: the numerical value of the flux ratio as defined above -- example:
%0.5

%Outputs:
%model: cobra-style model with added flux ratio constraint

%check if inputs are sized correctly
if size(top_rx) ~= size(top_dir)
    fprintf(['Error. The top reaction and direction do not agree.' '\n']);
    return
end

if size(bottom_rx) ~= size(bottom_dir)
    fprintf(['Error. The bottom reactions and directions do not agree.' '\n']);
    return
end

%check that the top reaction is also in the bottom reaction list
a=ismember(top_rx,bottom_rx);
for i=1:size(a,2)
    if a(1,i)==0
        fprintf(['Error. The top reaction must also appear in the bottom reaction vector.' '\n']);
        return
    end
end

%add new constraint row to stoichiometric matrix
new_row=size(model.S,1)+1;
model.S(new_row,:)=zeros(1,size(model.S,2));

%modify the row to constrain the flux ratio
for i=1:size(top_rx,2)
    model.S(new_row,top_rx(1,i))=top_dir(1,i).*(1-ratio);
end

for i=1:size(bottom_rx,2)
    if ismember(bottom_rx(1,i),top_rx)==0
        model.S(new_row,bottom_rx(1,i))=-1.*bottom_dir(1,i).*ratio;
    end
end

%add information to the metabolite fields
model.mets{new_row}=['ratio_' model.rxns{top_rx} ':' model.rxns{bottom_rx}];
model.metNames{new_row}=['ratio_' model.rxns{top_rx} ':' model.rxns{bottom_rx}];

%updating additional fields - different models contain different fields
%(this handles them all)
fnames=fieldnames(model);

%S and mets are already processed - remove them from the existing list
fnames(strmatch('mets',fnames,'exact'),:)=[];
fnames(strmatch('S',fnames,'exact'),:)=[];

for i=1:size(fnames,1)
    
    fget=getfield(model,fnames{i,1});
    
    if iscell(fget)==1
        
        if size(fget,1)==new_row-1
            fget{new_row,1}='';
        end
        
    end
    
    if isnumeric(fget)==1
        
        if size(fget,1)==new_row-1
            fget(new_row,1)=0;
        end
    end
    
    model=setfield(model,fnames{i,1},fget);
    
end
        

end