function model=add_exchange(model,metabolite,lb,ub)

%written by the Senger Research Group, 5/23/2012; last updated 10/22/2017

%this fucntion is an easy way to add an additional exchange reaction to a
%model.

%Inputs:
%model_in: The COBRA style model
%metabolite: The metabolite to be exchanged - input the corresponding row
%in the "mets" field.
%lb: the lower bound of the exchange reaction (set to 0 for export only)
%ub: the upper bound of the exchange reaction (set to 0 for import only)


%add the new column and exchange reaction
new_col=size(model.S,2)+1;
model.S(:,new_col)=zeros(1,size(model.S,1));

model.S(metabolite,new_col)=-1;

%add the exchange reaction to the rxns list
model.rxns{new_col,1}=['EX_' model.mets{metabolite,1}];

%add the constraints
model.lb(new_col,1)=lb;
model.ub(new_col,1)=ub;

%add the objective field (set to zero by default)
model.c(new_col,1)=0;

%reversibility
if lb<0 && ub>0
    rev_temp=1;
else
    rev_temp=0;
end
if isfield(model,'rev')
    model.rev(new_col,1)=rev_temp;
end

%subsystesm
if isfield(model,'subSystems')
    model.subSystems{new_col,1}='Exchange (added)';
end

%reaction-gene matrix
if isfield(model,'rxnGeneMat')
    model.rxnGeneMat(new_col,:)=0;
end

%updating additional fields - different models contain different fields
%(this handles them all)
fnames=fieldnames(model);

%S, rxns, c, lb, ub, rev, sub-systems are already processed - remove them from the existing list
fnames(strmatch('rxns',fnames,'exact'),:)=[];
fnames(strmatch('S',fnames,'exact'),:)=[];
fnames(strmatch('c',fnames,'exact'),:)=[];
fnames(strmatch('lb',fnames,'exact'),:)=[];
fnames(strmatch('ub',fnames,'exact'),:)=[];
fnames(strmatch('rev',fnames,'exact'),:)=[];
fnames(strmatch('subSystems',fnames,'exact'),:)=[];
fnames(strmatch('rxnGeneMat',fnames,'exact'),:)=[];

if size(fnames,1)>0
    
    for i=1:size(fnames,1)
        
        fget=getfield(model,fnames{i,1});
        
        if iscell(fget)==1
            
            if size(fget,1)==new_col-1
                fget{new_col,1}='';
            end
            
        end
        
        if isnumeric(fget)==1
            
            if size(fget,1)==new_col-1
                fget(new_col,1)=0;
            end
        end
        
        model=setfield(model,fnames{i,1},fget);
        
    end
    
end


end

