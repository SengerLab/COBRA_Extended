function model=constrain_ratio(model_in,top_col,bottom_col,ratio,top_dir,bottom_dir)

%Written by Senger Research Group 7/13/2011, last updated 5/24/2012

%modeln=constrain_ratio(modeli,top_row,bottom_row,ratio) adds a constraint
%of the ratio between two reactions to be a certain value to a cobra-style
%modeln
%
%Inputs:
% model_in      cobra-style model
% top_row     column(s) of the stoichiometric matrix pertaining to 
%             reaction(s) on the top of the ratio (row vector)
% bottom_row  column(s) of the stoichiometric matrix pertaining to 
%             reaction(s) on the bottom of the ratio (row vector)
% ratio       value to constrain the ratio (sum of top rxns to sum of 
%             bottom rxns) to
% top_dir     row vector/value containing directionality info for top rxns (1
%             for forward, -1 for reverse) (default 1 for all)
% bottom_dir  row vector/value containing directionality info for bottom rxns
%             (1 for forward, -1 for reverse) (default 1 for all)
%
%Outputs:
% model       final cobra-style model with added constraint




model=model_in;

%set defaults if necessary
if nargin<6
    bottom_dir=ones(1,size(bottom_col,2));
end
if nargin<5
    top_dir=ones(1,size(top_col,2));
end

%add new constraint row to stoichiometric matrix
new_row=size(model.S,1)+1;
model.S(new_row,:)=zeros(1,size(model.S,2));

%modify the row to constrain as needed
model.S(new_row,top_col)=-1.*top_dir;
model.S(new_row,bottom_col)=ratio.*bottom_dir;

%add information to the metabolite fields
model.mets{new_row}=['ratio_' model.rxns{top_col} ':' model.rxns{bottom_col} '[c]'];
model.metNames{new_row}=['ratio_' model.rxns{top_col} ':' model.rxns{bottom_col}];

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