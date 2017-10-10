function model=delete_reaction(model,reaction_line)

%written by the Senger Research Group, 9/23/2012

%this function deletes an exchange reaction (or any other reaction) from a
%COBRA style model.

%Inputs
%model: the COBRA style model
%reaction_line: the ID of the reaction to be deleted (this is the row of
%the reaction in model.rxns)

%get all fields for updating
fnames=fieldnames(model);

%delete reaction information from all fields
if size(model.S,1)~=size(model.S,2)
    
    for i=1:size(fnames,1)
        
        %extract the field as a variable
        fget=getfield(model,fnames{i,1});
        
        %check the size and make modifications if necessary
        if size(fget,1)==size(model.S,2)
            
            %delete the information
            fget(reaction_line,:)=[];
            
            %replace the field
            model=setfield(model,fnames{i,1},fget);
            
        end
        
    end
    
    %delete the column in the stoichoimetric matrix
    model.S(:,reaction_line)=[];
    
else
    
    fprintf(['The model has the same number of metabolites and reactions. Deletion must be done manually.' '\n']);
    return
    
end
