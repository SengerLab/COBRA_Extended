function numatoms=count_atoms(formula, element)
%Count the number of atoms in a molecular formula
%
%[numatoms]=count_atoms(formula)
%
% Written by the Senger Research Group, 3/1/2013
% Last updated by the Senger Research Group, 6/9/2015 
%
%inputs:
%formula - the molecular formula of a molecule (e.g., H2O, Cu, C4H2O4)
%element - a specific element (e.g., C, N, Cu) (optional - leave empty to
%count all atoms in the formula)
%
%outputs:
%numatoms - the number of atoms in the molecular formula

%identify ions
s7=regexp(formula, '+|-');

%remove the ionic charge from the formula for atom counting
if size(s7,2)>0
    
    %accomodate single charge (simple) (e.g., Fe+3 is counted as Fe (1 atom), but
    %Fe3+ is counted as Fe3 (3 atoms). This is because of molecules like PO4-, where
    %the charge is -1 not -4. May return erroneous result based on what
    %nomenclature is used.
    formula=formula(1,1:s7(1,1)-1);

end

%locate all capital letters
if nargin == 1
    
    s1=regexp(formula,'[A-Z]'); %count all atoms
    
else
    
    s1=regexp(formula, element);    %count a specified element
    
end

%locate all lower case letters
s2=regexp(formula,'[a-z]');

%locate all single numbers 
s3=regexp(formula,'\D\d');

%locate all numbers 10-99 
s4=regexp(formula,'\D\d\d');

%locate all numbers 100-999
s5=regexp(formula,'\D\d\d\d');

%locate all numbers 1000-9999
s6=regexp(formula,'\D\d\d\d\d');



%initiate the molecular matrix, quick return if specific element not found
if size(s1,2)>0
    
    mol=zeros(size(s1,2),3);
    
else
    
    numatoms=0;
    
    return
    
end

%populate the capital letter ID (col 1), initialize lower-case ID (col2), initialize stoichiometry (col3):
for i=1:size(s1,2)
    
    mol(i,1)=s1(1,i);
    mol(i,2)=s1(1,i);
    
end

mol(:,3)=1;

%populate the lower-case letter ID (col 2):
if size(s2,2)>0
    
    for i=1:size(s1,2)
        
        for j=1:size(s2,2)
            
            if s2(1,j)-s1(1,i)==1
                
                mol(i,2)=s2(1,j);
                
            end
        end
    end
end

%adding in stoichiometry 1-9
if size(s3,2)>0
    
    for i=1:size(s3,2)
        
        stoich=str2num(formula(1,s3(1,i)+1));
        
        index=find(mol(:,2)==s3(1,i));
        
        mol(index,3)=stoich;
        
    end
end

%adding in stoichiometry 10-99 
if size(s4,2)>0
    
    for i=1:size(s4,2)
        
        stoich=str2num([formula(1,s4(1,i)+1) formula(1,s4(1,i)+2)]);
        
        index=find(mol(:,2)==s4(1,i));
        
        mol(index,3)=stoich;
        
    end
end

%adding in stoichiometry 100-999
if size(s5,2)>0

    for i=1:size(s5,2)
        
        stoich=str2num([formula(1,s5(1,i)+1) formula(1,s5(1,i)+2) formula(1,s5(1,i)+3)]);
        
        index=find(mol(:,2)==s5(1,i));
        
        mol(index,3)=stoich;
        
    end
end

%adding in stoichiometry 1000-9999
if size(s6,2)>0

    for i=1:size(s6,2)
        
        stoich=str2num([formula(1,s6(1,i)+1) formula(1,s6(1,i)+2) formula(1,s6(1,i)+3) formula(1,s6(1,i)+4)]);
        
        index=find(mol(:,2)==s6(1,i));
        
        mol(index,3)=stoich;
        
    end
end

%removing stoichiometry associated with lower-case letters
if nargin > 1
    
    for i=1:size(mol,1)
    
        if mol(i,2)-mol(i,1)==size(element,2)
            
            mol(i,3)=0;
            
        end
    end
end

%determing the total number of atoms
numatoms=sum(mol(:,3));

end

        
        
                
                
