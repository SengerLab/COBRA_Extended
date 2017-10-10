function r=import_constraints(r,data)

%written by the Senger Research Group, 4/24/2012, last updated 9/26/2012

%this simple program allows me to copy constraints and the objective function column from excel, import, and
%immediately incorporate into the model.

%copy paste 2 columns: col1 - LB; col2 - ub
%or copy paste 3 columns: col1 - LB, col2 -ub, col3 - obj

r.lb=data(:,1);
r.ub=data(:,2);

if size(data,2)==3
    r.c(:,1)=data(:,3);
end

fprintf(['Constraints and objectives are up to date. Save manually.' '\n']);