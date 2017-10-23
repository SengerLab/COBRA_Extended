function model = InVitroModelSetup(model)

%Written by Ryan Senger, 10/22/2017

%This function creates a model for in vitro metabolic engineering by
%1. disabling all membrane transport and exchange reactions
%2. creating new exchange reactions for all cytosolic compounds

%note: this function uses the findExcRxns.m function of the COBRA toolbox,
%the findTransRxns2.m function, and the add_exchange function from the Senger Lab


%find transport reactions
isTrans = findTransRxns2(model);

%find exchange reactions
selExc = findExcRxns(model);

%disable all membrane transport and exchange reactions
for i = 1:size(model.rxns,1)
    
    if isTrans(i,1) == 1    %it is a transport reaction
        
        model.lb(i,1) = 0;
        model.ub(i,1) = 0;
        
    end
    
    if selExc(i,1) == 1     %it is an exchange reaction
        
        model.lb(i,1) = 0;
        model.ub(i,1) = 0;
        
    end
    
end

%add exchange reactions for all cytosolic compounds "[c]"
for i = 1:size(model.mets,1)
    
    fprintf(['Working on reaction ' num2str(i) ' of ' num2str(size(model.mets,1)) '\n']);

    if isempty(findstr('[c]',model.mets{i,1})) == 0     %finds the [c] compartment
        
        model=add_exchange(model,i,-1000,1000); %adds the exchange reaction
        
    end
    
end

fprintf(['done!' '\n']);

