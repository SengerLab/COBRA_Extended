function m = constrainOutFlux(m0,node,target,ratio)

%written by the Senger Research Group, 9/21/2012

% m0 - the model
% node - the met id where node is considered
% target - the rxn id in which the ratio is targeted
% ratio - value between 0 and 1

rxns = find(m0.S(node,:) < 0);
m = m0;

if length(rxns) == 1
    error('Only one outgoing rxns\n');
end

if isempty(find(rxns == target))
    error('wrong target\n');
end

n = length(m0.mets) + 1;
for i = 1:length(rxns)
    if rxns(i) == target
        m.S(n,rxns(i)) = 1-ratio;
    else
        m.S(n,rxns(i)) = -ratio;
    end
end

m.mets(n) = {['Ratio_' num2str(node)]};
m.metNames(n) = {['Ratio_' num2str(node)]};
