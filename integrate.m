function [y,opt4biomass,opt4target,opt4f,opt5biomass,blockedRxns] =...
    integrate(model2,model3,targetRID,biomassLB,biomassUB,targetLB,targetUB,n,biomassRxnID,lp)

model4=model3;
lp.lb(biomassRxnID)=biomassLB;
lp.ub(biomassRxnID)=biomassUB;
lp.lb(targetRID)=targetLB;
lp.ub(targetRID)=targetUB;
options = optimoptions('linprog','Display','none','Algorithm','dual-simplex');
[opt4.x, opt4.f, opt4.stat,opt4.output] = cplexlp(lp.f, lp.A, lp.b, lp.Aeq, lp.beq, lp.lb, lp.ub,options);
blockedRxns=[];
if opt4.stat~=1
    y=0;
    opt4biomass=-99999;
    opt4target=-99999;
    opt4f=-99999;
    opt5biomass=-99999;
    rxnList={};
    return
end
opt4biomass=opt4.x(biomassRxnID);
opt4target=opt4.x(targetRID);
opt4f=opt4.f;
usedRxns=find(abs(opt4.x)>=0.0000001);

for i=1:n
    if ismember(i,  usedRxns)==0
        blockedRxns=union(blockedRxns,i);
    end
end
blockedRxns=(blockedRxns)';

rxnList=model4.rxns(blockedRxns);
model5=changeRxnBounds(model2,rxnList,0,'b');
opt5=optimizeCbModel(model5);
if opt5.stat==1
    [minFlux5,maxFlux5]=fluxVariability(model5,100,'max',model5.rxns(targetRID));
end

switch opt5.stat
    case 1
        y=minFlux5;
        opt5biomass=opt5.x(biomassRxnID);
    otherwise
        y=0;
        opt5biomass=-999999;
end
return
end

