function [blockedRxns,  biomass,minFlux]...
    = FMDL(model,targetMet,P,minGR,minPR)

biomassRxnID=find(model.c);
[minFlux,Brange,Trange,biomass,TMPR,MB,blockedRxns,extype,targetRID,model2]=...
    constraintSearch(model,targetMet,minGR,minPR,biomassRxnID,P);

save('FMDL.mat');

end

