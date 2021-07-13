function [production,Brange,Trange,biomass,TMPR,MB,blockedRxns,extype,targetRID,model2,redundant]=...
    constraintSearch(model,targetMet,minGR,minPR,biomassRxnID,P)


target=findMetIDs(model,targetMet);
m=size(model.mets,1);
n=size(model.rxns,1);
if isempty(find(strcmp(model.rxns,strcat('EX_',targetMet))))==0
    targetRID=find(strcmp(model.rxns,strcat('EX_',targetMet)));
    model2=model;
    extype=1;
elseif isempty(find(strcmp(model.rxns,strcat('DM_',targetMet))))==0
    targetRID=find(strcmp(model.rxns,strcat('DM_',targetMet)));
    model2=model;
    extype=2;
else
    [model2,rxnIDexists]=addReaction(model,'Transport',targetMet,[-1]);
    m=size(model2.mets,1);
    n=size(model2.rxns,1);
    model2.S(target,n)=-1;
    model2.ub(n)=999999;
    model2.lb(n)=0;
    model2.rev(n)=0;
    targetRID=n;
    extype=3;
end

opt2=optimizeCbModel(model2);
MB=opt2.f;
model3=model2;
model3.c(biomassRxnID)=0;
model3.c(targetRID)=1;
model3.lb(biomassRxnID)=minGR;
opt2=optimizeCbModel(model3);
TMPR=opt2.f;

model3=model2;
model3.lb(biomassRxnID)=minGR;
lp.f=[zeros(n,1); ones(n,1)];
lp.A=[-eye(n) -eye(n);eye(n) -eye(n)];
lp.Aeq=[model3.S zeros(m,n)];
lp.beq=zeros(m,1);
lp.lb=[model3.lb; zeros(n,1)];

lp.ub=[model3.ub; 999999*ones(n,1)];
tableBest=0;
biomass(1,1)=0;
Brange=0;Trange=0;
production=0;
blockedRxns={};
rList={};
x=1;redundant(1,1)=0;
mmm=0.00001;
for i=1:P
    biomassLB=mmm+((MB-mmm)/P)*(i-1);
    biomassUB=mmm+((MB-mmm)/P)*i;
    for j=1:P
        targetLB=mmm+((TMPR-mmm)/P)*(j-1);
        targetUB=mmm+((TMPR-mmm)/P)*j;

        lp.b=[zeros(2*n,1)];
        [table(i,j),opt4biomass,opt4target,opt4f,opt5biomass,rList]=...
            integrate(model2,model3,targetRID,biomassLB,biomassUB,targetLB,targetUB,n,biomassRxnID,lp);

        if (opt5biomass >= minGR)&&(table(i,j)>=minPR)
            flag=0;
            for w=1:x-1
                if isequal(rList,blockedRxns{w,1})==1
                    flag=1;
                end
            end
            if flag==0
                biomass(x,1)=opt5biomass;
                production(x,1)=table(i,j);
                Brange(x,1)=i;Trange(x,1)=j;
                blockedRxns{x,1}=rList;
                x=x+1
            end
        end
    end
end
end

