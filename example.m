function [outputArg1,outputArg2] = example(inputArg1,inputArg2)

load('e_coli_core.mat');
P=10;
minGR=0.00001;
minPR=0.00001;
[blockedRxns,  biomass,minFlux]=...
    FMDL(e_coli_core,{'pyr_e'},P,minGR,minPR);
save('example1.mat');


%load('iML1515.mat');
%P=5;
%minGR=0.00001;
%minPR=0.00001;
%[blockedRxns,  biomass,minFlux]=...
    %FMDL(iML1515,{'pyr_e'},P,minGR,minPR);
%save('example2.mat');


scatter(biomass,minFlux)
xlim([0 1])
ylim([0 20])
set(gca,'FontSize',20)

end

