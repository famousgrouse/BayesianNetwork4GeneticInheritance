% You can create a network and convert it to a format that can be viewed in
% Samiam by running this script.

pedigree = struct('parents', [0,0;1,3;0,0;1,3;2,6;0,0;2,6;4,9;0,0]);
pedigree.names = {'Ira','James','Robin','Eva','Jason','Rene','Benjamin','Sandra','Aaron'};
alleleFreqs = [0.1; 0.2;0.7];
alleleList = {'F', 'f','n'};
alphaList = [0.1; 0.1; 0.1;0.2;0.3;0.4;0.5;0.6;0.8];
phenotypeList = {'CysticFibrosis', 'NoCysticFibrosis'};
positions = [520, 600, 520, 500; 650, 400, 650, 300; 390, 600, 390, 500; 260, 400, 260, 300; 780, 200, 780, 100; 1040, 400, 1040, 300; 910, 200, 910, 100; 130, 200, 130, 100; 0, 400, 0, 300];

% This will construct a Bayesian network and convert it into a file that 
% can be viewed in SamIam.

factorList = constructGeneticNetwork(pedigree, alleleFreqs, alphaList);
sendToSamiam(pedigree, factorList, alleleList, phenotypeList, positions, 'cysticFibrosisBayesNet');