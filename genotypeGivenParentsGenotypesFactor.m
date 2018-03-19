function genotypeFactor = genotypeGivenParentsGenotypesFactor(numAlleles, genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo)
% This function computes a factor representing the CPD for the genotype of
% a child given the parents' genotypes.

% THE VARIABLE TO THE LEFT OF THE CONDITIONING BAR MUST BE THE FIRST
% VARIABLE IN THE .var FIELD FOR GRADING PURPOSES

% When writing this function, make sure to consider all possible genotypes 
% from both parents and all possible genotypes for the child.

% Input:
%   numAlleles: int that is the number of alleles
%   genotypeVarChild: Variable number corresponding to the variable for the
%   child's genotype (goes in the .var part of the factor)
%   genotypeVarParentOne: Variable number corresponding to the variable for
%   the first parent's genotype (goes in the .var part of the factor)
%   genotypeVarParentTwo: Variable number corresponding to the variable for
%   the second parent's genotype (goes in the .var part of the factor)
%
% Output:
%   genotypeFactor: Factor in which val is probability of the child having 
%   each genotype (note that this is the FULL CPD with no evidence 
%   observed)

% The number of genotypes is (number of alleles choose 2) + number of 
% alleles -- need to add number of alleles at the end to account for homozygotes

genotypeFactor = struct('var', [], 'card', [], 'val', []);

% Each allele has an ID.  Each genotype also has an ID.  We need allele and
% genotype IDs so that we know what genotype and alleles correspond to each
% probability in the .val part of the factor.  For example, the first entry
% in .val corresponds to the probability of having the genotype with
% genotype ID 1, which consists of having two copies of the allele with
% allele ID 1, given that both parents also have the genotype with genotype
% ID 1.  There is a mapping from a pair of allele IDs to genotype IDs and 
% from genotype IDs to a pair of allele IDs below; we compute this mapping 
% using generateAlleleGenotypeMappers(numAlleles). (A genotype consists of 
% 2 alleles.)

[allelesToGenotypes, genotypesToAlleles] = generateAlleleGenotypeMappers(numAlleles);

% One or both of these matrices might be useful.
%
%   1.  allelesToGenotypes: n x n matrix that maps pairs of allele IDs to 
%   genotype IDs, where n is the number of alleles -- if 
%   allelesToGenotypes(i, j) = k, then the genotype with ID k comprises of 
%   the alleles with IDs i and j
%
%   2.  genotypesToAlleles: m x 2 matrix of allele IDs, where m is the 
%   number of genotypes -- if genotypesToAlleles(k, :) = [i, j], then the 
%   genotype with ID k is comprised of the allele with ID i and the allele 
%   with ID j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Fill in genotypeFactor.var.  This should be a 1-D row vector.
genotypeFactor.var = [genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo];

% Fill in genotypeFactor.card.  This should be a 1-D row vector.
num = length(genotypeFactor.var);
indiceslength = length(genotypesToAlleles(:,1));
genotypeFactor.card = linspace(indiceslength,indiceslength,num);

genotypeFactor.val = zeros(1, prod(genotypeFactor.card));
% Replace the zeros in genotypeFactor.val with the correct values. get all
% assignment of those 27 pair genetype
assign_ = IndexToAssignment(1:prod(genotypeFactor.card),genotypeFactor.card);
% actully total difference of genetype are three 1=FF,2=Ff,3=ff
%thinking way: becasue we only have three differnt type and the child gene
%type also only have three, if we list all different combination of gene 
%type of parent1,2, the child gene type should be in the combined list, for
%each pair parent, the gene type only have four possible combination such
%as Ff mix with Ff they are 
% FF, Ff, Ff, ff so correpondent gene type index is 1, 2, 2, 3 because we
% are going to loop the whole assignement, so the child gene type must be
% in side of those four combination, and each type has an ID k, which we
% can count each child type probability in all possible parent
% mix.condition: sum(child type) = sum(parent type)
for i = 1:length(assign_)
     GeneType_child = assign_(i,1);
     GeneType_par1  = assign_(i,2);
     GeneType_par2  = assign_(i,3);
     Type_child = genotypesToAlleles(GeneType_child,:);
     Type_par1 = genotypesToAlleles(GeneType_par1,:);
     Type_par2 = genotypesToAlleles(GeneType_par2,:);
     combination = [];
     length_tmp_allpar = length(Type_par1);
     for j = 1:length_tmp_allpar
         for k = 1:length_tmp_allpar
             combination = [combination;Type_par1(j) Type_par2(k)];
         end
     end
     num_count = 0;
     length_tmp_combination = length(combination);
     for j = 1:length_tmp_combination
         combin_row = combination(j,:);
         if sort(Type_child) == sort(combin_row)
             num_count = num_count+1;
         end
     end
     %genotypeFactor=SetValueOfAssignment(genotypeFactor,assign_(i,:),num_count/length_tmp_combination);
     genotypeFactor.val(i)= num_count/length_tmp_combination;
end
 i= 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%