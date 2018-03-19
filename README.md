# Bayes Nets for Genetic Inheritance

Using bayesian network to model the genetic Suppose a genetic pedigree. From this
pedigree, we can construct a Bayesian network to model the process of genetic inheritance. The network consists three person pedigree with six nodes, which means each person will have two nodes. For each person, there is a factor (factor type 1) for P(person’s phenotype | person’s genotype). Each person’s genotype is determined by that person’s parents’ genotypes. Each person’s genotype is determined by that person’s parents’ genotypes, so, for each person, there is a factor (factor type 3) for P(person’s genotype | genotype of person’s first
parent, genotype of person’s second parent). The genotype factor for a person whose parents are not specified is just the prior P(person’s genotype) (factor type 2), and its values are based on the allele frequencies in the population.
