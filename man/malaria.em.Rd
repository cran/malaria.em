\name{malaria.em}
\alias{malaria.em}
\alias{geno}
\alias{Ginv}
\alias{findCandidatesForAllGenotypes}
\alias{findHaplotypeCandidateOneGenotype}
\alias{funOfLambda}
\alias{genAllHaploSetList}
\alias{genAllHaploSetList.C}
\alias{genAllHaploSetList.C.alldata}
\alias{genAllHaploSetList.c}
\alias{genAllHaploSetList.c.i}
\alias{genAllHaploSetListForOneGenotype}
\alias{getAllCoefficientQij.C}
\alias{getAllHaplotypes}
\alias{getAllQij.C}
\alias{getAllScore.C.em}
\alias{getCoefficientQij.C}
\alias{getDesignforAllGenotypes.additive.allhaplotype}
\alias{getDesignforOneGenotype.additive}
\alias{getDesignforOneGenotype.additive.allhaplotype}
\alias{getHaplotypesWithLoci}
\alias{getLocitypeFromOneGenotype}
\alias{getNextHaploSet}
\alias{getNextIndex}
\alias{getOneHaplotype}
\alias{getOneNumset.C}
\alias{getOneQij.C} 
\alias{getOneScore.C.em}
\alias{getOneScore.haplotype.em}
\alias{getStd.em}
\alias{getUniqueGenotype}
\alias{getUniqueGenotypeCount}
\alias{getVector.C}
\alias{getallConfig.C}
\alias{getallWeight.C.em}
\alias{getoneHaploSet.vector.C}
\alias{getoneL.C.em}
\alias{getoneQij.vector.C}
\alias{getoneWeight.C.em}
\alias{haploSet2Genotype}
\alias{isLastSet}
\alias{logL}
\alias{newton.glm}
\alias{num.unique.allele}
\alias{num.unique.allele.locus}
\alias{print.malaria.em}
\alias{sortGenotypes}
\alias{toGenotypes}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM Estimation of Malaria Haplotype Probabilities and Number of Infections 

from Multiply Infected Human Blood Samples
}
\description{
Compute the maximum likelihood estimates for malaria haplotype probabilities and number of infections based on the malaria genotype data using Expectation-Maximization approach. }
\usage{
malaria.em(geno, sizes=c(2), locus.label=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{
matrix of genotypes. Each column represents the alleles observed in each of the locus. The observed alleles in each locus are separated by a space, and the order of columns corresponds to the order of loci on a chromosome. If there are K loci, then ncol(geno) = K. Rows represent the genotypes for each subject.
}
  \item{sizes}{
An integer or a vector of possible number of parasite strains within individuals in observed data. For example, sizes=c(1:6) means the possible number of parasite strains is ranging from 1 to 6. If the length of sizes is greater than 1, the estimation will assume zero truncated Poisson distribution on the number of parasite strains.   }
  \item{locus.label}{
vector of labels for loci.
}
}
\details{
This program is an extension of haplo.em() function in haplo.stats package, which computes the maximum likelihood estimates (MLEs) of haplotype probabilities in diploid population. Since in malaria setting, the number of the malaria parasite strains is unknown due to the possible multiple infections, we impose a probability distribution on this number. Parameter sizes determines the underlying probability distribution. Length of sizes more than one will lead to the estimation assuming a truncated Poisson distribution. In this case, the estimation of the Poisson rate and number of infections will be reported. Otherwise, if sizes is defined as an integer- C for example, the estimation will assume the number of strains is fixed at C. Although this assumption is not applicable for malaria setting, it is very useful for any genetic data, in which the number of chromosomes is same across all samples.  By default, sizes=2, which reduces to haplo.em() and give the MLEs of the haplotype probabilities in diploid population.
}

\value{
\item{haplo.prob.tab}{matrix of unique haplotypes, MLEs of estimated haplotype probabilities, and their standard errors.}
\item{haplotype}{matrix of unique haplotypes. Each row represents a unique haplotype, and the number of columns is the number of loci.}
\item{haplo.prob}{vector of MLEs of haplotype probabilities. The ith element of hap.prob corresponds to the ith row of haplotype.}
\item{haplo.prob.std}{standard error of the estimated haplotype frequencies. }
\item{lambda}{estimated Poisson parameter.}
\item{NumofInfection}{estimated number of infections.}
\item{haplo.sets}{ List of all possible haplotype combinations and their posterior probability per subject. The first column named ids is a vector for row index of subjects after expanding to all possible haplotype combinations for each subject. If ids=i, then i is the ith row of geno. If the ith subject has n possible haplotype combinations that correspond to their marker genotype, then i is repeated n times. The value in the second column is the row numbers of the unique haplotypes in the returned haplotype matrix.}
\item{n.haplo.set}{vector of maximum number of haplotype combinations per subject that are consistent with their marker data in the matrix geno. The length of n.haplo.set = nrow(geno).}
\item{pred.haplo.set}{Predicted haplotype combination that is consistent with their marker data for each subject. The values in pred.haplo.set are the row numbers of the unique haplotypes in the returned haplotype matrix. }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Li, X., Foulkes, A.S., Yucel, R. and Rich, S.M. (2007) An expectation maximization approach to estimate malaria haplotype frequencies in multiply infected children, Statistical Applications in Genetics and Molecular Biology, Vol. 6 : Iss. 1, Article 33.
}
\author{
Xiaohong Li
}
\examples{
data(geno)

sizes<-c(2)
ret1<-malaria.em(geno,sizes,locus.label=c("DQB","DRB") )

sizes<-c(1:3)
ret2<-malaria.em(geno,sizes,locus.label=c("DQB","DRB") )

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


