gc()
### get all haplotypes from a list of unique genotypes
### genotypes is a list of genotype strings
getAllHaplotypes<- function(genotypes, allPrefix ){
  allHaplotypes <- list()
  ng <- length(genotypes )
  for( i in 1:ng ){
    g <- genotypes [[i]]
    if( is.na(g) ){
      next
    }
    
    loci <- getLocitypeFromOneGenotype(g)
    hs <- getHaplotypesWithLoci(loci, allPrefix)
    allHaplotypes <- c(allHaplotypes, hs)
  }
  allHaplotypes <- unique(allHaplotypes)
  allHaplotypes
}

### given a genotype that is a string,
### return a list of locitypes
getLocitypeFromOneGenotype <- function(g){
  gs <- strsplit(g, split=",")[[1]]
  genoSize <- length(gs)
  
  lociSet <- list()
  for( i in 1: genoSize ) {
    ss <- gs[i]
    slist <- strsplit(ss, split=" ")[[1]]
    hu <- sort(unique(slist))
    lociSet <- c(lociSet, list(hu))
  }
  lociSet
}

### generate all haplotypes based on all locitypes
### suppose locus i have ni values, and there are k differnt loci,
### then the number haplotypes is n1.n2.....nk
getHaplotypesWithLoci <- function(allLoci, allPrefix){
  nloci <- length(allLoci)
  index <- rep(1, times=nloci)
  
  total <- 1
  for( i in 1:nloci ){
    total <- total * length(allLoci[[i]])
  }
  
  allHaplotypes <- list()
  for( i in 1:total ){
    haplo <- getOneHaplotype(allLoci, index)
    if( length(allPrefix)>0 ){
      hp <- haplo[1:nloci-1]
      hprefix<- apply(as.matrix(hp), 2, paste, collapse="")
      if(any(i <- grep(hprefix,allPrefix))){
        allHaplotypes <- c(allHaplotypes, list(haplo))
      }
    }
    else{
      allHaplotypes <- c(allHaplotypes, list(haplo))
    }
    
    index <- getNextIndex(allLoci, index)
  }
  
  allHaplotypes
}


### get one haplotype based on index
getOneHaplotype <- function(allLoci, index){
  haplo <- c()
  nl <- length(index)
  for( i in 1:nl ){
    haplo <- c(haplo, allLoci[[i]][index[i]])
  }
  haplo
}

### get the next index
getNextIndex <- function(allLoci, index){
  ret <- index
  nloci<-length(allLoci)
  for( i in nloci:1 ){
    ni <- length(allLoci[[i]])
    if( index[i]<ni ){
      ret[i] <- index[i] + 1
      break
    }
    else{
      ret[i] <- 1
    }
  }
  ret
}

##########################################################
##sort genoData by genotypes
## genoData is a matrix

sortGenotypes<- function(genoData){
  allGenotypes <- c()
  ng <- dim(genoData)[1]
  nloci <- dim(genoData)[2]
  for( i in 1:ng ) {
    if( is.na(genoData[i, 1])) {
      next
    }
    g <- "";
    
    for( j in 1:nloci ){
      if( is.na(genoData[i, j])) {
        next
      }
      
      t <- strsplit(genoData[i, j], split=" ")
      ss<- sort(unique(t[[1]]))
      ss <- apply(as.matrix(ss), 2, paste, collapse=" ")
      if( j==1 ){
        g <- ss
      }
      else{
        g <- paste(g, ss, sep=",")
      }
    }
    allGenotypes <- c(allGenotypes, g)
  }
  allGenotypes <- sort(allGenotypes)
  allGenotypes 
}

## find all haplotypes that can be part of the given genotype
## return a vector of indexes of such haplotypes
## genotype is a vector of strings
findHaplotypeCandidateOneGenotype<- function(allHaplotypes, genotypeString){
  candidate <- c()
  nH <- length(allHaplotypes)
  genotype <- strsplit(genotypeString, split=",")[[1]]
  nLoci <- length(genotype)
  for( i in 1:nH ){
    h <- allHaplotypes[[i]]
    isCandidate <- TRUE
    for( j in 1:nLoci ){
      itmp<-grep(h[j], genotype[j])
      if( any(itmp) ){
      }
      else{
        isCandidate <- FALSE
        break
      }
    }
    
    if( isCandidate ){
      candidate <- c(candidate, i)
    }
  }
  candidate	
}

findCandidatesForAllGenotypes<- function(allHaplotypes, genotypes){
  nG <- length(genotypes)
  candidates <- vector("list", nG)
  for( gi in 1:nG ){
    cgi <- findHaplotypeCandidateOneGenotype(allHaplotypes, genotypes[[gi]])
    if( length(cgi)>0 ){
      candidates[[gi]] <- cgi
    }
  }
  candidates 
}

## generate all haplo set of a given size for one genotypes
## setSize is the size of the set
genAllHaploSetListForOneGenotype<- function( allHaplotypes,
                                             setSize,
                                             genotypeString,
                                             NLoci,
                                             candidate ){
  nH <- length(allHaplotypes)
  n.candidates <- length(candidate)
  haploSetList <- list()

  if( n.candidates>0 ){
    h.candidate <- rep(1, times=setSize)
    while( TRUE ){
      h <- h.candidate
      for( hi in 1:setSize){
        h[hi] <- candidate[h.candidate[hi]]
      }

      g <- haploSet2Genotype (h, allHaplotypes, NLoci)
      if( g==genotypeString ){
        haploSetList <- c( haploSetList, list(h))
      }

      if( isLastSet(n.candidates, h.candidate) ){
        break
      }

      h.candidate <- getNextHaploSet (n.candidates, h.candidate)
    }
  }
  
  haploSetList
}


## Given a set of genotypes and all haplotypes, find
## all haplotype set for all geneotype
genAllHaploSetList <- function( allHaplotypes, setSize, genotypes, NLoci, candidates ){
  nG <- length(genotypes)
  
  ##generate nG list, each of which store haploset for a genotype
  haploSetList <- vector("list", nG)
  for( gi in 1:nG ){
    haploSetList[[gi]] <- list()
  }
  for( gi in 1:nG ){
    ##note: candidates[[gi]] can be null
    if( !is.null(candidates[[gi]]) ){
      lengi<-length(candidates[[gi]])
      ##cat("\n for genotype " , gi, "; candidate size = ", lengi)
      if( lengi>0 ){
        haploSetList[[gi]] <- genAllHaploSetListForOneGenotype( allHaplotypes,
                                                               setSize,
                                                               genotypes[[gi]],
                                                               NLoci,
                                                               candidates[[gi]] )
      }
      ##cat("\n", gi, ",", lengi, ",", length(haploSetList[[gi]]), "  ;  ")
    }
  }
  ##cat ("\n")
  haploSetList
}


## to compare if the haplo set s is the last set
## the last set = (nH, nH, ... nH)
isLastSet <- function( nH, s ){
  n <- length(s)
  ret <- TRUE
  for( i in 1:n ) {
    if( s[i]!=nH ){
      ret <- FALSE
      break
    }
  }
  ret
}

## get the next set based on the current set
## nH is the largest index of haplotype
## e.g.: {2,2,1} --next--> {2,2,2} assuming nH=2
getNextHaploSet <- function( nH, cur){
  nset <- cur
  
  n <- length(cur)
  for( i in 1:n ){
    if( nset[i]<nH ){
      nset[i] <- nset[i]+1
      for( j in 1:(i-1) ){
        nset[j] <- nset[i]
      }
      break
    }
  }
  
  nset
}

haploSet2Genotype <- function( haploIndexSet, allHaplotypes, NLoci){
  ns <- length(haploIndexSet)
  g <- rep("", times=NLoci)
  for( i in 1:NLoci ){
    lociSet <- c()
    for( j in 1:ns ){
      index <- haploIndexSet[j]
      h <- allHaplotypes[[index]]
      lociSet <- c(lociSet, h[i])
    }
    ul <- sort(unique(lociSet))
    ss <- apply(as.matrix(ul), 2, paste, collapse=" ")
    g[i] <- ss
  }
  g
  g <- apply(as.matrix(g), 2, paste, collapse=",")
  g
}


## generate all haploset with sizes in the given vector
## sizes: the given size vector
## nH is the number of haplotypes
genAllHaploSetList.C<-function(allHaplotypes, genotypes, sizes, NLoci){
  candidates <- findCandidatesForAllGenotypes (allHaplotypes, genotypes)
  cat("finding possible haplotype combinations... \n")
  
  myHaploSetList.C<-list()
  sl <- length(sizes)
  
  for (i in 1:sl) {
    setSize <- sizes[i]
    cat("haploset size = " , setSize , "\n")
    newSetList<-genAllHaploSetList(allHaplotypes,
                                   setSize,
                                   genotypes,
                                   NLoci,
                                   candidates) 
    myHaploSetList.C<-c(myHaploSetList.C,list(newSetList))
  }

  myHaploSetList.C
}

toGenotypes <- function(genoData){
	allGenotypes <- c()
	ng <- dim(genoData)[1]
	nloci <- dim(genoData)[2]
	for( i in 1:ng ) {
		if( is.na(genoData[i, 1])) {
			next
		}
		g <- "";

		for( j in 1:nloci ){
			if( is.na(genoData[i, j])) {
				next
			}

			t <- strsplit(genoData[i, j], split=" ")
			ss<- sort(unique(t[[1]]))
			ss <- apply(as.matrix(ss), 2, paste, collapse=" ")
			if( j==1 ){
				g <- ss
			}
			else{
				g <- paste(g, ss, sep=",")
			}
		}
		allGenotypes <- c(allGenotypes, g)
	}
	allGenotypes 
}
getDesignforOneGenotype.additive<-function(HaploSet, nhaplo){
	nset<-length(HaploSet)
	if (nset==0){
		xmat<-mat.or.vec(nset, nhaplo)
		xmat
	}
	else{
		xmat<-mat.or.vec(nset, nhaplo)
	
		nsize<-length(HaploSet[[1]])

		for (i in 1:nset){
			for(j in 1:nsize){
				hIndex<-HaploSet[[i]][j]
				xmat[i,hIndex]<-xmat[i,hIndex]+1
			}
		}
		xmat	
	}
}


getDesignforOneGenotype.additive.allhaplotype<-function (allHaploSet.C.alldata, i,nhaplo){
	nc<-length(allHaploSet.C.alldata)
	#get the design matrix for the i individual with the first haploset size
	xmat<-getDesignforOneGenotype.additive(allHaploSet.C.alldata[[1]][[i]], nhaplo)
	if (nc<2) {
		xmat<-xmat
		}
	else {	
		for(c in 2:nc){
			xmat.new<-getDesignforOneGenotype.additive(allHaploSet.C.alldata[[c]][[i]], nhaplo)
			xmat<-rbind(xmat,xmat.new)
			}
		}
		xmat
	}	
	
getDesignforAllGenotypes.additive.allhaplotype<-function(allHaploSet.C.alldata, allHaplotypes){
	nhaplo<-length(allHaplotypes)
	ng<-length(allHaploSet.C.alldata[[1]])
	
	xmat<-getDesignforOneGenotype.additive.allhaplotype(allHaploSet.C.alldata, 1, nhaplo)

	for(i in 2:ng){
		xmat.new<- getDesignforOneGenotype.additive.allhaplotype(allHaploSet.C.alldata, i, nhaplo)
		xmat<-rbind(xmat,xmat.new)
		}
		xmat
	}
getVector.C<-function(xmat.additive){
	n<-dim(xmat.additive)[1]
	cvector<-rep(0,n)
	for(i in 1:n){
		cvector[i]<-sum(xmat.additive[i,])
	}
	cvector
}

#####################################################################
##### E step
# allHaploSet.c 
# 1st dimension: (c), set size
# 2nd dimension: (i), genotype
# 3rd dimension: (j), haplotype combination
getAllCoefficientQij.C <- function(allHaploSet.C)
{
  allQijList.C <- list()
  cLen <- length(allHaploSet.C)

  for( c in 1:cLen )
  {
	facSize <- factorial(c)
	qiList<-list()
	iLen <- length(allHaploSet.C[[c]])
	for(i in 1:iLen)
	{	
		qicList <- c()
		hic <- allHaploSet.C[[c]][[i]]
 		jLen <- length(hic)
		if( jLen >0 )
		{
		  for ( j in 1: jLen )
		  {
	  		qicj <- getCoefficientQij.C(hic[[j]], facSize)
	  		qicList <- c(qicList, qicj)
		  }
		}
		else
		{
			#qicList <- c(qicList, 0)
		}

		qiList<-c(qiList,list(qicList))
	}

	allQijList.C <- c( allQijList.C, list(qiList))
  }
  allQijList.C
}



#oneHaploSet is a haplo set, like (5, 3, 3)
# oneHaploSet is a vector
# for (5,3,3), the returned value will be: facSize / (1! . 2!)
getCoefficientQij.C <- function(oneHaploSet, facSize)
{
 t <- table(oneHaploSet, deparse.level=0)
 tlen <- length(t)
 dd <- 1
 for( i in 1:tlen )
 {
	tmp <- factorial(as.numeric(t[i]))
	dd <- dd * tmp
 }
 qicj <- facSize/dd
 qicj
}



#oneHaploSet is a haplo set, like (5, 3, 3)
# oneHaploSet is a vector
getOneQij.C <- function(oneHaploSet, allEstiPs)
{
 len <- length(oneHaploSet)
 dd <- 1
 for( i in 1:len )
 {
	dd <- dd * allEstiPs[oneHaploSet[i]]
 }
 dd
}




getAllQij.C <- function(allQijCoefficients.C, 
				allHaploSet.C, 
				allEstiPs,
				EstiLambda,
				nObs,
				genoCountList)
{
	allQij.C<- allQijCoefficients.C
  	n<-sum(genoCountList)
  	n0<-nObs-n
   	cLen <- length(allHaploSet.C)
  	for( c in 1:cLen )
  	{
 		iLen <- length(allHaploSet.C[[c]])
	for( i in 1:iLen )
	{
		jLen <- length(allHaploSet.C[[c]][[i]])
		if( jLen>0 ){
			for ( j in 1: jLen  ){
		  		allQij.C[[c]][[i]][j]<- allQijCoefficients.C[[c]][[i]][j] * getOneQij.C(allHaploSet.C[[c]][[i]][[j]], allEstiPs) * funOfLambda(EstiLambda, c, n0)
		 	}
		}
	}
  }
allQij.C
}

# Find the corresponding haplo set with size c for the ith obs 
genAllHaploSetList.c.i<-function(allHaploSet.C, genotypes, genotypes.alldata, c, i){
	ng.unique<-length(genotypes)
	index<-0
	for (g in 1:ng.unique){
		if (genotypes.alldata[i]==genotypes[g]){
			index<-g
			break
		}
	}
	myHaploSetList.c.i<-allHaploSet.C[[c]][[index]]
	myHaploSetList.c.i
}

genAllHaploSetList.c<-function(allHaploSet.C, genotypes, genotypes.alldata, c){
	ni<-length(genotypes.alldata)
	myHaploSetList.c<-list()
	for (i in 1:ni){
		myHaploSetList.c.i<-genAllHaploSetList.c.i(allHaploSet.C, genotypes, genotypes.alldata, c, i)
		myHaploSetList.c<-c(myHaploSetList.c, list(myHaploSetList.c.i))	
	}
	myHaploSetList.c
}			
genAllHaploSetList.C.alldata<-function(allHaploSet.C, genotypes, genotypes.alldata){
	nc<-length(allHaploSet.C)
	myHaploSetList.C.alldata<-list()
	for (c in 1:nc){
		myHaploSetList.c<-genAllHaploSetList.c(allHaploSet.C, genotypes, genotypes.alldata, c)
		myHaploSetList.C.alldata<-c(myHaploSetList.C.alldata, list(myHaploSetList.c))
	}
	myHaploSetList.C.alldata
}

## Get the weight of each possible haplotype set for a given obs with genotype information
# return a vector of weight for each haplotype set 	
getoneWeight.C.em<-function(i, allHaploSet.C.alldata,allQij.C.alldata){
	nset<-getOneNumset.C(i, allHaploSet.C.alldata)	
	if (nset==0){
		w<-NULL
	}
	else{
		Qij.vector.C<-getoneQij.vector.C(allQij.C.alldata, i)
		w<-Qij.vector.C
		totalw<-sum(w)
		for(j in 1:nset){
			w[j]<-w[j]/totalw
		}
	}
	w
}


getallWeight.C.em<-function(allHaploSet.C.alldata, allQij.C.alldata){
#	nhaplo<-length(allHaplotypes)
	ng<-length(allHaploSet.C.alldata[[1]])
	w<-getoneWeight.C.em(1, allHaploSet.C.alldata, allQij.C.alldata)
	for (i in 2:ng){
		w<-c(w,getoneWeight.C.em(i, allHaploSet.C.alldata,allQij.C.alldata))
	}
	w
}

## Get the weight of each possible haplotype set for a given obs with genotype information
# return a vector of weight for each haplotype set 	
getoneL.C.em<-function(i, allHaploSet.C.alldata,allQij.C.alldata){
	#w<-getoneWeight.C.em(i, allHaploSet.C.alldata, allQij.C.alldata)
	nset<-getOneNumset.C(i, allHaploSet.C.alldata)	
	if (nset==0){
		logL<-NULL
	}else{
		Qij.vector.C<-getoneQij.vector.C(allQij.C.alldata, i)
		logL<-log(sum(Qij.vector.C))
	}
	logL
}


logL<-function(allHaploSet.C.alldata, allQij.C.alldata,EstiLambda,n0){
#	nhaplo<-length(allHaplotypes)
	ng<-length(allHaploSet.C.alldata[[1]])
	logL<--n0*EstiLambda

	for (i in 1:ng){
		logL<-logL+getoneL.C.em(i, allHaploSet.C.alldata,allQij.C.alldata)
	}
	logL
}


getOneNumset.C<-function(i, allHaploSet.C.alldata){
	nc<-length(allHaploSet.C.alldata)
	nset<-length(allHaploSet.C.alldata[[1]][[i]])
	if (nc<2){
		nset<-nset
	}else{
		for(c in 2:nc){
			newnset<-length(allHaploSet.C.alldata[[c]][[i]])
			nset<-nset+newnset		}
	}
	nset
}

getoneQij.vector.C<-function(allQij.C, i){
	nc<-length(allQij.C)
	Qij.vector.C<-allQij.C[[1]][[i]]
	if (nc<2){
		Qij.vector.C<-Qij.vector.C
		}else {
			for (c in 2:nc){
				Qij.vector.C<-c(Qij.vector.C, allQij.C[[c]][[i]]) 
			}
		}
	Qij.vector.C	
}

getallConfig.C<-function(allHaploSet.C.alldata,allQij.C.alldata){
	post.p<-NULL
	haplo.set<-NULL
	ids<-NULL
	n.haplo.set<-NULL
	pred.haplo.set<-NULL
	for(i in 1:length(allHaploSet.C.alldata[[1]])){
		nset<-getOneNumset.C(i, allHaploSet.C.alldata)	
		Qij.vector.C<-getoneQij.vector.C(allQij.C.alldata, i)
		w<-Qij.vector.C
		totalw<-sum(w)
		for(j in 1:nset){
			w[j]<-w[j]/totalw
		}
		o<-order(w,decreasing=TRUE)
		w<-w[o]
		post.p<-c(post.p,w)		
		haplo.set.i<-getoneHaploSet.vector.C(allHaploSet.C.alldata,i)
		haplo.set.i<-haplo.set.i[o]
		ids<-c(ids, rep(i, length(w)))
		haplo.set<-c(haplo.set, haplo.set.i)
		n.haplo.set[i]<-length(w)
		pred.haplo.set[i]<-haplo.set.i[1]
	}
	haplo.sets<-data.frame(ids, haplo.set, post.p)
	ret<-list(haplo.sets=haplo.sets, pred.haplo.set=pred.haplo.set, n.haplo.set=n.haplo.set)
	ret
}


getoneHaploSet.vector.C<-function(allHaploSet.C.alldata,i){
	nc<-length(allHaploSet.C.alldata)
	haplo.set.vector.C<-c()
	for(c in 1:nc){
		jLen <- length(allHaploSet.C.alldata[[c]][[i]])
		if( jLen<=0 ) {
			next
		}
		for ( j in 1: jLen  ) {
			hs <- allHaploSet.C.alldata[[c]][[i]][j]
			for(s in 1:length(hs[[1]])){
				if(s==1){
					hs.new<-hs[[1]][[s]]
				}else{
					hs.new<-paste(hs.new, hs[[1]][[s]],sep=" ")
				}
			}
			haplo.set.vector.C<-c(haplo.set.vector.C, hs.new)
		}
	}
	haplo.set.vector.C
}

funOfLambda<- function(EstiLambda, c,n0)
{
	#If n0>0,use unconditional poisson;
	#Otherwise, use conditional poisson


	r <- 1
	if( EstiLambda<0 )
	{
		r <- 1
	}
	else
	{
		if(n0>0){
			r <- exp(-EstiLambda)*(EstiLambda)^c/factorial(c)
			}
		else{
			r <- ((EstiLambda)^c/factorial(c))/(exp(EstiLambda)-1)
			}
	}
	r
}
##Newton's method to find the estimate of lambda in glm setting (is OK to use in runEM.new)
newton.glm<-function(w, cvector, EstiLambda)
{
	fLambda<-sum(w*cvector/EstiLambda-w*exp(EstiLambda)/(exp(EstiLambda)-1))
	flambda1<-sum(w*exp(EstiLambda)/(exp(EstiLambda)-1)^2-w*cvector/(EstiLambda^2))
	EstiLambda<<-EstiLambda-fLambda/flambda1
}

################################################################################################################
#get the score function for haplotype for one haplo set of one observation
# return a vector w[i]*d(logLi)/dP
getOneScore.haplotype.em<-function(w,i, xmat.additive, allEstiPs){
	nhaplo<-dim(xmat.additive)[2]
	si.haplotype<-rep(0, nhaplo-1)
	for (j in 1:nhaplo-1){
		si.haplotype[j]<-w[i]*(xmat.additive[i,j]/allEstiPs[j]-xmat.additive[i,nhaplo]/allEstiPs[nhaplo])
		}
	si.haplotype
}

	
## Get the empirical observed information matrix for the ith observation
# 	Reture a matrix 
	
	
getOneScore.C.em<-function(i, nhaplo,allHaploSet.C.alldata,allQij.C.alldata, allEstiPs,EstiLambda,n0){
	nset<-getOneNumset.C(i, allHaploSet.C.alldata)
	if (nset>0){
		w<-getoneWeight.C.em(i, allHaploSet.C.alldata, allQij.C.alldata)
		xmat.additive<-getDesignforOneGenotype.additive.allhaplotype(allHaploSet.C.alldata, i,nhaplo)
		cvector<-getVector.C(xmat.additive)		
		si.haplotype<- getOneScore.haplotype.em(w,1, xmat.additive, allEstiPs)
		#si.lambda<-t(w)%*%cvector
		if(n0>0){
			si.lambda<-t(w)%*%(cvector/EstiLambda-1)
		}else{
			si.lambda<-t(w)%*%(cvector/EstiLambda-exp(EstiLambda)/(exp(EstiLambda)-1))
		}
		if (nset>1){
		for(i in 2:nset){
			si.haplotype<-si.haplotype + getOneScore.haplotype.em(w,i, xmat.additive, allEstiPs)
			}
		}
	}else{
		si.haplotype<-xmat.additive[1,1:nhaplo-1]*0
		si.lambda<-0
		
	}
	oneScore<-c(si.haplotype,si.lambda)%*%t(c(si.haplotype,si.lambda))
	oneScore
}

## Get the empirical observed information matrix for all observations

getAllScore.C.em<-function(nhaplo,allHaploSet.C.alldata,allQij.C.alldata, allEstiPs,EstiLambda,n0){
	ng<-length(allQij.C.alldata[[1]])
	score<-getOneScore.C.em(1, nhaplo,allHaploSet.C.alldata,allQij.C.alldata, allEstiPs,EstiLambda,n0)
	for(i in 2:ng){
		score<-score + getOneScore.C.em(i, nhaplo,allHaploSet.C.alldata,allQij.C.alldata,allEstiPs,EstiLambda,n0)
	}
	score
}

getStd.em<-function(allScore){
	nb<-dim(allScore)[1]
	invScore<-Ginv(allScore)$Ginv
	
	lambda.std<-sqrt(invScore[nb,nb])
	
	invScore.haplo<-as.matrix(invScore[-nb,-nb])

	haplo.std<-rep(0,nb)
	for(i in 1:(nb-1)){
		haplo.std[i]<-sqrt(invScore.haplo[i,i])
	}

	l<-rep(1,nb-1)
	haplo.std[nb]<-sqrt(t(l)%*%invScore.haplo%*%l)
	param.std<-c(haplo.std, lambda.std)
	param.std
}

Ginv<-function (x, eps = 1e-06) 
{
    if (any(is.na(x))) 
        stop("NA(s) encountered in matrix, cannot invert.")
    if (length(x) > 1) {
        savesvd <- svd(x, LINPACK = FALSE)
        U.svd <- savesvd$u
        V.svd <- savesvd$v
        d.svd <- savesvd$d
        maxd <- max(d.svd)
        w <- ifelse((d.svd/maxd) < eps, rep(0, length(d.svd)), 
            1/d.svd)
        rank <- sum(d.svd/maxd >= eps)
        Ginv <- V.svd %*% diag(w) %*% t(U.svd)
    }
    else {
        Ginv <- ifelse(x < eps, 0, 1/x)
        rank <- ifelse(x < eps, 0, 1)
    }
    list(Ginv = Ginv, rank = rank)
}



num.unique.allele.locus<-function(x){
	max(unlist(lapply(x, num.unique.allele)))
}

num.unique.allele<-function(x){
	xx<-unlist(strsplit(x,split=" "))
	length(unique(xx))
}

malaria.em <- function(geno, sizes=c(2), locus.label=NA){
	if(length(sizes)<1){
		stop("Sizes information is incorrect! Length of sizes should be at least 1"); 
	}
	
	max.num.unique.allele<-max(apply(geno, 2,num.unique.allele.locus))

	if(max.num.unique.allele>max(sizes)){
		stop(paste("Sizes information is incorrect! The largest element in sizes should be at least", max.num.unique.allele, "!"));
	}

	time1<-Sys.time()	
	allPrefix<-c()
	nObs <- dim(geno)[1]
	NLoci<-dim(geno)[2]
	
	
	sg <- sortGenotypes(geno)
	genotypes <- unique (sg)
	genoCountList <- table (sg)
	n0<-nObs-sum(genoCountList)
	#cat("computing haplotypes......\n")
	allHaplotypes <- getAllHaplotypes(genotypes,allPrefix)
	nhaplo <- length(allHaplotypes)
	if(nhaplo>10){
	cat("The program is slow, please be patient......\n")	 
	}	

	#cat("numbe of haplotypes = ", nhaplo , "\n")
	#cat("computing haploset......\n")

	#Get the haplotype set for each unique genotype
 	allHaploSet.C <- genAllHaploSetList.C(allHaplotypes, genotypes, sizes, NLoci)
	#Get the haplotype set for each observation

  	genotypes.alldata<-toGenotypes(geno)
  	allHaploSet.C.alldata<-genAllHaploSetList.C.alldata(allHaploSet.C, genotypes, genotypes.alldata)
 
	 ##Get additive model matrix
  	xmat.additive <-getDesignforAllGenotypes.additive.allhaplotype(allHaploSet.C.alldata, allHaplotypes)
  	cvector<-getVector.C(xmat.additive)
	
	allQijCoefficients.C.alldata<-getAllCoefficientQij.C(allHaploSet.C.alldata)
	allEstiPs<-rep(1/nhaplo, times=nhaplo)
	EstiLambda<-0.1
	if( length(sizes)==1 ) {
		EstiLambda<- -1
	}
	allQij.C.alldata<-getAllQij.C(allQijCoefficients.C.alldata, allHaploSet.C.alldata, 	allEstiPs,EstiLambda,nObs,genoCountList)
	cat("Starting em......\n")	 
	ll<-logL(allHaploSet.C.alldata, allQij.C.alldata,EstiLambda,n0)
#	cat("logL=", ll, "\n")
	lllast<-ll
	diffL<-1000
	
	while( (diffL>0.00001) &(sum(allEstiPs==0)==0)){
		allEstiPslast<-allEstiPs
		EstiLambdalast<-EstiLambda 
		w<-getallWeight.C.em (allHaploSet.C.alldata,allQij.C.alldata)
		newnobs<-as.numeric(t(as.matrix(w))%*%cvector)
		allEstiPs <-t(as.matrix(w))%*%as.matrix(xmat.additive)/newnobs
		if( EstiLambda >0 )
		{
			if(n0>0){
				EstiLambda <- newnobs/ nObs
			}else{
				EstiLambda<-newton.glm(w, cvector, EstiLambda)
			}
		}
		allQij.C.alldata<-getAllQij.C(allQijCoefficients.C.alldata, allHaploSet.C.alldata, 	allEstiPs,EstiLambda,nObs,genoCountList)
		ll<-logL(allHaploSet.C.alldata, allQij.C.alldata,EstiLambda,n0)
#		cat("logL=", ll, "\n")
		diffL<-ll-lllast
		lllast<-ll
	}
	allEstiPs<-allEstiPslast
	EstiLambda<-EstiLambdalast 
	if( EstiLambda >0 )
	{
		if(n0>0){
			EstiNumInfection <- EstiLambda 
		}else{
			EstiNumInfection <- EstiLambda/(1-exp(-EstiLambda)) 
		}
	}else{
			EstiNumInfection<-sizes[1]
	}

	allQij.C.alldata<-getAllQij.C(allQijCoefficients.C.alldata, allHaploSet.C.alldata, 	allEstiPs,EstiLambda,nObs,genoCountList)
	w<-getallWeight.C.em(allHaploSet.C.alldata,allQij.C.alldata)
	
	allScore<-getAllScore.C.em(nhaplo,allHaploSet.C.alldata,allQij.C.alldata,allEstiPs,EstiLambda,n0)
	p.std<-getStd.em(allScore)

	if(EstiLambda>0){
		lambda.std<-p.std[length(p.std)]
	}else{
		EstiLambda<-lambda.std<-NA
	}

	ret<-getallConfig.C(allHaploSet.C.alldata,allQij.C.alldata)
	haplo.sets<-ret$haplo.sets
	pred.haplo.set<-ret$pred.haplo.set
	n.haplo.set<-ret$n.haplo.set
	haplotype<-matrix(NA, nrow=length(allHaplotypes), ncol=length(allHaplotypes[[1]]))
	for(i in 1:length(allHaplotypes)){
		haplotype[i,]<-allHaplotypes[[i]]
	}

	haplo.prob.tab<-cbind(haplotype,as.vector(allEstiPs), p.std[1:length(allHaplotypes)])
	
	if(is.na(locus.label)[1]){
	colnames(haplo.prob.tab)<-c(paste("loc-", rep(1:length(allHaplotypes[[1]])),sep=""),"hap.prob","hap.prob.std")
	}else{
			colnames(haplo.prob.tab)<-c(locus.label, "hap.prob","hap.prob.std")

	}

	est <- list(haplo.prob.tab=haplo.prob.tab,haplotype=haplotype, haplo.prob=as.vector(allEstiPs), 
			haplo.prob.std=p.std[1:length(allHaplotypes)],lambda=EstiLambda,lambda.std=lambda.std,
			NumofInfection=EstiNumInfection, haplo.sets=haplo.sets,n.haplo.set=n.haplo.set, pred.haplo.set=pred.haplo.set)
	class(est)<-"malaria.em"
	cat("Done......\n")
	est
}

print.malaria.em<-function(x,...){
	cat("Frenquency Estimate:\n")
	cat("======================================================================:\n")
	print(x$haplo.prob.tab)
	cat("======================================================================:\n")
	cat("Estimated Mean Number of Infections:\n")
	print(x$NumofInfection)
}

