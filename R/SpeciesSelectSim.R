#' DirectSelectSim: Simulate population data on 3-tip species tree
#'
#' This function evolves population genomic data through a 3-tip species tree forward in time
#' @param theta theta value for all populations = 4*Ne*u.site
#' @param u.site mutation rate per site per generation
#' @param tau divergence time between species = tau/u.site
#' @param l number of sites in alignment
#' @param startp starting allele frequency of selected site in ancestor ABC
#' @keywords population genetics, coalescent models, multispecies forward-time
#' @export
#' @examples
#' y <- DirectSelectSim(theta = 0.01, u.site = 2*10^-8, tau = 0.0000001, l = 1000, startp = 0.5, w.list = c(1, 1, 0.9))
#' 
#' 

################################################
################ Mulitspecies simulation model
################################################
library(stringr)

DirectSelectSim<- function(N.diploid, u.site, l, ngen, startp, w.list){
	u.locus <- u.site*l
	wAA <- w.list[1]
	wAa <- w.list[2]
	waa <- w.list[3]
	ABC.freq <- c(startp^2,2*startp*(1-startp),(1-startp)^2)
	names(ABC.freq) <- c("A:A", "A:a", "a:a")
	ABC.pop <- rmultinom(1, size = N.diploid, prob = ABC.freq)/N.diploid
	AB.pop <- rmultinom(1, size = N.diploid, prob = ABC.pop[,1])/N.diploid
	C.pop <- rmultinom(1, size = N.diploid, prob = ABC.pop[,1])/N.diploid
	
	for (generation in 1:ngen){
		AB.pop <- rmultinom(1, size = N.diploid, prob = AB.pop[,1])/N.diploid
		C.pop <- rmultinom(1, size = N.diploid, prob = C.pop[,1])/N.diploid

		#mutate AB alleles
		mutate.AB <- runif(n = 2*N.diploid, min = 0, max = 1)
		mutated.AB <- mutate.AB[mutate.AB<u.locus]
		for (mutated.locus in 1:length(mutated.AB)){
			this <- sample(x = rownames(AB.pop), size = 1, prob = AB.pop)
			this.allele <- strsplit(x = this, split = ":")
			which.allele <- sample(x = c(1, 2), size = 2, replace = F)
			mutate.allele <- this.allele[[1]][which.allele[1]]
			same.allele <- this.allele[[1]][which.allele[2]]
			random.position = sample(x=seq(1, l, 1), size = 1)
			mutate.this.allele = paste(mutate.allele, random.position, sep = ".")
			if (which.allele[1] == 1){
				mutate.this.allele = paste(mutate.this.allele, same.allele, sep = ":")
			}
			if (which.allele[1] == 2){
				mutate.this.allele = paste(same.allele, mutate.this.allele, sep = ":")
			}
			f <- 1/N.diploid
			AB.pop <- rbind(AB.pop,f)
			rownames(AB.pop) <-  gsub(pattern = "f", replacement = mutate.this.allele, x = rownames(AB.pop))
			AB.pop[rownames(AB.pop) == this] <- ((AB.pop[rownames(AB.pop) == this]*N.diploid) - 1)/N.diploid
			AB.pop <- matrix(AB.pop[AB.pop>0], dimnames = list(rownames(AB.pop)[AB.pop > 0], 1))
		}
		
		#mutate C alleles
		mutate.C <- runif(n = 2*N.diploid, min = 0, max = 1)
		mutated.C <- mutate.C[mutate.C<u.locus]
		for (mutated.locus in 1:length(mutated.C)){
			this <- sample(x = rownames(C.pop), size = 1, prob = C.pop)
			this.allele <- strsplit(x = this, split = ":")
			which.allele <- sample(x = c(1, 2), size = 2, replace = F)
			mutate.allele <- this.allele[[1]][which.allele[1]]
			same.allele <- this.allele[[1]][which.allele[2]]
			random.position = sample(x=seq(1, l, 1), size = 1)
			mutate.this.allele = paste(mutate.allele, random.position, sep = ".")
			if (which.allele[1] == 1){
				mutate.this.allele = paste(mutate.this.allele, same.allele, sep = ":")
			}
			if (which.allele[1] == 2){
				mutate.this.allele = paste(same.allele, mutate.this.allele, sep = ":")
			}
			f <- 1/N.diploid
			C.pop <- rbind(C.pop,f)
			rownames(C.pop) <-  gsub(pattern = "f", replacement = mutate.this.allele, x = rownames(C.pop))
			C.pop[rownames(C.pop) == this] <- ((C.pop[rownames(C.pop) == this]*N.diploid) - 1)/N.diploid
			C.pop <- matrix(C.pop[C.pop>0], dimnames = list(rownames(C.pop)[C.pop > 0], 1))
		}
		
	}
	return(C.pop)
}


y <- DirectSelectSim(N.diploid = 125000, u.site = 2e-8, ngen = 10000, l = 1000, startp = 0.5, w.list = c(1, 1, 1))
sum(y)
