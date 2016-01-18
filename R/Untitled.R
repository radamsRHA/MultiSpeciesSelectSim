#library(stringr)
#library(compiler)

DirectSelectSim<- function(theta, u.site, tau, l, startp, w.list){
	u.locus <- u.site*l
	N.diploid <- theta/(4*u.site)
	ngen <- tau/u.site
	wAA <- w.list[1]
	wAa <- w.list[2]
	waa <- w.list[3]
	ABC.freq <- c(startp^2,2*startp*(1-startp),(1-startp)^2)
	names(ABC.freq) <- c("A:A", "A:a", "a:a")
	ABC.pop <- rmultinom(1, size = N.diploid, prob = ABC.freq)/N.diploid
	AB.pop <- rmultinom(1, size = N.diploid, prob = ABC.pop[,1])/N.diploid

	for (generation in 1:ngen){
		AB.pop <- rmultinom(1, size = N.diploid, prob = AB.pop[,1])/N.diploid
		
		#Frequencies before Selection
		fAA <- AB.pop[str_count(rownames(AB.pop), pattern = "A") == 2]
		fAa <- AB.pop[str_count(rownames(AB.pop), pattern = "A") == 1]
		faa <- AB.pop[str_count(rownames(AB.pop), pattern = "A") == 0]

		#Frequencies after selection
		w.bar <- fAA*wAA+fAa*wAa+faa*waa
		fAA.prime <- fAA*wAA/w.bar
		fAa.prime <- fAa*wAa/w.bar
		faa.prime <- faa*waa/w.bar
		
		diff.AA <- round(((fAA.prime*N.diploid)-(fAA*N.diploid)), digits = 0)
		diff.Aa <- round(((fAa.prime*N.diploid)-(fAa*N.diploid)), digits = 0)
		diff.aa <- round(((faa.prime*N.diploid)-(faa*N.diploid)), digits = 0)
		
		
		AA.pop <- AB.pop[str_count(rownames(AB.pop), pattern = "A") == 2]*N.diploid
		names(AA.pop) <- rownames(AB.pop)[str_count(rownames(AB.pop), pattern = "A") == 2]
		
		print(remove.individuals)

		#mutate AB alleles
		#mutate.AB <- runif(n = 2*N.diploid, min = 0, max = 1)
		#mutated.AB <- mutate.AB[mutate.AB<u.locus]
		#mutated.copies <- match(mutated.AB, mutate.AB)
		#for (mutated.locus in mutated.copies){
	#		mutate.this.allele = AB.pop[mutated.locus]
	#		random.position = sample(x=seq(2, 1000, 1), size = 1)
	#		mutate.this.allele = paste(mutate.this.allele, random.position, sep = ".")
	#		AB.pop[mutated.locus] <- mutate.this.allele
	#	}
	}
	return(fAA)
}

AA.pop <- AB.pop[str_count(rownames(AB.pop), pattern = "A") == 2]*N.diploid
names(AA.pop) <- rownames(AB.pop)[str_count(rownames(AB.pop), pattern = "A") == 2]
chosen.allele <- sample(x = names(AA.pop), size = 0, replace = F, prob = AA.pop/sum(AA.pop))
AA.pop+chosen.allele

y <- DirectSelectSim(theta = 0.01, u.site = 2*10^-8, tau = 0.0000001, l = 1000, startp = 0.5, w.list = c(0.999, 1, 1))
