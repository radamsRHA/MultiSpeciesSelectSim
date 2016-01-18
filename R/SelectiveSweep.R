#' SelectiveSweep: Simulate natural selection at removing, adding individuals at random
#'
#'	This function returns a multinomial population object after a selective event 
#' @param population population genotype frequencies in multinomial format
#' @param w.list list of selection relative fitness
#' @param N. diploid # diploid individials
#' @keywords population genetics, coalescent models, multispecies forward-time
#' @export
#' @examples
#' DirectSelectSim(theta = 0.01, u.site = 2*10^-8, tau = 0.00001, l = 1000, startp = 0.5)
#' y <- SelectiveSweep(population = AB.pop, w.list = c(.9, .9, 0.999), N.diploid = 125000)
#' sum(y) <- 0

################################################
################ Mulitspecies simulation model
################################################


SelectiveSweep <- function(population, w.list, N.diploid){
	wAA <- w.list[1]
	wAa <- w.list[2]
	waa <- w.list[3]
	
	#Frequencies before Selection
	fAA <- sum(population[str_count(rownames(population), pattern = "A") == 2])
	fAa <- sum(population[str_count(rownames(population), pattern = "A") == 1])
	faa <- sum(population[str_count(rownames(population), pattern = "A") == 0])
	
	#Frequencies after selection
	w.bar <- fAA*wAA+fAa*wAa+faa*waa
	fAA.prime <- fAA*wAA/w.bar
	fAa.prime <- fAa*wAa/w.bar
	faa.prime <- faa*waa/w.bar
	
	# Selected individuals
	diff.AA <- round(((fAA.prime*N.diploid)-(fAA*N.diploid)), digits = 0)
	diff.Aa <- round(((fAa.prime*N.diploid)-(fAa*N.diploid)), digits = 0)
	diff.aa <- round(((faa.prime*N.diploid)-(faa*N.diploid)), digits = 0)
	
	diff.list <- c(diff.AA, diff.Aa, diff.aa)
	if (diff.AA != 0){
		AA.pop <- population[str_count(rownames(population), pattern = "A") == 2]*N.diploid
		names(AA.pop) <- rownames(population)[str_count(rownames(population), pattern = "A") == 2]
		selected.sample <- sample(x = names(AA.pop), size = abs(diff.AA), replace = T, prob = AA.pop/sum(AA.pop))
		for (i in unique(selected.sample)){
			n.selected <- length(selected.sample[selected.sample == i])
			if (diff.AA < 0){
				n.selected <- -n.selected
			}
			population[rownames(population) == i] <- (population[rownames(population) == i]*N.diploid+n.selected)/N.diploid
		}
	}
	if (diff.Aa != 0){
		Aa.pop <- population[str_count(rownames(population), pattern = "A") == 1]*N.diploid
		names(Aa.pop) <- rownames(population)[str_count(rownames(population), pattern = "A") == 1]
		selected.sample <- sample(x = names(Aa.pop), size = abs(diff.Aa), replace = T, prob = Aa.pop/sum(Aa.pop))
		for (i in unique(selected.sample)){
			n.selected <- length(selected.sample[selected.sample == i])
			if (diff.Aa < 0){
				n.selected <- -n.selected
			}
			population[rownames(population) == i] <- (population[rownames(population) == i]*N.diploid+n.selected)/N.diploid
		}		
	}
	if (diff.aa != 0){
		aa.pop <- population[str_count(rownames(population), pattern = "A") == 0]*N.diploid
		names(aa.pop) <- rownames(population)[str_count(rownames(population), pattern = "A") == 0]
		selected.sample <- sample(x = names(aa.pop), size = abs(diff.aa), replace = T, prob = aa.pop/sum(aa.pop))
		for (i in unique(selected.sample)){
			n.selected <- length(selected.sample[selected.sample == i])
			if (diff.aa < 0){
				n.selected <- -n.selected
			}
			population[rownames(population) == i] <- (population[rownames(population) == i]*N.diploid+n.selected)/N.diploid
		}		
	}
	
	return(population)
}

y <- SelectiveSweep(population = AB.pop, w.list = c(1,1,0.9), N.diploid = 125000)
sum(y)
