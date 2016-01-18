# MultiSpeciesSelectSim

Add on package to simulate genomic datasets with varying degree and strengths of loci under selection.

1) Uses the MCcoal software provided in BPP to simulat m-x loci under the neutral coalescent
2) Uses a forward-time algorithm to simulate x loci under selective pressure


- For the x selected loci, we apply a Wright-Fisher forward time algorithm to simulate loci under selective pressure through a species tree that is parameterized by theta and divergence times
- These parameters are converted in to disrecte variables (divergence time/mutation rate = # generations)(theta = 4*mutationrate*#diploid)
