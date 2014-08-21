




testProfile <- c(0,8,0,0,0,2,0,2,0,0,0,0)  # tau = 0.95


## Tau from Yanai et al. (2005) Bioinformatics.
tau <- function(profile)  {
	(sum(1 - (profile/max(profile)))) / (length(profile) -1 )
	
} 


(tau(testProfile) )

testProfile <- rep(1,20)

# testProfile <- rep(0,20)



testProfile <- c(0,2,0,0,0,2,0,2,0,0,0,0)

