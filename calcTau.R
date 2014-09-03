




testProfile <- c(0,8,0,0,0,2,0,2,0,0,0,0)  # tau = 0.95
testProfile2 <- c(0.8, 0.8, 0.8, 0.1, 0.8, 0.8, 0.8)  # tau = 0.14  inv.tau= 0.87

## Tau from Yanai et al. (2005) Bioinformatics.
tau <- function(profile)  {
	(sum(1 - (profile/max(profile)))) / (length(profile) -1 )
	
} 

inv.tau <- function(profile)  {
    (sum(1 - (min(profile)/profile))) / (length(profile) -1 )
    
}




(tau(testProfile) )


(tau(testProfile2) )
(inv.tau(testProfile2) )

#testProfile <- rep(1,20)

# testProfile <- rep(0,20)



#testProfile <- c(0,2,0,0,0,2,0,2,0,0,0,0)

