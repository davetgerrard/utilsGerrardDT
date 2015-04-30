




testProfile <- c(0,8,0,0,0,2,0,2,0,0,0,0)  # tau = 0.95
testProfile2 <- c(0.8, 0.8, 0.8, 0.1, 0.8, 0.8, 0.8)  # tau = 0.14  inv.tau= 0.87
testProfile3 <- rep(0,10)

## Tau from Yanai et al. (2005) Bioinformatics.
tau <- function(profile)  {
  if(min(profile) == max(profile))  {
    tau <- 0 
  } else {
    tau <- (sum(1 - (profile/max(profile)))) / (length(profile) -1 )
  }
  return(tau)
} 

inv.tau <- function(profile)  {
    profile <- ifelse(profile==0,  1e-22, profile)  # cannot use zeroes
    if(min(profile) == max(profile))  {
      tau <- 0 
    } else { 
      tau <-  (sum(1 - (min(profile)/profile))) / (length(profile) -1 )
    }
    return(tau)
}

# takes a vector of grouping elements and first reduces the profile to a set of within-group summary stats using "method"
group.tau <- function(profile, groups, method=mean, use.inv=FALSE) {
  sub.profile <- tapply(profile, as.factor(as.character(groups)), method)
  
  if(min(profile) == max(profile))  {
    value <- 0
  } else {
    if(use.inv) {
      value <- inv.tau(sub.profile)
    } else {
      value <- tau(sub.profile)
    }
  }
  return(value)
}


(tau(testProfile) )


(tau(testProfile2) )
(inv.tau(testProfile2) )
(inv.tau(testProfile) )

(tau(testProfile3) )
(inv.tau(testProfile3) )

(tau(-testProfile2) )  # don't work well for negative numbers.
(inv.tau(-testProfile2) )

testProfile4 <- c(0,8,2,4,6,2,6)
tau(testProfile4)
group.tau(testProfile4, c("A", "A", "B", "B", "B", "C","C"))    # group means all== 4 , tau = 0


#testProfile <- rep(1,20)

# testProfile <- rep(0,20)



#testProfile <- c(0,2,0,0,0,2,0,2,0,0,0,0)

