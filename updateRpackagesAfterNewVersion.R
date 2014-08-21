# after updating R, update packages that were used before.


# start a clean R session 
# On Windows "Run as administrator"? 
# DO NOT LOAD ANY PACKAGES - these will fail to update and may influence other packages which depend upon them.

# I use bioconductor and CRAN. The function biocinstallRepos() will return a vector of bioconductor repos AND CRAN.
# checkBuilt=TRUE tells R that packages installed in earlier versions of R should be added to the 'old' packages list.
# The bioc installer package needs to be updated first, so that it can decide which versions of bioconductor packages to install.

getOption("repos")

# set the bioconductor repo
source("http://bioconductor.org/biocLite.R")
biocinstallRepos()

biocLite("BiocUpgrade",ask=FALSE)		# upgrade bioconductor (and pre-installed packages)

# update all other packages.
update.packages(repos=biocinstallRepos(), checkBuilt=TRUE, ask=FALSE)


