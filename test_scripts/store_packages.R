# store_packages.R
# from https://hlplab.wordpress.com/2012/06/01/transferring-installed-packages-between-different-installations-of-r/
# stores a list of your currently installed packages
.libPaths(c(.libPaths(), "C:/Users/Dave/Documents/R/win-library/3.1"))
tmp = installed.packages()

installedpackages = as.vector(tmp[is.na(tmp[,"Priority"]), 1])
save(installedpackages, file="~/Desktop/installed_packages.rda")