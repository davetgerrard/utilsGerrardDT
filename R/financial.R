

monthly <- function(principal, annual.interest, years=25, payments=12*years, monthly.interest=annual.interest/12) {
  
  payment <- principal * ( (monthly.interest * ((1 + monthly.interest)^payments))  / (((1+monthly.interest)^payments)-1))
  return(payment)
}


monthly(100000, 0.03)
monthly(100000, 0.05)
monthly(200000, 0.05)
