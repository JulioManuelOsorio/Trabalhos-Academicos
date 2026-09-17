##########################################################################################
#
# reorder_fund_shks.R
#
# Code to reorder the fundamental shocks if required
# Philip Barrett, Washington DC
# First version: 20jan2024
#
##########################################################################################


m.A <- est.A$A
reorder.from.idx <- n.fcast + 1
reorder.to.idx <- n.fcast * 2 + which(y==order.first)

m.trans <- diag(nrow(m.A)) ; 
m.trans[reorder.from.idx,reorder.from.idx] <- m.trans[reorder.to.idx,reorder.to.idx] <- 0 ; 
m.trans[reorder.from.idx,reorder.to.idx] <- m.trans[reorder.to.idx,reorder.from.idx] <- 1
m.trans.sub <- m.trans[-(1:n.fcast),-(1:n.fcast)]

m.Sigma <- l.var.coefs$Sigma
m.Sigma.trans <- m.trans %*% m.Sigma %*% t(m.trans)
  # The target for the outcome
m.A.trans <- m.trans %*% m.A
m.A.trans.check <- m.A.trans %*% t(m.A.trans) - m.Sigma.trans
  # The initial transposed A matrix

m.D <- m.A.trans[-(1:n.fcast),-(1:n.fcast)]
m.D.sq <- m.D %*% t(m.D)
m.D.tilde <- chol( m.D.sq ) %>% t()
m.D.tilde.check <- m.D.sq - ( m.D.tilde %*% t(m.D.tilde) )
  # Partial out the new D
m.B <- m.A.trans[1,-1]
m.B.tilde <- (solve( m.D.tilde ) %*% m.D %*% m.B) %>% t
m.B.tilde.check <- m.B.tilde ^2 %>% sum() - m.B ^2 %>% sum()

m.A.alt <- m.A.trans ; m.A.alt[-(1:n.fcast),-(1:n.fcast)] <- m.D.tilde ; m.A.alt[1,-1] <- m.B.tilde
m.Sigma.trans.check <- m.A.alt %*% t(m.A.alt) - m.Sigma.trans
  
m.A.final <- m.trans %*% m.A.alt %>%
  set_colnames(colnames(m.A)) %>%
  set_rownames(rownames(m.A))
m.A.final.check <- (m.A.final %*% t(m.A.final) - m.Sigma) %>% abs %>% max
est.A <- list( A=m.A.final, sigma.err=m.A.final.check ) 
