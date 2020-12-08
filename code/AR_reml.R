
# Define time trend-finding function (created by Tony Ives)
AR_reml <- function(formula, data = list()) {

	AR_reml_funct <- function(par, x, u) {
		b <- par
		n.obs <- length(x)
		q <- dim(u)[2]
		B <- diag(n.obs)
		diag(B[-1, ]) <- -b

		iS <- diag(n.obs)
		iS[1, 1] <- (1 - b^2)
		iV <- t(B) %*% iS %*% B
		logdetV <- -determinant(iV)$modulus[1]

		coef <- solve(t(u) %*% iV %*% u, t(u) %*% iV %*% x)
		H <- x - u %*% coef

		s2 <- (t(H) %*% iV %*% H)/(n.obs - q)
		LL <- 0.5 * ((n.obs - q) * log(s2) + logdetV + determinant(t(u) %*% iV %*% u)$modulus[1] + (n.obs - q))
		#show(c(LL,b))
		return(LL)
	}


	mf <- model.frame(formula = formula, data = data)
	u <- model.matrix(attr(mf, "terms"), data = mf)
	x <- model.response(mf)

	q <- dim(u)[2]

	opt <- optim(fn = AR_reml_funct, par = 0.2, method = "Brent", upper = 1, lower = -1, control = list(maxit = 10^4), x = x, u = u)
	b <- opt$par

	n.obs <- length(x)
	q <- dim(u)[2]
	B <- diag(n.obs)
	diag(B[-1, ]) <- -b

	iS <- diag(n.obs)
	iS[1, 1] <- (1 - b^2)
	iV <- t(B) %*% iS %*% B
	logdetV <- -determinant(iV)$modulus[1]

	coef <- solve(t(u) %*% iV %*% u, t(u) %*% iV %*% x)
	H <- x - u %*% coef

	MSE <- as.numeric((t(H) %*% iV %*% H)/(n.obs - q))
	s2coef <- MSE * solve(t(u) %*% iV %*% u)
	Pr <- 1:q
	for (i in 1:q) Pr[i] <- 2 * pt(abs(coef[i])/s2coef[i, i]^0.5, df = n.obs - q, lower.tail = F)

	logLik <- 0.5 * (n.obs - q) * log(2 * pi) + determinant(t(u) %*% u)$modulus[1] - opt$value

	return(list(MSE = MSE, b = b, coef = coef, s2coef = s2coef, Pr = Pr, logLik = logLik))
}
