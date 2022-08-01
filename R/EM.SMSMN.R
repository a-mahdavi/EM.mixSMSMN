EM.mixSMSMN <- function(y, M=1, w=1, xi=NULL, S=NULL, la=NULL, nu=NULL,
 family="SMSTN.mix", get.init = TRUE, iter.max=50, tol=10^-6, eqnu=FALSE, CML=TRUE, group=TRUE){
 # M: number of components to fit
 # w: the vector of components probabilites.
 # xi: the list of vector of location parameters.
 # S: the list of cov-variance matrices.
 # la: the list of vector of shape parameters.
 # nu: the list of flatness parameters.
 # family: distribution family to be used in fitting ("SMSTN.mix", "SMSSLN.mix",
 # "SMSCN.mix").
 # get.init: if TRUE, the initial values are generated.
 # iter.max: the maximum number of iterations of the EM algorithm. Default= 50.
 # tol: the covergence maximum error. Default= 10^-6.
 # eqnu: if TRUE the DF parameters assumed to be equal.
 # group: if TRUE it returns the group id of each observation
	require(orthoDr)
	if ((family != "SMSTN.mix") && (family != "SMSSLN.mix") && (family != 
        "SMSCN.mix") )  stop(paste("Family", family, "not recognized.\n", sep = " "))
	y <- as.matrix(y) ; p <- ncol(y); n <- dim(y)[1]
	if ( p <= 1 ) 
        stop("This function is developed for multivariate case.\n ")
	EM.SMSTN.mix <- function(y, M, w=NULL, xi=NULL, S=NULL, la=NULL, nu=NULL, iter.max=50, tol=10^-6,get.init = TRUE, CML=TRUE, group=FALSE, eqnu=FALSE){
	 begin <- proc.time()[3] 
	MSTN=function(y, xi, S, la, nu){
	y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
	S <- as.matrix(S) ; xi.mat <- matrix(xi,n,p,byrow=T)
	eig.S = eigen(S) ; 	y.xi <- y-xi.mat 
	D <- eig.S$vec ; A=eig.S$val
	f <- A^(-1/2)*dt(apply(y.xi,1,function(x)A^(-1/2)*t(D)%*%x),nu)
	2*apply(f,2,function(x) prod(x))*pnorm( apply(y.xi,1,function(x) t(la)%*%solve(D%*%diag(A^(1/2)))%*%x))
		 }
	SMSTN.mix=function(y, M, w, xi, S, la, nu){
	dens <- 0
	for(m in 1:M)
	dens <- dens + w[m]*MSTN(y, xi[[m]], S[[m]], la[[m]], nu[[m]])
	return(dens)
	}
	y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
        dif <- 1 ;  count <- 0 ; 	LL <- 1 
	if (get.init == TRUE) {
               init <- kmeans(y, M,  algorithm="Hartigan-Wong")
               w <- init$size/nrow(y)
               xi <- as.vector(init$centers)
               s <- sqrt(init$withinss/init$size)	
			xi <-  la <- S <- nu <- list() 
			for(m in 1:M){
			xi [[m]] <- init$centers[m, ]
			S[[m]] <- var(y[init$cluster == m, ])
              	 la[[m]] <-  sign(apply((y[init$cluster == 
                    m, ] - matrix(rep(xi[[m]], nrow(y[init$cluster == 
                    m, ])), nrow = nrow(y[init$cluster == m, 
                    ]), ncol = p, byrow = TRUE))^3, 2, sum))
				if(eqnu)
				nu[[m]] <- ceiling(runif(1,0,9))
				else
              	 nu[[m]] <- ceiling(runif(p,0,9))
          		}     }
	while ((dif > tol) && (count <= iter.max)) {
	# E step
	gam.hat <-s.hat <- array(0,c(n,p,M)); z.hat <- W.hat <- matrix(0,n,M)
	Del.gam <- array(0, c(p,p,n,M)) ; SMSTN.total <- SMSTN.mix(y,M,w,xi,S,la,nu)
	for (m in 1:M){
	z.hat[,m] <- w[m]*MSTN(y,xi[[m]],S[[m]],la[[m]],nu[[m]])/SMSTN.total
	Sm <- as.matrix(S[[m]]) ; xi.mat <- matrix(xi[[m]],n,p,byrow=T)
	y.xi <- y-xi.mat 
	eig.S = eigen(Sm)
	D <- eig.S$vec ; A=diag(eig.S$val); P=D%*%(A^(1/2))
	alpha.T <- t(la[[m]])%*%solve(P); alpha <- t(alpha.T); nu.m <- nu[[m]]
	if (eqnu)
	nu.m <- rep(nu.m,p)
	for (i in 1:n){
	W.hat[i,m] <- as.numeric(z.hat[i,m]*alpha.T%*%y.xi[i,])
	if (pnorm(alpha.T%*%y.xi[i,])!=0)
	W.hat[i,m] <- W.hat[i,m]+ z.hat[i,m]*dnorm(alpha.T%*%y.xi[i,])/pnorm(alpha.T%*%y.xi[i,])
	for(j in 1:p){
	del2 <- (A[j,j]^(-1/2)*(t(D)%*%y.xi[i,])[j])^2
	gam.hat[i,j,m] <- z.hat[i,m]*(nu.m[j]+1)/(nu.m[j]+del2) 
	s.hat[i,j,m] <- z.hat[i,m]*(digamma((nu.m[j]+1)/2)-log((nu.m[j]+del2)/2))}
	Del.gam[,,i,m] <- diag(as.vector(gam.hat[i,,m]))
	}
      # MCE steps
	w[m] <- sum(z.hat[,m])/n
	a <- 0 ; b <- 0
	for(i in 1:n){
	a <- a+ solve(t(P))%*%Del.gam[,,i,m]%*%solve(P)+z.hat[i,m]*alpha%*%alpha.T
	b <- b + solve(t(P))%*%Del.gam[,,i,m]%*%solve(P)%*%y[i,]+z.hat[i,m]*alpha%*%alpha.T%*%y[i,]-alpha*W.hat[i,m]
	}
	xi[[m]] <- solve(a)%*%b
	xi.mat <- matrix(xi[[m]],n,p,byrow=T)
	y.xi <- y-xi.mat
	aa <- 0
	for(i in 1:n)
	aa <- aa + Del.gam[,,i,m]%*%t(D)%*%y.xi[i,]%*%t(y.xi[i,])%*%D
	A <- diag(diag(aa),p,p)/sum(z.hat[,m])
	D <- ortho_optim(D,function(B){
	c <- 0
	for(i in 1:n)
	c <- c + t(y.xi[i,])%*%B%*%Del.gam[,,i,m]%*%solve(A)%*%t(B)%*%y.xi[i,]
	return(c)
	},maxitr = 200)$B
	 P=D%*%(A^(1/2))
	S[[m]] <-  D%*%A%*%t(D)
	cc <- 0 ; dd <- 0
	for(i in 1:n){
	cc <- cc + z.hat[i,m]*y.xi[i,]%*%t(y.xi[i,])
	dd <- dd + W.hat[i,m]*y.xi[i,] }
	alpha <- solve(cc)%*% dd
	alpha.T <- t(alpha)
	la[[m]] <- t(P)%*%alpha
	if (CML==F && eqnu==F){
	for(j in 1:p) 
	nu[[m]][j] <- uniroot(function(x) (log(x/2)+1)-digamma(x/2)+sum(s.hat[,j,m]-gam.hat[,j,m])/sum(z.hat[,m])
	,c(.001,100),maxiter = 10000)$root
 	}	}
	if(eqnu){
		nu0 <- nu
	nu <- optim(unlist(nu),function(x){
		for (m in 1:M)
		nu0[[m]] <- x[m]
		-sum(log(SMSTN.mix(y,M,w, xi, S, la, nu0))) 
		},method="L-BFGS-B",lower=.1,upper=100)$par
				} else{
	nu0 <- nu
	nu <- optim(unlist(nu),function(x){
		x <- cbind(x,rep(1:M,each=p)) 
		sp <- split(x[,1], x[,2])
		for (m in 1:M)
		nu0[[m]] <- sp[[m]]
		-sum(log(SMSTN.mix(y,M,w, xi, S, la, nu0))) 
		},method="L-BFGS-B",lower=.1,upper=100)$par
			nu <- cbind(nu,rep(1:M,each=p)) 
			sp <- split(nu[,1], nu[,2])
		for (m in 1:M)
		nu0[[m]] <- sp[[m]]
		nu <- nu0
		}
	LL.new <- sum(log(SMSTN.mix(y,M, w, xi, S, la, nu))) # log-likelihood function
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	LL <- LL.new
	}
	if(eqnu){
	aic <- -2 * LL.new + 2 * (M*(2*p+p*(p+1)/2+1)+M-1)
	bic <- -2 * LL.new + log(n) * (M*(2*p+p*(p+1)/2+1)+M-1)
	} else{
	aic <- -2 * LL.new + 2 * (M*(3*p+p*(p+1)/2)+M-1)
	bic <- -2 * LL.new + log(n) * (M*(3*p+p*(p+1)/2)+M-1) }
	end <- proc.time()[3]
	time <- end-begin
	obj.out <- list(w=w, xi=xi, S=S, la=la , nu=nu, loglik=LL.new, AIC=aic, BIC=bic, iter=count,elapsed=as.numeric(time),group = apply(z.hat, 
                1, which.max))
	if (group==FALSE)
	obj.out <- obj.out[-length(obj.out)]
	obj.out
	}
	EM.SMSSLN.mix <- function(y, M, w, xi, S, la, nu, iter.max=50, tol=10^-6, get.init=TRUE, CML=TRUE, eqnu=F, group=FALSE){
 	begin <- proc.time()[3] 
	f.s <- function(x,nu){
	if (norm(x)==0)
		nu/((nu+1/2)*sqrt(2*pi))
	else
		nu*gamma(nu+1/2)*2^nu/(sqrt(pi))*pgamma(x^2/2,nu+1/2)/abs(x)^(2*nu+1)
		}
	MSSL=function(y, xi, S, la, nu){
	y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
	S <- as.matrix(S) ; xi.mat <- matrix(xi,n,p,byrow=T)
	eig.S = eigen(S) ; 	y.xi <- y-xi.mat 
	D <- eig.S$vec ; A=eig.S$val
	f <- A^(-1/2)*f.s(apply(y.xi,1,function(x)A^(-1/2)*t(D)%*%x),nu)
	2*apply(f,2,function(x) prod(x))*pnorm( apply(y.xi,1,function(x) t(la)%*%solve(D%*%diag(A^(1/2)))%*%x))
		 }
	SMSSLN.mix=function(y,M, w, xi, S, la, nu){
		 dens <- 0
	for(m in 1:M)
	dens <- dens + w[m]*MSSL(y, xi[[m]], S[[m]], la[[m]], nu[[m]])
	return(dens)
	}
	y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
        dif <- 1 ;  count <- 0 ; 	LL <- 1 
	if (get.init == TRUE) {
               init <- kmeans(y, M,  algorithm="Hartigan-Wong")
               w <- init$size/nrow(y)
               xi <- as.vector(init$centers)
               s <- sqrt(init$withinss/init$size)	
			xi <-  la <- S <- nu <- list() 
			for(m in 1:M){
			xi [[m]] <- init$centers[m, ]
			S[[m]] <- var(y[init$cluster == m, ])
              	 la[[m]] <-  sign(apply((y[init$cluster == 
                    m, ] - matrix(rep(xi[[m]], nrow(y[init$cluster == 
                    m, ])), nrow = nrow(y[init$cluster == m, 
                    ]), ncol = p, byrow = TRUE))^3, 2, sum))
			if(eqnu)
				nu[[m]] <- ceiling(runif(1,0,9))
				else
              	 nu[[m]] <-  ceiling(runif(p,0,9))
          		}     }
	while ((dif > tol) && (count <= iter.max)) {
	# E step
	gam.hat <-s.hat <- array(0,c(n,p,M)); z.hat <- W.hat <- matrix(0,n,M)
	Del.gam <- array(0, c(p,p,n,M)) ; SMSSLN.total <- SMSSLN.mix(y,M,w,xi,S,la,nu)
	for (m in 1:M){
	z.hat[,m] <- w[m]*MSSL(y,xi[[m]],S[[m]],la[[m]],nu[[m]])/SMSSLN.total
	Sm <- as.matrix(S[[m]]) ; xi.mat <- matrix(xi[[m]],n,p,byrow=T)
	y.xi <- y-xi.mat 
	eig.S = eigen(Sm)
	D <- eig.S$vec ; A=diag(eig.S$val); P=D%*%(A^(1/2))
	alpha.T <- t(la[[m]])%*%solve(P); alpha <- t(alpha.T); nu.m <- nu[[m]]
		if (eqnu)
	nu.m <- rep(nu.m,p)
	for (i in 1:n){
	W.hat[i,m] <- as.numeric(z.hat[i,m]*alpha.T%*%y.xi[i,])
	if (pnorm(alpha.T%*%y.xi[i,])!=0)
	W.hat[i,m] <- W.hat[i,m]+ z.hat[i,m]*dnorm(alpha.T%*%y.xi[i,])/pnorm(alpha.T%*%y.xi[i,])
	for(j in 1:p){
	del2 <- (A[j,j]^(-1/2)*(t(D)%*%y.xi[i,])[j])^2
	if (del2==0)
	{
	gam.hat[i,j,m] <- z.hat[i,m]*(2*nu.m[j]+1)/(2*nu.m[j]+3)
	s.hat[i,j,m] <- -1/(nu.m[j]+1/2)*z.hat[i,m]} else{
	gam.hat[i,j,m] <- z.hat[i,m]*(2*nu.m[j]+1)/(del2)*pgamma(del2/2,nu.m[j]+3/2)/pgamma(del2/2,nu.m[j]+1/2) 
	s.hat[i,j,m] <- z.hat[i,m]*(log(2/del2)+integrate(function(x)log(x)*x^(nu.m[j]-1/2)*exp(-x),lower=0,upper=del2/2)$value/(gamma(nu.m[j]+1/2)*pgamma(del2/2,nu.m[j]+1/2)))
	}}	
	Del.gam[,,i,m] <- diag(as.vector(gam.hat[i,,m]))
	}
	# MCE steps
	w[m] <- sum(z.hat[,m])/n
	a <- 0 ; b <- 0
	for(i in 1:n){
	a <- a+ solve(t(P))%*%Del.gam[,,i,m]%*%solve(P)+z.hat[i,m]*alpha%*%alpha.T
	b <- b + solve(t(P))%*%Del.gam[,,i,m]%*%solve(P)%*%y[i,]+z.hat[i,m]*alpha%*%alpha.T%*%y[i,]-alpha*W.hat[i,m]
	}
	xi[[m]] <- solve(a)%*%b
	xi.mat <- matrix(xi[[m]],n,p,byrow=T)
	y.xi <- y-xi.mat
	aa <- 0
	for(i in 1:n)
	aa <- aa + Del.gam[,,i,m]%*%t(D)%*%y.xi[i,]%*%t(y.xi[i,])%*%D
	A <- diag(diag(aa),p,p)/sum(z.hat[,m])
	D <- ortho_optim(D,function(B){
	c <- 0
	for(i in 1:n)
	c <- c + t(y.xi[i,])%*%B%*%Del.gam[,,i,m]%*%solve(A)%*%t(B)%*%y.xi[i,]
	return(c)
	},maxitr = 200)$B
	 P=D%*%(A^(1/2))
	S[[m]] <-  D%*%A%*%t(D)
	cc <- 0 ; dd <- 0
	for(i in 1:n){
	cc <- cc + z.hat[i,m]*y.xi[i,]%*%t(y.xi[i,])
	dd <- dd + W.hat[i,m]*y.xi[i,] }
	alpha <- solve(cc)%*% dd
	alpha.T <- t(alpha)
	la[[m]] <- t(P)%*%alpha
	if (CML==F && eqnu==F){
	for(j in 1:p) 
	nu[[m]][j] <- -sum(z.hat[,m])/sum(s.hat[,j,m]) 	
	}	}
	if(eqnu){
		nu0 <- nu
	nu <- optim(unlist(nu),function(x){
		for (m in 1:M)
		nu0[[m]] <- x[m]
		-sum(log(SMSSLN.mix(y,M,w, xi, S, la, nu0))) 
		},method="L-BFGS-B",lower=.1,upper=100)$par
				} else{
	nu0 <- nu
	nu <- optim(unlist(nu),function(x){
		x <- cbind(x,rep(1:M,each=p)) 
		sp <- split(x[,1], x[,2])
		for (m in 1:M)
		nu0[[m]] <- sp[[m]]
		-sum(log(SMSSLN.mix(y,M,w, xi, S, la, nu0))) 
		},method="L-BFGS-B",lower=0.1,upper=100)$par
			nu <- cbind(nu,rep(1:M,each=p)) 
			sp <- split(nu[,1], nu[,2])
		for (m in 1:M)
		nu0[[m]] <- sp[[m]]
	nu <- nu0
		}
	LL.new <- sum(log(SMSSLN.mix(y,M, w, xi, S, la, nu))) # log-likelihood function
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	LL <- LL.new
	}
	if(eqnu){
	aic <- -2 * LL.new + 2 * (M*(2*p+p*(p+1)/2+1)+M-1)
	bic <- -2 * LL.new + log(n) * (M*(2*p+p*(p+1)/2+1)+M-1)
	} else{

	aic <- -2 * LL.new + 2 * (M*(3*p+p*(p+1)/2)+M-1)
	bic <- -2 * LL.new + log(n) * (M*(3*p+p*(p+1)/2)+M-1) }
	end <- proc.time()[3]
	time <- end-begin
	obj.out <- list(w=w, xi=xi, S=S, la=la , nu=nu, loglik=LL.new, AIC=aic, BIC=bic, iter=count,elapsed=as.numeric(time),group = apply(z.hat, 
                1, which.max))
	if (group==FALSE)
	obj.out <- obj.out[-length(obj.out)]
	obj.out
	}

	EM.SMSCN.mix <- function(y, M, w, xi, S, la, nu, iter.max=50, tol=10^-6, get.init=TRUE, eqnu=FALSE, group=FALSE){
 	begin <- proc.time()[3] 
	MSCN=function(y, xi, S, la, nu){
	y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
	S <- as.matrix(S) ; xi.mat <- matrix(xi,n,p,byrow=T)
	eig.S = eigen(S) ; 	y.xi <- y-xi.mat 
	nu <- matrix(unlist(nu),p,2); nu1<-nu[,1] ; nu2=nu[,2]
	D <- eig.S$vec ; A=eig.S$val
	f <- A^(-1/2)*( nu1*nu2^(1/2)*dnorm(apply(y.xi,1,function(x)nu2^(1/2)*A^(-1/2)*t(D)%*%x))+(1-nu1)*dnorm(apply(y.xi,1,function(x)A^(-1/2)*t(D)%*%x)))
	2*apply(f,2,function(x) prod(x))*pnorm( apply(y.xi,1,function(x) t(la)%*%solve(D%*%diag(A^(1/2)))%*%x))
	}
	SMSCN.mix=function(y,M, w, xi, S, la, nu){
	 dens <- 0
	for(m in 1:M)
	dens <- dens + w[m]*MSCN(y, xi[[m]], S[[m]], la[[m]], nu[[m]])
	return(dens)
	}
	y <- as.matrix(y) ; n <- dim(y)[1] ; p <- dim(y)[2]
        dif <- 1 ;  count <- 0 ; 	LL <- 1 
	if (get.init == TRUE) {
               init <- kmeans(y, M,  algorithm="Hartigan-Wong")
               w <- init$size/nrow(y)
               xi <- as.vector(init$centers)
               s <- sqrt(init$withinss/init$size)	
			xi <-  la <- S <- list() ; nu <- list(list()) 
			for(m in 1:M){
			xi [[m]] <- init$centers[m, ]
			S[[m]] <- var(y[init$cluster == m, ])
              	 la[[m]] <-  sign(apply((y[init$cluster == 
                    m, ] - matrix(rep(xi[[m]], nrow(y[init$cluster == 
                    m, ])), nrow = nrow(y[init$cluster == m, 
                    ]), ncol = p, byrow = TRUE))^3, 2, sum))
			if(eqnu)
			 nu[[m]] <- list(runif(2,0,1))
			else
			for(j in 1:p)
              	 nu[[m]][j] <- list(runif(2,0,1))
          		}     }
	while ((dif > tol) && (count <= iter.max)) {
	# E step
	gam.hat <-s.hat <- array(0,c(n,p,M)); z.hat <- W.hat <- matrix(0,n,M)
	Del.gam <- array(0, c(p,p,n,M)) ; SMSCN.total <- SMSCN.mix(y,M,w,xi,S,la,nu)
	for (m in 1:M){
	z.hat[,m] <- w[m]*MSCN(y,xi[[m]],S[[m]],la[[m]],nu[[m]])/SMSCN.total
	Sm <- as.matrix(S[[m]]) ; xi.mat <- matrix(xi[[m]],n,p,byrow=T)
	y.xi <- y-xi.mat 
	eig.S = eigen(Sm)
	D <- eig.S$vec ; A=diag(eig.S$val); P=D%*%(A^(1/2))
	alpha.T <- t(la[[m]])%*%solve(P); alpha <- t(alpha.T); nu.m <- nu[[m]]
	if (eqnu)
	nu.m <- rep(nu.m,p)
	for (i in 1:n){
	W.hat[i,m] <- as.numeric(z.hat[i,m]*alpha.T%*%y.xi[i,])
	if (pnorm(alpha.T%*%y.xi[i,])!=0)
	W.hat[i,m] <- W.hat[i,m]+ z.hat[i,m]*dnorm(alpha.T%*%y.xi[i,])/pnorm(alpha.T%*%y.xi[i,])
	for(j in 1:p){
	nuj1 <- as.vector(nu.m[[j]])[1] ; nuj2 <- as.vector(nu.m[[j]])[2]
	del2 <- (A[j,j]^(-1/2)*(t(D)%*%y.xi[i,])[j])^2
	gam.hat[i,j,m] <- z.hat[i,m]*(1-nuj1+nuj1*nuj2^(3/2)*exp((1-nuj2)*del2/2))/(1-nuj1+nuj1*nuj2^(1/2)*exp((1-nuj2)*del2/2))
	}
	Del.gam[,,i,m] <- diag(as.vector(gam.hat[i,,m]))
	}
	# MCE steps
	w[m] <- sum(z.hat[,m])/n
	a <- 0 ; b <- 0
	for(i in 1:n){
	a <- a+ solve(t(P))%*%Del.gam[,,i,m]%*%solve(P)+z.hat[i,m]*alpha%*%alpha.T
	b <- b + solve(t(P))%*%Del.gam[,,i,m]%*%solve(P)%*%y[i,]+z.hat[i,m]*alpha%*%alpha.T%*%y[i,]-alpha*W.hat[i,m]
	}
	xi[[m]] <- solve(a)%*%b
	xi.mat <- matrix(xi[[m]],n,p,byrow=T)
	y.xi <- y-xi.mat
	aa <- 0
	for(i in 1:n)
	aa <- aa + Del.gam[,,i,m]%*%t(D)%*%y.xi[i,]%*%t(y.xi[i,])%*%D
	A <- diag(diag(aa),p,p)/sum(z.hat[,m])
	D <- ortho_optim(D,function(B){
	c <- 0
	for(i in 1:n)
	c <- c + t(y.xi[i,])%*%B%*%Del.gam[,,i,m]%*%solve(A)%*%t(B)%*%y.xi[i,]
	return(c)
	},maxitr = 200)$B
	 P=D%*%(A^(1/2))
	S[[m]] <-  D%*%A%*%t(D)
	cc <- 0 ; dd <- 0
	for(i in 1:n){
	cc <- cc + z.hat[i,m]*y.xi[i,]%*%t(y.xi[i,])
	dd <- dd + W.hat[i,m]*y.xi[i,] }
	alpha <- solve(cc)%*% dd
	alpha.T <- t(alpha)
	la[[m]] <- t(P)%*%alpha
		}
	if(eqnu){
		nu0 <-nu.0 <- nu 
	nu <- optim(unlist(nu),function(x){
		for (m in 1:M)
		nu0[[m]] <- c(x[2*m-1],x[2*m])
		-sum(log(SMSCN.mix(y,M,w, xi, S, la, nu0))) 
		},method="L-BFGS-B",lower=.01,upper=.99)$par
		for (m in 1:M)
		nu.0[[m]][[1]] <- c(nu[2*m-1],nu[2*m])
		nu <- nu.0
				} else{
	nu0 <- nu
	nu <- optim(unlist(nu),function(x){
		x <- cbind(x,rep(1:(p*M),each=2)) 
		sp <- split(x[,1], x[,2])
		for (m in 1:M)
		for (j in 1:p)
		nu0[[m]][[j]] <- sp[[m+j-1]]
		-sum(log(SMSCN.mix(y,M,w, xi, S, la, nu0))) 
		},method="L-BFGS-B",lower=.01,upper=.99)$par
			nu <- cbind(nu,rep(1:(p*M),each=2)) 
			sp <- split(nu[,1], nu[,2])
		for (m in 1:M)
		for (j in 1:p)
		nu0[[m]][[j]] <- sp[[m+j-1]]
	nu <- nu0
	}
	LL.new <- sum(log(SMSCN.mix(y,M, w, xi, S, la, nu))) # log-likelihood function
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	LL <- LL.new
	}
	if(eqnu){
	aic <- -2 * LL.new + 2 * (M*(3*p+p*(p+1)/2)+M-1)
	bic <- -2 * LL.new + log(n) * (M*(3*p+p*(p+1)/2)+M-1)
	} else {
	aic <- -2 * LL.new + 2 * (M*(4*p+p*(p+1)/2)+M-1)
	bic <- -2 * LL.new + log(n) * (M*(4*p+p*(p+1)/2)+M-1)}
	end <- proc.time()[3]
	time <- end-begin
	obj.out <- list(w=w, xi=xi, S=S, la=la , nu=nu, loglik=LL.new, AIC=aic, BIC=bic, iter=count,elapsed=as.numeric(time),group = apply(z.hat, 
                1, which.max))
	if (group==FALSE)
	obj.out <- obj.out[-length(obj.out)]
	obj.out
	}
	if (family=='SMSTN.mix')
	fit <- EM.SMSTN.mix(y, M=M, w=w, xi=xi, S=S, la=la , nu=nu, iter.max=iter.max, tol=tol, CML=CML, group=group )
	if (family=='SMSSLN.mix')
	fit <- EM.SMSSLN.mix(y, M=M, w=w, xi=xi, S=S, la=la , nu=nu, iter.max=iter.max, tol=tol, CML=CML, group=group )
	if (family=='SMSCN.mix')
	fit <- EM.SMSTN.mix(y, M=M, w=w, xi=xi, S=S, la=la , nu=nu, iter.max=iter.max, tol=tol, group=group )
	return(list(family=family,fit=fit))
	}




