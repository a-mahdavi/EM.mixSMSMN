r.mixSMSMN <- function(n , w=1, xi, S, la, nu=NULL, family="SMSTN.mix"){
	if ((family != "SMSTN.mix") && (family != "SMSSLN.mix") && (family != 
        "SMSCN.mix") )  stop(paste("Family", family, "not recognized.\n", sep = " "))
		p <- length(xi[[1]])
	if( p!=ncol(S[[1]]) || p!=length(la[[1]]) )
	 stop("The size of the parametrs vectors are not compatibles.\n")
		r.SMSTN.mix <- function(n, w , xi, S, la, nu){
	r.SMSTN <- function(n , xi, S, la, nu){
	p <- length(xi) ;	S <- as.matrix(S) ; eig.S = eigen(S) 
	D <- eig.S$vec ; A <- diag(eig.S$val); P <- D%*%A^(1/2)
	 y<-matrix(0,n,p)
		for(i in 1:n){
		del.tau <- 0
		for(j in 1:p)
		del.tau[j] <- rgamma(1,nu[j]/2,nu[j]/2)^(-1/2)
		Z <- rnorm(p)
		Z.tau <- Z*del.tau 
		Z0 <- rnorm(1)
		if(Z0 < (t(la)%*%Z.tau)){
		y[i,] <- Z.tau
		y[i,] <- xi+P%*%y[i,]
			}
		else{
		y[i,] <- -Z.tau
		y[i,] <- xi+P%*%y[i,]
			}	}
		return(y)
		}
		M<- length(xi) ; Z <- rmultinom(n,size=1,prob=w); n <- rowSums(Z)
		p <- length(xi) ;  y<-0
		for( m in 1:M)
		y <- rbind(y,r.SMSTN(n[m], xi[[m]], S[[m]], la[[m]], nu[[m]]))
		y <- y[-1,]
		return(y)
		}
		r.SMSSLN.mix <- function(n, w , xi, S, la, nu){
	r.SMSSLN <- function(n , xi, S, la, nu){
	p <- length(xi) ;	S <- as.matrix(S) ; eig.S = eigen(S) 
	D <- eig.S$vec ; A <- diag(eig.S$val); P <- D%*%A^(1/2)
	 y<-matrix(0,n,p)
		for(i in 1:n){
		del.tau <- 0
		for(j in 1:p)
		del.tau[j] <- rbeta(1,nu[j],1)^(-1/2)
		Z <- rnorm(p)
		Z.tau <- Z*del.tau 
		Z0 <- rnorm(1)
		if(Z0 < (t(la)%*%Z.tau)){
		y[i,] <- Z.tau
		y[i,] <- xi+P%*%y[i,]
			}
		else{
		y[i,] <- -Z.tau
		y[i,] <- xi+P%*%y[i,]
			}	}
		return(y)
		}
		M<- length(xi) ; Z <- rmultinom(n,size=1,prob=w); n <- rowSums(Z)
		p <- length(xi) ;  y<-0
		for( m in 1:M)
		y <- rbind(y,r.SMSSLN(n[m], xi[[m]], S[[m]], la[[m]], nu[[m]]))
		y <- y[-1,]
		return(y)
		}
	r.SMSCN.mix <- function(n, w , xi, S, la, nu){
	r.SMSCN <- function(n , xi, S, la, nu){
	p <- length(xi) ;	S <- as.matrix(S) ; eig.S = eigen(S) 
	D <- eig.S$vec ; A <- diag(eig.S$val); P <- D%*%A^(1/2)
	 y<-matrix(0,n,p)
		for(i in 1:n){
		del.tau <- 0
		for(j in 1:p){
		nuj1 <- as.vector(nu[[j]])[1] ; nuj2 <- as.vector(nu[[j]])[2]
		del.tau[j] <- sample(c(nuj2,1),size=1,prob=c(nuj1,(1-nuj1)))^(-1/2) 
		}
		Z <- rnorm(p)
		Z.tau <- Z*del.tau 
		Z0 <- rnorm(1)
		if(Z0 < (t(la)%*%Z.tau)){
		y[i,] <- Z.tau
		y[i,] <- xi+P%*%y[i,]
			}
		else{
		y[i,] <- -Z.tau
		y[i,] <- xi+P%*%y[i,]
			}	}
		return(y)
		}
		M<- length(xi) ; Z <- rmultinom(n,size=1,prob=w); n <- rowSums(Z)
		p <- length(xi) ;  y<-0
		for( m in 1:M)
		y <- rbind(y,r.SMSCN(n[m], xi[[m]], S[[m]], la[[m]], nu[[m]]))
		y <- y[-1,]
		return(y)
		}
	if (family=='SMSTN.mix')
	y <- r.SMSTN.mix(n, w=w, xi=xi, S=S, la=la , nu=nu )
	if (family=='SMSSLN.mix')
	y <- r.SMSSLN.mix(n, w=w, xi=xi, S=S, la=la , nu=nu )
	if (family=='SMSCN.mix')
	y <- r.SMSCN.mix(n, w=w, xi=xi, S=S, la=la , nu=nu  )
	return(y)
	}









