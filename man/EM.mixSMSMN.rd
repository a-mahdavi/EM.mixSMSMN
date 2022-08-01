\name{EM.mixSMSMN}
\alias{EM.mixSMSMN}
\title{EM.mixSMSMN}
\usage{
EM.mixSMSMN(y, M, w, xi, S, la, nu, family="SMSTN.mix", get.init = TRUE, iter.max=50,
	 tol=10^-6, eqnu=FALSE, CML=TRUE, group=TRUE)
}
\description{
 Fit the mixures of SMSMN distributions using EM-algorithm
  M: number of components to fit
  w: the vector of components probabilites.
  xi: the list of vector of location parameters.
  S: the list of cov-variance matrices.
  la: the list of vector of shape parameters.
  nu: the list of flatness parameters.
  family: distribution family to be used in fitting ("SMSTN.mix", "SMSSLN.mix",
  "SMSCN.mix").
  get.init: if TRUE, the initial values are generated.
  iter.max: the maximum number of iterations of the EM algorithm. Default= 50.
  tol: the covergence maximum error. Default= 10^-6.
  eqnu: if TRUE the DF parameters assumed to be equal.
  group: if TRUE it returns the group id of each observation
 }
\examples{

 #  Example 1:
 # Simulating 100 samples from one component SMSTN distribution:
 	y <- r.mixSMSMN(n=100, xi=list(c(0,5)), S=list(matrix(c(1,.4,.4,4),2,2)), la=list(c(-2,3)),
 	nu=list(c(5,3)), family="SMSTN.mix")
 # n: the number of random samples
 # EM output with specific initial values: 
 EM.mixSMSMN(y, xi=list(c(0,5)), S=list(matrix(c(1,.4,.4,4),2,2)), la=list(c(-2,3)), nu=list(c(5,3)),
 		family="SMSTN.mix", get.init=FALSE)
 # EM output without specific initial values: 
 EM.mixSMSMN(y, family="SMSTN.mix", get.init=TRUE)

 # Example 2:
 # Simulating 100 samples from mixtures of SMSTN distributions:
	y <-  r.mixSMSMN(n=200,w=c(.7,.3), xi=list(c(0,5),c(-10,15)), S=list(matrix(c(1,.4,.4,4),2,2),diag(2)),
	 la=list(c(-2,3),c(7,-5)) , nu=list(c(3,10),c(15,5)), family="SMSTN.mix" )

 # EM output with specific initial values: 
 EM.mixSMSMN(y, M=2, w=c(.7,.3), xi=list(c(0,5),c(-10,15)), S=list(matrix(c(1,.4,.4,4),2,2),diag(2)),
	 la=list(c(-2,3),c(7,-5)) , nu=list(c(3,10),c(15,5)), family="SMSTN.mix", get.init=FALSE)
  # EM output without specific initial values: 
  EM.mixSMSMN(y, M=2, family="SMSTN.mix")

 # Example 3:
 # Simulating 100 samples from mixtures of SMSCN distributions:
	y <- r.mixSMSMN(n=100,w=c(.7,.3), xi=list(c(0,5),c(-2,3)), S=list(matrix(c(1,.4,.4,4),2,2),diag(2)),
	 la=list(c(-2,3),c(7,-5)) , nu=list(list(c(.3,.1),c(.5,.7)),list(c(.7,.5),c(.1,.4))), family="SMSCN.mix" )
  # EM output without specific initial values: 
 EM.mixSMSMN(y, M=2, family="SMSCN.mix", get.init=TRUE)
   }

