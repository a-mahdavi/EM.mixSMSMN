\name{r.mixSMSMN}
\alias{r.mixSMSMN}
\title{r.mixSMSMN function}
\usage{
r.mixSMSMN(n ,w, xi, S, la, nu=NULL, family="SMSTN.mix")
}
\description{
 Generating random samples from mixtures of SMSMN distributions 
 n: number of samples.
 w: the vector of probability of each component.
 xi: the list of vector of location parameters.
 S: the list of cov-variance matrices.
 la: the list of vector of shape parameters.
 nu: the list of vector of flatness parameters.
 family: distribution family to be used in fitting ("SMSTN.mix", "SMSSLN.mix", "SMSCN.mix").
	}
\examples{
 # Example 1:
 # Simulating 100 samples from one component SMSTN distribution:
 	y <- r.mixSMSMN(n=100, xi=list(c(0,5)), S=list(matrix(c(1,.4,.4,4),2,2)), la=list(c(-2,3)),
 	nu=list(c(5,3)), family="SMSTN.mix")

 # Example 2:
 # Simulating 100 samples from mixtures of SMSTN distributions:
	y <-  r.mixSMSMN(n=200,w=c(.7,.3), xi=list(c(0,5),c(-10,15)), S=list(matrix(c(1,.4,.4,4),2,2),diag(2)),
	 la=list(c(-2,3),c(7,-5)) , nu=list(c(3,10),c(15,5)), family="SMSTN.mix" )

# Example 3:
 # Simulating 100 samples from mixtures of SMSCN distributions:
	y <- r.mixSMSMN(n=100,w=c(.7,.3), xi=list(c(0,5),c(-2,3)), S=list(matrix(c(1,.4,.4,4),2,2),diag(2)),
	 la=list(c(-2,3),c(7,-5)) , nu=list(list(c(.3,.1),c(.5,.7)),list(c(.7,.5),c(.1,.4))), family="SMSCN.mix" )
   }


