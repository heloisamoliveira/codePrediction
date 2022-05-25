library(dispersion)
library(latex2exp)

K <- 3
N <- 26

load(file="AP.RData")
ap <- Ftreat$XBEST
Ftreat$criterion.opt
(dfap <- Ftreat$criterion.opt[13])
AP <- variance.dispersion(N=26,K=3,ap,REGION=2, SEARCH=3)

# D 
load(file="D.RData")
d <- Ftreat$XBEST
(dfd <- Ftreat$criterion.opt[13])
D <- variance.dispersion(N=26,K=3,d,REGION=2, SEARCH=3)

# DP
load(file="DP.RData")
dp <- Ftreat$XBEST
Ftreat$criterion.opt
(dfdp <- Ftreat$criterion.opt[13])
DP <- variance.dispersion(N=26,K=3,dp,REGION=2, SEARCH=3)

# I 
load(file="I.RData")
i <- Ftreat$XBEST
Ftreat$criterion.opt
(dfi <- Ftreat$criterion.opt[13])
I <- variance.dispersion(N=26,K=3,i,REGION=2, SEARCH=3)

# IP
load(file="IP.RData")
ip <- Ftreat$XBEST
Ftreat$criterion.opt
(dfip <- Ftreat$criterion.opt[13])
IP <- variance.dispersion(N=26,K=3,ip,REGION=2, SEARCH=3)

# I_D 
load(file="ID.RData")
id <- Ftreat$XBEST
(dfid <- Ftreat$criterion.opt[13])
ID <- variance.dispersion(N=26,K=3,id,REGION=2, SEARCH=3)


# I_D_Paper (Borroti, Sambo, Mylona, Gilmour, 2016)
#load("C:\\Users\\Heloisa\\Documents\\Luzia_paper_dic\\techometrics_programas\\mand\\ID_paper.RData")
id_paper <- matrix(scan(,n=K*N),nc=3,byrow=T)
0	-1	1
1	1	-1
0	0	0
-1	-1	1
-1	1	1
1	1	1
-1	-1	1
0	1	0
-1	1	-1
1	-1	1
-1	1	-1
0	0	1
1	0	-1
-1	0	0
0	1	-1
1	-1	-1
-1	1	1
1	1	1
-1	-1	-1
1	0	1
1	-1	0
-1	-1	0
-1	0	-1
0	-1	-1
0	0	0
1	1	0

(dfid_paper <- 5)
ID_paper <- variance.dispersion(N=26,K=3,id_paper,REGION=2, SEARCH=3)


# I_DP
load(file="IDP.RData")
id_p <- Ftreat$XBEST
Ftreat$criterion.opt
(dfid_p <- Ftreat$criterion.opt[13])
ID_P <- variance.dispersion(N=26,K=3,id_p,REGION=2, SEARCH=3)


# Compound D50I50 
load(file="D5+I.RData")
d_i <- Ftreat$XBEST
Ftreat$criterion.opt
(dfd_i <- Ftreat$criterion.opt[13])
D_I <- variance.dispersion(N=26,K=3,d_i,REGION=2, SEARCH=3)

# DP50I50
load(file="DP5+I.RData")
dp_i <- Ftreat$XBEST
Ftreat$criterion.opt
(dfdp_i <- Ftreat$criterion.opt[13])
DP_I <- variance.dispersion(N=26,K=3,dp_i,REGION=2, SEARCH=3)

load("DP5+ID.RData")
dp_id <- Ftreat$XBEST
Ftreat$criterion.opt
(dfdp_id <- Ftreat$criterion.opt[13])
DP_ID <- variance.dispersion(N=26,K=3,dp_id,REGION=2, SEARCH=3)

BBD <- matrix(scan(,n=K*N),nc=K,byrow=TRUE)
-1 -1 0
-1 -1 0
1 -1 0
1 -1 0
-1 1 0
-1 1 0
1 1 0
1 1 0
-1 0 -1
-1 0 -1
1 0 -1
1 0 -1
-1 0 1
-1 0 1
1 0 1
1 0 1
0 -1 -1
0 -1 -1
0 1 -1
0 1 -1
0 -1 1
0 -1 1
0 1 1
0 1 1
0 0 0
0 0 0


bbd <- variance.dispersion(N=26,K=3,BBD,REGION=2,FRACTION=sqrt(2/3),SEARCH=3)
Efficiency(BBD)

ccd <- matrix(scan(,n=K*N),nc=K,byrow=TRUE)
-1 -1 -1
-1 -1  1
-1  1 -1
-1  1  1
1 -1 -1
1 -1  1
1  1 -1
1  1  1
-1 -1 -1
-1 -1  1
-1  1 -1
-1  1  1
1 -1 -1
1 -1  1
1  1 -1
1  1  1
-1  0  0
1  0  0
0 -1  0
0  1  0
0  0 -1
0  0  1
0  0  0
0  0  0
0  0  0
0  0  0


x11()
CCD <- variance.dispersion(N=26,K=3,ccd,REGION=2,SEARCH=3)
Efficiency(ccd)

# TO SEE IF df__ = zero  PE
dfd 
dfdp
dfi 
dfip
dfid
dfid_paper
dfid_p
dfdp_i
(dfbbd <- 13)
(dfccd <- 11)


##############################################
##############################################
# Functions for relative volume for the cube
# Valid for K = 3 only

region_1 <- function(rho, d)
{
  x <- 0.5*(d^2*sqrt(rho^2-2*d^2))
  return(x)
}

region_1_2_theta <- function(rho, d)
{
  x <- acos(d/sqrt(rho^2 - d^2))
  return(x) 
}

region_2_integrand <- function(theta, rho, d)
{ 
  x <- sqrt((cos(theta))^2-(d/rho)^2)/(cos(theta))^3
  return(x)
}

region_2 <- function(rho, d)
{
  i4 <- (d^3/6)*(rho^2/d^2 - 1) * (pi/4 - region_1_2_theta(rho, d))
  i3 <- (d^2*rho/3)*integrate(region_2_integrand, region_1_2_theta(rho,d), pi/4,rho,d)$value
  x <- i4 + i3
  return(x)
}

#spherical region
region_3_integrand <- function(theta, rho, d)
{
  x <- sqrt((cos(theta))^2 - (d/rho)^2)/cos(theta)
  return(x)
}

region_3 <- function(rho, d)
{
  x <- (rho^3/3)*((d/rho)*(pi/4-region_1_2_theta(rho, d))-
                    integrate(region_3_integrand, region_1_2_theta(rho, d), pi/4, rho, d)$value)
  return(x)
}

calc_volume <- function(rho, d)
{
  alpha <- rho/d
  if(alpha<=1)
    vol <- (4/3)*pi*(rho)^3
  if(alpha<=sqrt(2)&alpha>1)
    vol <- (4/3)*pi*(rho)^3-6*pi*(2*(rho)^3 - 3*d*(rho)^2 + d^3)/3
  if(alpha<sqrt(3)&alpha>sqrt(2))
    vol <- 16*(region_1(rho,d) + region_2(rho, d) + region_3(rho,d))
  if(alpha>=sqrt(3))
    vol <- 8*d^3
  return(vol)
}

rho <- seq(0.05, 1, 0.05)
vol <- calc_volume(rho[1], sqrt(1/3))
for(i in 2:20)
{
  x <- calc_volume(rho[i], sqrt(1/3))
  vol <- c(vol,x)
}

vr <- c(0,vol)/vol[20]

# VDGpe

X11()
par(mfrow=c(1,1),bty="n")
plot((D$VDG[,2]*qf(.95,1,dfd))~D$VDG[,1],ylim=c(0.5,3.5),type="n",xlab="Distance",ylab=" " )
#mtext(expression(sigma^-1*(Var(hat(y)(x))*F[1][","][d][";" ][.95])),side=2,line=2,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}$'),cex=.8,side=2,line=2,lwd=1)

lines((I$VDG[,2]*qf(.95,1,dfi))~I$VDG[,1],col=4,lty=1,lwd=2)
#points((I$VDG[,2]*qf(.95,1,dfi))~I$VDG[,1],col=4,pch=8,cex=.5)
lines((I$VDG[,3]*qf(.95,1,dfi))~I$VDG[,1],col=4,lty=1,lwd=2)
#points((I$VDG[,3]*qf(.95,1,dfi))~I$VDG[,1],col=4,pch=8,cex=.5)
lines((IP$VDG[,2]*qf(.95,1,dfip))~IP$VDG[,1],col=1,lty=1,lwd=2)
#points((IP$VDG[,2]*qf(.95,1,dfip))~IP$VDG[,1],col=1,pch=19,cex=.8)
lines((IP$VDG[,3]*qf(.95,1,dfip))~IP$VDG[,1],col=1,lty=1,lwd=2)
#points((IP$VDG[,3]*qf(.95,1,dfip))~IP$VDG[,1],col=1,pch=19,cex=.8)
lines((ID$VDG[,2]*qf(.95,1,dfid))~ID$VDG[,1],col=2,lty=1,lwd=2)
#points((ID$VDG[,2]*qf(.95,1,dfid))~ID$VDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((ID$VDG[,3]*qf(.95,1,dfid))~ID$VDG[,1],col=2,lty=1,lwd=2)
#points((ID$VDG[,3]*qf(.95,1,dfid))~ID$VDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((ID_paper$VDG[,2]*qf(.95,1,dfid_paper))~ID_paper$VDG[,1],col=1,lty=2,lwd=2)
lines((ID_paper$VDG[,3]*qf(.95,1,dfid_paper))~ID_paper$VDG[,1],col=1,lty=2,lwd=2)
lines((ID_P$VDG[,2]*qf(.95,1,dfid_p))~ID_P$VDG[,1],col=6,lty=1,lwd=2)
#points((ID_P$VDG[,2]*qf(.95,1,dfid_p))~ID_P$VDG[,1],col=6,pch=8,cex=.5)
lines((ID_P$VDG[,3]*qf(.95,1,dfid_p))~ID_P$VDG[,1],col=6,lty=1,lwd=2)
#points((ID_P$VDG[,3]*qf(.95,1,dfid_p))~ID_P$VDG[,1],col=6,pch=8,cex=.5)
lines((DP_ID$VDG[,2]*qf(.95,1,dfdp_id))~DP_ID$VDG[,1],col=3,lty=1,lwd=2)
#points((DP_ID$VDG[,2]*qf(.95,1,dfdp_id))~DP_ID$VDG[,1],col=3,pch=24,cex=.8,bg="white")
lines((DP_ID$VDG[,3]*qf(.95,1,dfdp_id))~DP_ID$VDG[,1],col=3,lty=1,lwd=2)
#points((DP_ID$VDG[,3]*qf(.95,1,dfdp_id))~DP_ID$VDG[,1],col=3,pch=24,cex=.8,bg="white")


#legend(0.0,3.5," ",col=2,pch=8,cex=.8,bty="n")
legend(0.0,3.5,"I",lty=1,col=4,bty="n",cex=.8,lwd=2)
#legend(0.25,3.5," ",col=1,pch=19,cex=.8,bty="n")
legend(0.25,3.5,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
#legend(0.50,3.5," ",col=2,pch=22,cex=.8,bg="white",bty="n")
legend(0.50,3.5,expression(I[D]),lty=1,col=2,bty="n",cex=.8,lwd=2)
#legend(0,3.3," ",col=6,pch=8,cex=.8,bty="n")
legend(0,3.3,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
#legend(0.25,3.3," ",col=3,pch=24,cex=.8,bg="white",bty="n")
legend(0.25,3.3,expression(kappa[1]==~kappa[7]==0.5),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0,3.1,expression(I[D]~D[S]~A[S]-sym),lty=2,col=1,bty="n",cex=.8,lwd=2)


######################################################################
# RELATIVE VOLUME

# VDGpe RV

X11()
par(mfrow=c(1,1),bty="n")
plot((D$VDG[,2]*qf(.95,1,dfd))~vr,ylim=c(0.5,3.5),type="n",xlab="Relative Volume",ylab=" " ,
     col=2,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}$'),cex=.8,side=2,line=2,lwd=1)

lines((I$VDG[,2]*qf(.95,1,dfi))~vr,col=4,lty=1,lwd=2)
#points((I$VDG[,2]*qf(.95,1,dfi))~vr,col=1,pch=8,cex=.8)
lines((I$VDG[,3]*qf(.95,1,dfi))~vr,col=4,lty=1,lwd=2)
#points((I$VDG[,3]*qf(.95,1,dfi))~vr,col=1,pch=8,cex=.8)
lines((IP$VDG[,2]*qf(.95,1,dfip))~vr,col=1,lty=1,lwd=2)
#points((IP$VDG[,2]*qf(.95,1,dfip))~vr,col=1,pch=19,cex=.8)
lines((IP$VDG[,3]*qf(.95,1,dfip))~vr,col=1,lty=1,lwd=2)
#points((IP$VDG[,3]*qf(.95,1,dfip))~vr,col=1,pch=19,cex=.8)
lines((ID$VDG[,2]*qf(.95,1,dfid))~vr,col=2,lty=1,lwd=2)
#points((ID$VDG[,2]*qf(.95,1,dfid))~vr,col=2,pch=22,cex=.8,bg="white")
lines((ID$VDG[,3]*qf(.95,1,dfid))~vr,col=2,lty=1,lwd=2)
#points((ID$VDG[,3]*qf(.95,1,dfid))~vr,col=2,pch=22,cex=.8,bg="white")
lines((ID_paper$VDG[,2]*qf(.95,1,dfid_paper))~vr,col=1,lty=2,lwd=2)
lines((ID_paper$VDG[,3]*qf(.95,1,dfid_paper))~vr,col=1,lty=2,lwd=2)
lines((ID_P$VDG[,2]*qf(.95,1,dfid_p))~vr,col=6,lty=1,lwd=2)
#points((ID_P$VDG[,2]*qf(.95,1,dfid_p))~vr,col=6,pch=8,cex=1)
lines((ID_P$VDG[,3]*qf(.95,1,dfid_p))~vr,col=6,lty=1,lwd=2)
#points((ID_P$VDG[,3]*qf(.95,1,dfid_p))~vr,col=6,pch=8,cex=1)
lines((DP_ID$VDG[,2]*qf(.95,1,dfdp_id))~vr,col=3,lty=1,lwd=2)
#points((DP_ID$VDG[,2]*qf(.95,1,dfdp_id))~vr,col=3,pch=24,cex=.8,bg="white")
lines((DP_ID$VDG[,3]*qf(.95,1,dfdp_id))~vr,col=3,lty=1,lwd=2)
#points((DP_ID$VDG[,3]*qf(.95,1,dfdp_id))~vr,col=3,pch=24,cex=.8,bg="white")


legend(0.0,3.5,"I",lty=1,col=4,bty="n",cex=.8,lwd=2)
#legend(0.25,3.5," ",col=1,pch=19,cex=.8,bty="n")
legend(0.25,3.5,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
#legend(0.50,3.5," ",col=2,pch=22,cex=.8,bg="white",bty="n")
legend(0.50,3.5,expression(I[D]),lty=1,col=2,bty="n",cex=.8,lwd=2)
#legend(0,3.3," ",col=6,pch=8,cex=.8,bty="n")
legend(0,3.3,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
#legend(0.25,3.3," ",col=3,pch=24,cex=.8,bg="white",bty="n")
legend(0.25,3.3,expression(kappa[1]==~kappa[7]==0.5),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0,3.1,expression(I[D]~D[S]~A[S]-sym),lty=2,col=1,bty="n",cex=.8,lwd=2)




############################################################################


################
# DVDG
X11()
lim <- (c(min(D$DVDG[,3],DP$DVDG[,3],I$DVDG[,3],IP$DVDG[,3],ID$DVDG[,3],ID_paper$DVDG[,3],
              ID_P$DVDG[,3],DP_ID$DVDG[,3],bbd$DVDG[,2]),
          max(D$DVDG[,2],DP$DVDG[,2],I$DVDG[,2],IP$DVDG[,2],ID$DVDG[,2],ID_paper$DVDG[,2],
              ID_P$DVDG[,2],DP_ID$DVDG[,2],bbd$DVDG[,2])))
par(mfrow=c(1,1),bty="n")
plot((D$DVDG[,2])~D$DVDG[,1],ylim=c(0,1),type="n",xlab="Distance",ylab=" " ,col=2,lwd=2)

mtext(TeX('$\\sigma^{-2}\\{\\textit{Var} \\left(\\hat{y}(\\mathbf{x})-\\hat{y}(\\mathbf{x}_{0})\\right)}$'),cex=.8,side=2,line=2,lwd=1)

lines((I$DVDG[,2])~I$DVDG[,1],col=4,lty=1,lwd=2)
#points((I$DVDG[,2])~I$DVDG[,1],col=4,pch=8,cex=.8)
lines((I$DVDG[,3])~I$DVDG[,1],col=4,lty=1,lwd=2)
#points((I$DVDG[,3])~I$DVDG[,1],col=4,pch=8,cex=.8)
lines((IP$DVDG[,2])~IP$DVDG[,1],col=1,lty=1,lwd=2)
#points((IP$DVDG[,2])~IP$DVDG[,1],col=1,pch=19,cex=.8)
lines((IP$DVDG[,3])~IP$DVDG[,1],col=1,lty=1,lwd=2)
#points((IP$DVDG[,3])~IP$DVDG[,1],col=1,pch=19,cex=.8)
lines((ID_paper$DVDG[,2])~ID_paper$DVDG[,1],col=1,lty=2,lwd=2)
lines((ID_paper$DVDG[,3])~ID_paper$DVDG[,1],col=1,lty=2,lwd=2)
lines((ID$DVDG[,2])~ID$DVDG[,1],col=2,lty=1,lwd=2)
#points((ID$DVDG[,2])~D$DVDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((ID$DVDG[,3])~ID$DVDG[,1],col=2,lty=1,lwd=2)
#points((ID$DVDG[,3])~D$DVDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((ID_P$DVDG[,2])~ID_P$DVDG[,1],col=6,lty=1,lwd=2)
#points((ID_P$DVDG[,2])~ID_P$DVDG[,1],col=6,pch=8,cex=1)
lines((ID_P$DVDG[,3])~ID_P$DVDG[,1],col=6,lty=1,lwd=2)
#points((ID_P$DVDG[,3])~ID_P$DVDG[,1],col=6,pch=8,cex=1)
lines((DP_ID$DVDG[,2])~DP_ID$DVDG[,1],col=3,lty=1,lwd=2)
#points((DP_ID$DVDG[,2])~DP_ID$DVDG[,1],col=3,pch=24,cex=.8,bg="white")
lines((DP_ID$DVDG[,3])~DP_ID$DVDG[,1],col=3,lty=1,lwd=2)


legend(0,1,"I",lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0.25, 1,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0.5,1,expression(I[D]),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(0.72,1,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)

legend(0,.95,expression(kappa[1]==~kappa[7]==0.5),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0.25,.95,expression(I[D]~D[S]~A[S]-sym),lty=2,col=1,bty="n",cex=.8,lwd=2)

# DVDG pe
par(mfrow=c(1,1))
X11()
plot((D$DVDG[,2]*qf(.95,1,dfd))~D$DVDG[,1],ylim=c(0,6),type="l",xlab="Distance",ylab=" " ,
     col=2,lwd=2)
#mtext(expression(sigma^-1*(Var(hat(y)-hat(y)[0])*F[.95][";" ][1][";"][dfPE])),side=2,line=2,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})-\\hat{y}(\\mathbf{x}_{0})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}$'),cex=.8,side=2,line=2,lwd=1)
points((D$DVDG[,2]*qf(.95,1,dfd))~D$DVDG[,1],col=2,pch=21,cex=.8,bg="white")
lines((D$DVDG[,3]*qf(.95,1,dfd))~D$DVDG[,1],col=2,lty=1,lwd=2)
points((D$DVDG[,3]*qf(.95,1,dfd))~D$DVDG[,1],col=2,pch=21,cex=.8,bg="white")
lines((DP$DVDG[,2]*qf(.95,1,dfdp))~DP$DVDG[,1],col=4,lty=3,lwd=2)
lines((DP$DVDG[,3]*qf(.95,1,dfdp))~DP$DVDG[,1],col=4,lty=3,lwd=2)
lines((AP$DVDG[,2]*qf(.95,1,dfap))~AP$DVDG[,1],col=1,lty=1,lwd=2)
lines((AP$DVDG[,3]*qf(.95,1,dfap))~AP$DVDG[,1],col=1,lty=1,lwd=2)
lines((I$DVDG[,2]*qf(.95,1,dfi))~I$DVDG[,1],col=1,lty=1,lwd=2)
points((I$DVDG[,2]*qf(.95,1,dfi))~I$DVDG[,1],col=1,pch=8,cex=.8)
lines((I$DVDG[,3]*qf(.95,1,dfi))~I$DVDG[,1],col=1,lty=1,lwd=2)
points((I$DVDG[,3]*qf(.95,1,dfi))~I$DVDG[,1],col=1,pch=8,cex=.8)
lines((IP$DVDG[,2]*qf(.95,1,dfip))~IP$DVDG[,1],col=1,lty=1,lwd=2)
points((IP$DVDG[,2]*qf(.95,1,dfip))~IP$DVDG[,1],col=1,pch=19,cex=.8)
lines((IP$DVDG[,3]*qf(.95,1,dfip))~IP$DVDG[,1],col=1,lty=1,lwd=2)
points((IP$DVDG[,3]*qf(.95,1,dfip))~IP$DVDG[,1],col=1,pch=19,cex=.8)
lines((ID$DVDG[,2]*qf(.95,1,dfid))~ID$DVDG[,1],col=2,lty=1,lwd=2)
points((ID$DVDG[,2]*qf(.95,1,dfid))~ID$DVDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((ID$DVDG[,3]*qf(.95,1,dfid))~ID$DVDG[,1],col=2,lty=1,lwd=2)
points((ID$DVDG[,3]*qf(.95,1,dfid))~ID$DVDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((ID_paper$DVDG[,2]*qf(.95,1,dfid_paper))~ID_paper$DVDG[,1],col=1,lty=2,lwd=2)
lines((ID_paper$DVDG[,3]*qf(.95,1,dfid_paper))~ID_paper$DVDG[,1],col=1,lty=2,lwd=2)
lines((ID_P$DVDG[,2]*qf(.95,1,dfid_p))~ID_P$DVDG[,1],col=6,lty=1,lwd=2)
points((ID_P$DVDG[,2]*qf(.95,1,dfid_p))~ID_P$DVDG[,1],col=6,pch=8,cex=1)
lines((ID_P$DVDG[,3]*qf(.95,1,dfid_p))~ID_P$DVDG[,1],col=6,lty=1,lwd=2)
points((ID_P$DVDG[,3]*qf(.95,1,dfid_p))~ID_P$DVDG[,1],col=6,pch=8,cex=1)
lines((DP_ID$DVDG[,2]*qf(.95,1,dfdp_id))~DP_ID$DVDG[,1],col=3,lty=1,lwd=2)
points((DP_ID$DVDG[,2]*qf(.95,1,dfdp_id))~DP_ID$DVDG[,1],col=3,pch=24,cex=.8,bg="white")
lines((DP_ID$DVDG[,3]*qf(.95,1,dfdp_id))~DP_ID$DVDG[,1],col=3,lty=1,lwd=2)
points((DP_ID$DVDG[,3]*qf(.95,1,dfdp_id))~DP_ID$DVDG[,1],col=3,pch=24,cex=.8,bg="white")
lines((bbd$DVDG[,2]*qf(.95,1,dfbbd))~bbd$DVDG[,1],col=4,lty=4,lwd=2)
lines((bbd$DVDG[,3]*qf(.95,1,dfbbd))~bbd$DVDG[,1],col=4,lty=4,lwd=2)
lines((CCD$DVDG[,2]*qf(.95,1,dfccd))~CCD$DVDG[,1],col=7,lty=1,lwd=2)
lines((CCD$DVDG[,3]*qf(.95,1,dfccd))~CCD$DVDG[,1],col=7,lty=1,lwd=2)

legend(0,6," ",col=2,pch=21,cex=.8,bg="white",bty="n")
legend(0,6,expression(D[S]),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(0.25,6,expression((DP)[S]),lty=3,col=4,bty="n",cex=.8,lwd=2)
legend(0.50,6,expression((AP)[S]),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0.72,6," ",col=1,pch=8,cex=.8,bty="n")
legend(0.72,6,"I",lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0,5.6," ",col=1,pch=19,cex=.8,bty="n")
legend(0,5.6,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0.25,5.6," ",col=2,pch=22,cex=.8,bg="white",bty="n")
legend(0.25,5.6,expression(I[D]),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(0.50,5.6," ",col=6,pch=8,cex=.8,bty="n")
legend(0.50,5.6,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0,5.2," ",col=3,pch=24,cex=.8,bg="white",bty="n")
legend(0,5.2,expression(kappa[1]==~kappa[7]==0.5),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0.25,5.2,expression(I[D]~D[S]~A[S]-sym),lty=2,col=1,bty="n",cex=.8,lwd=2)
legend(0,4.8,expression(BBD),lty=4,col=4,bty="n",cex=.8,lwd=2)
legend(0.25,4.8,expression(CCD),lty=1,col=7,bty="n",cex=.8,lwd=2)

############################################################
# Relative volume
# DVDG
X11()
par(mfrow=c(1,1),bty="n")
lim <- (c(min(D$DVDG[,3],DP$DVDG[,3],I$DVDG[,3],IP$DVDG[,3],ID$DVDG[,3],ID_paper$DVDG[,3],
              ID_P$DVDG[,3],DP_ID$DVDG[,3],bbd$DVDG[,3]),
          max(D$DVDG[,2],DP$DVDG[,2],I$DVDG[,2],IP$DVDG[,2],ID$DVDG[,2],ID_paper$DVDG[,2],
              ID_P$DVDG[,2],DP_ID$DVDG[,2],bbd$DVDG[,2])))



plot((D$DVDG[,2])~vr,ylim=c(0,0.8),type="n",xlab="Relative Volume",ylab=" " ,col=2,lwd=2)
#mtext(expression(sigma^-2*(Var(hat(y)-hat(y)[0]))),side=2,line=2,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})-\\hat{y}(\\mathbf{x}_{0})\\right)$'),cex=.8,side=2,line=2,lwd=1)

lines((I$DVDG[,2])~vr,col=4,lty=1,lwd=2)
lines((I$DVDG[,3])~vr,col=4,lty=1,lwd=2)
lines((IP$DVDG[,2])~vr,col=1,lty=1,lwd=2)
lines((IP$DVDG[,3])~vr,col=1,lty=1,lwd=2)
lines((ID_paper$DVDG[,2])~vr,col=1,lty=2,lwd=2)
lines((ID_paper$DVDG[,3])~vr,col=1,lty=2,lwd=2)
lines((ID$DVDG[,2])~vr,col=2,lty=1,lwd=2)
lines((ID$DVDG[,3])~vr,col=2,lty=1,lwd=2)
lines((ID_P$DVDG[,2])~vr,col=6,lty=1,lwd=2)
lines((ID_P$DVDG[,3])~vr,col=6,lty=1,lwd=2)
lines((DP_ID$DVDG[,2])~vr,col=3,lty=1,lwd=2)
lines((DP_ID$DVDG[,3])~vr,col=3,lty=1,lwd=2)



legend(0,.8,"I",lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0.25,.8,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0.5,.8,expression(I[D]),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(0,0.75,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0.25,.75,expression(kappa[1]==~kappa[7]==0.5),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0,.7,expression(I[D]~D[S]~A[S]-sym),lty=2,col=1,bty="n",cex=.8,lwd=2)


# DVDG pe rv
X11()
par(mfrow=c(1,1),bty="n")
plot((D$DVDG[,2]*qf(.95,1,dfd))~vr,ylim=c(0,5),type="n",xlab="Relative Volume",ylab=" " ,
     col=2,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})-\\hat{y}(\\mathbf{x}_{0})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}$'),cex=.8,side=2,line=2,lwd=1)

lines((I$DVDG[,2]*qf(.95,1,dfi))~vr,col=4,lty=1,lwd=2)
lines((I$DVDG[,3]*qf(.95,1,dfi))~vr,col=4,lty=1,lwd=2)
lines((IP$DVDG[,2]*qf(.95,1,dfip))~vr,col=1,lty=1,lwd=2)
lines((IP$DVDG[,3]*qf(.95,1,dfip))~vr,col=1,lty=1,lwd=2)
lines((ID$DVDG[,2]*qf(.95,1,dfid))~vr,col=2,lty=1,lwd=2)
lines((ID$DVDG[,3]*qf(.95,1,dfid))~vr,col=2,lty=1,lwd=2)
lines((ID_paper$DVDG[,2]*qf(.95,1,dfid_paper))~vr,col=1,lty=2,lwd=2)
lines((ID_paper$DVDG[,3]*qf(.95,1,dfid_paper))~vr,col=1,lty=2,lwd=2)
lines((ID_P$DVDG[,2]*qf(.95,1,dfid_p))~vr,col=6,lty=1,lwd=2)
lines((ID_P$DVDG[,3]*qf(.95,1,dfid_p))~vr,col=6,lty=1,lwd=3)
lines((DP_ID$DVDG[,2]*qf(.95,1,dfdp_id))~vr,col=3,lty=1,lwd=2)
lines((DP_ID$DVDG[,3]*qf(.95,1,dfdp_id))~vr,col=3,lty=1,lwd=1)

legend(0,5,"I",lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0.25,5,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0.5,5,expression(I[D]),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(0,4.6,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0.25,4.6,expression(kappa[1]==~kappa[7]==0.5),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0,4.2,expression(I[D]~D[S]~A[S]-sym),lty=2,col=1,bty="n",cex=.8,lwd=2)

##################################################################

# FDS pe
X11()
par(mfrow=c(1,1),bty="n")
plot((D$FDS[,2]*qf(.95,1,dfd))~D$FDS[,1],ylim=c(0.5,3.5),type="l",
     xlab="Fraction of design space",ylab=" ",lty=1,col=2,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}$'),cex=.8,side=2,line=2,lwd=1)

lines((I$FDS[,2]*qf(.95,1,dfi))~I$FDS[,1],col=4,lty=1,lwd=2)
#points((I$FDS[i,2]*qf(.95,1,dfi))~I$FDS[i,1],col=1,pch=8,cex=.8)
lines((IP$FDS[,2]*qf(.95,1,dfip))~IP$FDS[,1],col=1,lty=1,lwd=2)
#points((IP$FDS[i,2]*qf(.95,1,dfip))~IP$FDS[i,1],col=1,pch=19,cex=.8)
lines((ID$FDS[,2]*qf(.95,1,dfid))~ID$FDS[,1],col=2,lty=1,lwd=2)
#points((ID$FDS[i,2]*qf(.95,1,dfid))~ID$FDS[i,1],col=2,pch=22,cex=.8,bg="white")
lines((ID_paper$FDS[,2]*qf(.95,1,dfid_paper))~ID_paper$FDS[,1],col=1,lty=2,lwd=2)
lines((ID_P$FDS[,2]*qf(.95,1,dfid_p))~ID_P$FDS[,1],col=6,lty=1,lwd=2)
#points((ID_P$FDS[i,2]*qf(.95,1,dfid_p))~ID_P$FDS[i,1],col=6,pch=8,cex=1)
lines((DP_ID$FDS[,2]*qf(.95,1,dfdp_id))~DP_ID$FDS[,1],col=3,lty=1,lwd=2)
#points((DP_ID$FDS[i,2]*qf(.95,1,dfdp_id))~DP_ID$FDS[i,1],col=3,pch=24,cex=.8,bg="white")

legend(0.0,3.5,"I",lty=1,col=4,bty="n",cex=.8,lwd=2)
#legend(0.25,3.5," ",col=1,pch=19,cex=.8,bty="n")
legend(0.25,3.5,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
#legend(0.50,3.5," ",col=2,pch=22,cex=.8,bg="white",bty="n")
legend(0.50,3.5,expression(I[D]),lty=1,col=2,bty="n",cex=.8,lwd=2)
#legend(0,3.3," ",col=6,pch=8,cex=.8,bty="n")
legend(0,3.3,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
#legend(0.25,3.3," ",col=3,pch=24,cex=.8,bg="white",bty="n")
legend(0.25,3.3,expression(kappa[1]==~kappa[7]==0.5),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0,3.1,expression(I[D]~D[S]~A[S]-sym),lty=2,col=1,bty="n",cex=.8,lwd=2)

# DFDS pe
X11()
par(mfrow=c(1,1),bty="n")
plot((D$FDS[,3]*qf(.95,1,dfd))~D$FDS[,1],ylim=c(0,4),type="n",
     xlab="Fraction of design space",ylab=" ",lty=1,col=2,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})-\\hat{y}(\\mathbf{x}_{0})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}$'),cex=.8,side=2,line=2,lwd=1)

lines((I$FDS[,3]*qf(.95,1,dfi))~I$FDS[,1],col=4,lty=1,lwd=2)
lines((IP$FDS[,3]*qf(.95,1,dfip))~IP$FDS[,1],col=1,lty=1,lwd=2)
lines((ID$FDS[,3]*qf(.95,1,dfid))~ID$FDS[,1],col=2,lty=1,lwd=2)
lines((ID_paper$FDS[,3]*qf(.95,1,dfid_paper))~ID_paper$FDS[,1],col=1,lty=2,lwd=2)
lines((ID_P$FDS[,3]*qf(.95,1,dfid_p))~ID_P$FDS[,1],col=6,lty=1,lwd=2)
lines((DP_ID$FDS[,3]*qf(.95,1,dfdp_id))~DP_ID$FDS[,1],col=3,lty=1,lwd=2)

legend(0,4,"I",lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0.25,4,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0.5,4,expression(I[D]),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(0,3.7,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0.25,3.7,expression(kappa[1]==~kappa[7]==0.5),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0,3.4,expression(I[D]~D[S]~A[S]-sym),lty=2,col=1,bty="n",cex=.8,lwd=2)
