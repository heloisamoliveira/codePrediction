#################################################################################################
######### PROGRAM FOR CALCULATING EFFICIENCIES AND GRAPHS  ###################
#################################################################################################

#setwd("C:/Users/Luzia/Desktop/DelineamentosEx2")

source('C:/Users/Luzia/Desktop/DelineamentosEx2/function_unblock_Tech.R', encoding = 'UTF-8')


# NUMBER OF FACTORS
K <- 5

# LEVELS OF EACH FACTOR                          
Levels <- list(1:3,1:3,1:3,1:3,1:3) 

# DESIGN REGION: FOR CUBIC ENTER 'Y'  or 'N' FOR SPHERICAL     
Cubic <- 'N'

# NUMBER OF RUNS                       
N <- 30

# COMPOUND CRITERIA
Kap0 <- 1                   # weight for D 
Kap1 <- 1                     # weight for D_PE
Kap2 <- 1                          # weight for CIs (A_PE) - Lpe
Kap3 <- 1                       # weight for point estimation (A) - L
Kap4 <- 1                        # weight for df efficiency (lof)
Kap5 <- 1             # weight for I-otimality  - criteI
Kap6 <- 1                       # crit. I with pure error - critPEI
Kap7 <- 1                         # crit. I - difference of prediction - criteId
Kap8 <- 1                         # crit. I with pure error for difference of prediction - critePEId
Kap9 <- 1                         # crit. I with pure error for difference of prediction - critePEId

Kappa <- c(Kap0,Kap1,Kap2,Kap3,Kap4,Kap5,Kap6,Kap7,Kap8,Kap9)
# INDICATORS OF TERMS IN THE SECOND ORDER MODEL
# ENTRIES IN THE ORDER: LINEAR, QUADRATIC, INTERACTIONS  
Terms <- matrix(c(rep(1,K),rep(1,K),rep(1,(K*(K-1)/2))),nr=1)
Npar <- sum(Terms)

# CONFIDENCE COEFFICIENT FOR DP
# (1-alpha1)
prob1 <- 0.95  


# CONFIDENCE COEFFICIENT FOR AP
# (1-alpha2)
prob2 <- 0.95

# CONFIDENCE COEFFICIENT FOR crit.I / critEPI / criteId / criteEPId
# (1-alpha3)
prob3 <- 0.95


# BONFERRONI'S CORRECTION OF ALPHA FOR AP CRITERION?                                  
MC <- 'N'      # 'Y'=yes or 'N'=no  
prob2 <- ifelse(MC=='Y',prob2^(1/Npar),prob2)

# WEIGHTS FOR PARAMETERS (W DIAGONAL)
# LEAVE AS IT IS FOR WEIGHTS SUGGESTED IN THE PAPER
# NOT USED FOR D, D_PE AND LOF CRITERIA
W <- matrix(rep(1,Npar),nr=1)   
DIV <- ifelse(Cubic=='N',1,4)
for(i in (K+1):(2*K))
{if(Terms[i]==1)
  W[sum(Terms[1:i])] <- W[sum(Terms[1:i])]/DIV
}
W <- matrix(W/sum(W),nr=1)


Terms2 <- matrix(c(rep(1,K),rep(1,K*(K-1)),rep(1,choose(K,3))),nr=1)


load("C:/Users/Luzia/Desktop/DelineamentosEx2/D.RData")
dfd <- Ftreat$criterion.opt[12]
d <- Ftreat$XBest_ord
critD <- Efficiency(d)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/DP.RData")
dp <- Ftreat$XBest_ord
critDP <- Efficiency(dp)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/A.RData")
a <- Ftreat$XBest_ord
critA <- Efficiency(a)

# AP
D50_IDP <- matrix(scan(,n=K*30),nc=K,byrow=TRUE)
-1 -1 -1 -1 -1
-1 -1 -1 1 1
-1 -1 1 -1 1
-1 -1 1 -1 1
-1 -1 1 1 -1
-1 1 -1 -1 1
-1 1 -1 1 -1
-1 1 1 -1 -1
-1 1 1 1 1
0 -2.236068 0 0 0
0 0 0 -2.236068 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 2.236068
0 0 2.236068 0 0
1 -1 -1 -1 1
1 -1 -1 1 -1
1 -1 -1 1 -1
1 -1 1 -1 -1
1 -1 1 1 1
1 1 -1 1 1
1 1 -1 1 1
1 1 1 -1 1
1 1 1 1 -1
1 1 1 1 -1
1.118034 1.118034 -1.118034 0 -1.118034
1.118034 1.118034 -1.118034 0 -1.118034
2.236068 0 0 0 0

ap <- D50_IDP
critAP <- Efficiency(ap)

APpaper <- matrix(scan(,n=150),ncol=5,byrow=T)
-1	-1	-1	1	-1
-1	-1	-1	-1	1
-1	-1	1	1	1
-1	-1	1	1	1
-1	-1	1	-1	-1
-1	1	-1	-1	-1
-1	1	-1	1	1
-1	1	-1	1	1
-1	1	1	-1	1
-1	1	1	1	-1
-1	1	1	1	-1
1	-1	-1	-1	-1
1	-1	-1	1	1
1	-1	1	-1	1
1	-1	1	1	-1
1	-1	1	-1	1
1	1	-1	1	-1
1	1	-1	-1	1
1	1	1	1	1
1	1	1	-1	-1
2.24	0	0	0	0
0	-2.24	0	0	0
0	0	-2.24	0	0
0	0	0	-2.24	0
0	0	0	2.24	0
0	0	0	0	-2.24
0	0	0	0	0
0	0	0	0	0
0	0	0	0	0
0	0	0	0	0



AP2 <- matrix(scan(,n=K*30),nc=K,byrow=TRUE)
-1.290994 0 1.290994 -1.290994 0
-1.290994 0 1.290994 -1.290994 0
-1.118034 -1.118034 -1.118034 0 -1.118034
-1.118034 -1.118034 0 1.118034 1.118034
-1.118034 1.118034 -1.118034 0 1.118034
-1.118034 1.118034 0 1.118034 -1.118034
0 -1.581139 0 -1.581139 0
0 -1.581139 0 -1.581139 0
0 -1.118034 1.118034 1.118034 -1.118034
0 0 -1.581139 -1.581139 0
0 0 -1.581139 1.581139 0
0 0 -1.581139 1.581139 0
0 0 0 -1.581139 -1.581139
0 0 0 -1.581139 -1.581139
0 0 0 -1.581139 1.581139
0 0 0 -1.581139 1.581139
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 1.118034 1.118034 1.118034 1.118034
0 1.581139 0 -1.581139 0
0 1.581139 0 -1.581139 0
1.118034 -1.118034 -1.118034 0 -1.118034
1.118034 -1.118034 -1.118034 0 1.118034
1.118034 -1.118034 1.118034 0 1.118034
1.118034 1.118034 -1.118034 0 -1.118034
1.118034 1.118034 -1.118034 0 1.118034
1.118034 1.118034 1.118034 0 -1.118034
1.581139 0 0 -1.581139 0
1.581139 0 0 1.581139 0

critAP2 <- Efficiency(AP2)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/APnaoOtimo.RData")
critAP0 <- Efficiency(Ftreat$XBest_ord)


# Not the best
load("C:/Users/Luzia/Desktop/DelineamentosEx2/I.RData")
i <- d
critI <- Efficiency(i)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/IP.RData")
ip <- Ftreat$XBest_ord
critIP <- Efficiency(ip)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/ID.RData")
id <- Ftreat$XBest_ord
critID <- Efficiency(id)


load("C:/Users/Luzia/Desktop/DelineamentosEx2/IDPnaoOtimo.RData")
critIDP0 <- Efficiency(Ftreat$XBest_ord)

# IDP
D30_IDP <- matrix(scan(,n=K*30),nc=K,byrow=TRUE)
-1 -1 -1 -1 1
-1 -1 -1 1 -1
-1 -1 -1 1 -1
-1 -1 1 -1 -1
-1 -1 1 1 1
-1 1 -1 -1 -1
-1 1 -1 1 1
-1 1 1 -1 1
-1 1 1 1 -1
-1 1 1 1 -1
0 -2.236068 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 0
0 0 0 0 2.236068
0 0 0 2.236068 0
0 0 2.236068 0 0
0 1.118034 -1.118034 1.118034 -1.118034
0 1.118034 -1.118034 1.118034 -1.118034
1 -1 -1 -1 -1
1 -1 -1 1 1
1 -1 1 -1 1
1 -1 1 -1 1
1 -1 1 1 -1
1 1 -1 -1 1
1 1 1 -1 -1
1 1 1 1 1
2.236068 0 0 0 0

idp <- D30_IDP
critIDP <- Efficiency(idp)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/DP30+ID.RData")
K1_30K7_70 <- Ftreat$XBest_ord
critK1_30K7_70 <- Efficiency(K1_30K7_70)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/DP10+ID90.RData")
K1_10K7_90 <- Ftreat$XBest_ord
critK1_10K7_90 <- Efficiency(K1_10K7_90)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/D90+IDP10.RData")
K0_90K8_10 <- Ftreat$XBest_ord
critK0_90K8_10 <- Efficiency(K0_90K8_10)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/D2DP2A3ID3_AP+IDPotimo.RData")
K0_2K1_2K2_3K7_3 <- Ftreat$XBest_ord
critK0_2K1_2K2_3K7_3 <- Efficiency(K0_2K1_2K2_3K7_3)

load("C:/Users/Luzia/Desktop/DelineamentosEx2/BIAS5F.RData")
critBias <- Efficiency(Bias5F$XBest_ord)

###################
# OTHERS
load("C:/Users/Luzia/Desktop/DelineamentosEx2/D10+IDP90.RData")
critD10IDP90 <- Efficiency(Ftreat$XBest_ord)


load("C:/Users/Luzia/Desktop/DelineamentosEx2/DP50+ID.RData")
critDP50ID <- Efficiency(Ftreat$XBest_ord)
###################

 crit <- rbind(critD,critDP,critA,critAP,critIP,critID,critIDP,critK1_30K7_70,critK1_10K7_90,critK0_90K8_10,
               critK0_2K1_2K2_3K7_3,critBias)
dim(crit)
maxcrit <- apply(crit[,-c(10,11)],2,max)
minS <- min(crit[,10])
maxcrit <- matrix(maxcrit,nr=nrow(crit),nc=9,byrow=T)


ef <- round(cbind(100*crit[,-c(10,11)]/maxcrit,100*minS/crit[,10],crit[,11]),2)
ef



#############################################
library(dispersion)

# Os delineamentos s?o para 5 fatores em regi?o esferica.
# K=5 e o modelo de segunda ordem.

# delineamento D/I
dfd <- 0
D <- variance.dispersion(N=30,K=5,d,REGION=1,ORDER=2,ITMAX= 10000, FTOL=0.000001,SEARCH=3)

# delineamento DP
dfdp <- 9
DP <- variance.dispersion(N=30,K=5,dp,REGION=1,ORDER=2,SEARCH=3)

# delineamento A
dfa <- 1
A <- variance.dispersion(N=30,K=5,a,REGION=1,ORDER=2,SEARCH=3)

# delineamento AP
dfap <- 8
AP <- variance.dispersion(N=30,K=5,ap,REGION=1,NLOOPS=50000,SEARCH=3)

# IP

dfip <- 8
IP <- variance.dispersion(N=30,K=5,ip,REGION=1,ORDER=2,SEARCH=3)

# ID
dfid <- 3
ID <- variance.dispersion(N=30,K=5,id,REGION=1,ORDER=2,SEARCH=3)

# IDP D30_IDP
dfidp <- 8
IDP <- variance.dispersion(N=30,K=5,idp,REGION=1,ORDER=2,SEARCH=3)

# DP30_ID70
dfdp30id <- 7
DP30ID <- variance.dispersion(N=30,K=5,K1_30K7_70,REGION=1,ORDER=2,SEARCH=3)

############################
# Not icluded
load("DP50+ID.RData")
dp50id <- Ftreat$XBEST
dfdp50id <- Ftreat$criterion.opt[12]
DP50ID <- variance.dispersion(N=30,K=5,dp50id,REGION=1,ORDER=2,SEARCH=3)
##############################

# DP10_ID90
dfdp10id <- 5
DP10ID <- variance.dispersion(N=30,K=5,K1_10K7_90,REGION=1,ORDER=2,NLOOPS=50000,SEARCH=3)

# D90_IDP10
dfd90idp <- 5
D90IDP <- variance.dispersion(N=30,K=5,K0_90K8_10,REGION=1,ORDER=2,FTOL=0.000001,ITMAX=50000,SEARCH=3)

df2233 <- 5
D2233 <- variance.dispersion(N=30,K=5,K0_2K1_2K2_3K7_3,REGION=1,ORDER=2,SEARCH=3)

##############################################################
# GRAPHS
vol <- ((pi^(K/2))/gamma(K/2+1))*DP$VDG[,1]^K
RV <- DP$VDG[,1]^K  # Relative volume

###############################################################
# LEGEND
# D: col=2,lty=1
# DP: col=2,lty=2
# A: col=1,lty=1
# AP: col=1,lty=2
# I: SAME D
# IP: col=2,lty=3
# ID: col=4,lty=1
# IDP: col=4,lty=2
# DP30ID: col=3,lty=1
# DP10ID: col=6, lty=1
# D90IDP: col=1,lty=3
# D2233=D2DP2A3ID3: col=3,lty=2

# VDG (ONLY FOR SUPPLEMENTARY MATERIAL)
X11()
par(mfrow=c(1,1))
yvdg=c(min((D$VDG[,3]),(DP$VDG[,3]),(A$VDG[,3]),(AP$VDG[,3]),(IP$VDG[,3]),
           (ID$VDG[,3]),(IDP$VDG[,3]),
           (DP30ID$VDG[,3]),(DP10ID$VDG[,3]),(D90IDP$VDG[,3]),(D2233$VDG[,3])),
       max((DP$VDG[,2]),(AP$VDG[,2]),(IP$VDG[,2]),(ID$VDG[,2]),(IDP$VDG[,2]),
           (DP30ID$VDG[,2]),(DP10ID$VDG[,2]),(D90IDP$VDG[,2]),(D2233$VDG[,2])))

plot((D$VDG[,2])~D$VDG[,1],ylim=yvdg,type="l",xlab="Distance",
     ylab=" ",col=2,lty=1,lwd=2.0)
mtext(TeX('$\\sigma^{-2}\\textit{Var}\\left(\\hat{y}(\\mathbf{x})\\right)}$'),cex=.8,side=2,line=2,lwd=1)

points((D$VDG[,2])~D$VDG[,1],col=2,pch=21,cex=1,bg="white")
lines((D$VDG[,3])~D$VDG[,1],col=2,lty=1,lwd=2)
points((D$VDG[,3])~D$VDG[,1],col=2,pch=21,cex=1,bg="white")
lines((DP$VDG[,2])~DP$VDG[,1],col=4,lty=3,lwd=2)
lines((DP$VDG[,3])~DP$VDG[,1],col=4,lty=3,lwd=2)
lines((AP$VDG[,2])~AP$VDG[,1],col="cyan",lty=1,lwd=2)
points((AP$VDG[,2])~D$VDG[,1],col="cyan",pch=23,cex=1.2,bg="white")
lines((AP$VDG[,3])~AP$VDG[,1],col="cyan",lty=1,lwd=2)
points((AP$VDG[,3])~D$VDG[,1],col="cyan", pch=23,cex=1.2,bg="white")
lines((IP$VDG[,2])~IP$VDG[,1],col=1,lty=1,lwd=2)
points((IP$VDG[,2])~IP$VDG[,1],col=1,pch=19,cex=.8)
lines((IP$VDG[,3])~IP$VDG[,1],col=1,lty=1,lwd=2)
points((IP$VDG[,3])~IP$VDG[,1],col=1,pch=19,cex=.8)
lines((A$VDG[,2])~A$VDG[,1],col=1,lty=2,lwd=2)
lines((A$VDG[,3])~A$VDG[,1],col=1,lty=2,lwd=2)
lines((ID$VDG[,2])~ID$VDG[,1],col=2,lty=1,lwd=2)      #CCD
points((ID$VDG[,2])~ID$VDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((ID$VDG[,3])~ID$VDG[,1],col=2,lty=1,lwd=2)
points((ID$VDG[,3])~ID$VDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((IDP$VDG[,2])~IDP$VDG[,1],col=6,lty=1,lwd=2)
points((IDP$VDG[,2])~IDP$VDG[,1],col=6,pch=8,cex=1)
lines((IDP$VDG[,3])~IDP$VDG[,1],col=6,lty=1,lwd=2)
points((IDP$VDG[,3])~IDP$VDG[,1],col=6,pch=8,cex=1)
lines((DP30ID$VDG[,2])~DP30ID$VDG[,1],col=3,lty=1,lwd=2)
points((DP30ID$VDG[,2])~DP30ID$VDG[,1],col=3,pch=24,cex=.8,bg="white")
lines((DP30ID$VDG[,3])~DP30ID$VDG[,1],col=3,lty=1,lwd=2)
points((DP30ID$VDG[,3])~DP30ID$VDG[,1],col=3,pch=24,cex=.8,bg="white")
lines((DP10ID$VDG[,2])~DP10ID$VDG[,1],col=4,lty=4,lwd=2)
lines((DP10ID$VDG[,3])~DP10ID$VDG[,1],col=4,lty=4,lwd=2)
lines((D90IDP$VDG[,2])~D90IDP$VDG[,1],col=7,lty=1,lwd=2)
lines((D90IDP$VDG[,3])~D90IDP$VDG[,1],col=7,lty=1,lwd=2)
lines((D2233$VDG[,2])~D2233$VDG[,1],col=1,lty=1,lwd=2)
points((D2233$VDG[,2])~D2233$VDG[,1],col=1,pch=8,cex=.8)
lines((D2233$VDG[,3])~D2233$VDG[,1],col=1,lty=1,lwd=2)
points((D2233$VDG[,3])~D2233$VDG[,1],col=1,pch=8,cex=.8)
legend(0,2," ",col=2,pch=21,cex=.8,bg="white",bty="n")
legend(0,2,expression(D[S]/I),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(.2,2,expression((DP)[S]),lty=3,col=4,bty="n",cex=.8,lwd=2)
legend(.4,2,expression(A[S]),lty=2,col=1,bty="n",cex=.8,lwd=2)
legend(.6,2," ",col="cyan",pch=23,cex=.8,bg="white",bty="n")
legend(.6,2,expression((AP)[S]),lty=1,col="cyan",bty="n",cex=.8,lwd=2)
legend(0,1.9," ",col=1,pch=19,cex=.8,bty="n")
legend(0,1.9,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(.2,1.9," ",col=2,pch=22,cex=.8,bg="white",bty="n")
legend(.2,1.9,expression(I[D]/CCD),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(.4,1.9," ",col=6,pch=8,cex=.8,bty="n")
legend(.4,1.9,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0,1.8," ",col=3,pch=24,cex=.8,bg="white",bty="n")
legend(0,1.8,expression(8:kappa[1]==0.3~kappa[7]==0.7),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0.4,1.8,expression(9:kappa[1]==0.1~kappa[7]==0.9),lty=4,col=4,bty="n",cex=.8,lwd=2)
legend(0,1.7,expression(10:kappa[0]==0.9~kappa[8]==0.1),lty=1,col=7,bty="n",cex=.8,lwd=2)
legend(0.4,1.7," ",col=1,pch=8,cex=.8,bty="n")
legend(0.4,1.7,expression(11:kappa[0]==~kappa[1]==0.2~kappa[4]==~kappa[7]==0.3),lty=1,col=1,bty="n",cex=.8,lwd=2)


###############################################################
# VDG PE (TO BE IN THE PAPER)
X11()
par(mfrow=c(1,1),bty="n")
yvdg=c(min((DP$VDG[,3]*qf(.95,1,dfdp)),(AP$VDG[,3]*qf(.95,1,dfap)),
           (IP$VDG[,3]*qf(.95,1,dfip)),(ID$VDG[,3]*qf(.95,1,dfid)),(IDP$VDG[,3]*qf(.95,1,dfidp)),
           (DP30ID$VDG[,3]*qf(.95,1,dfdp30id)),(DP10ID$VDG[,3]*qf(.95,1,dfdp10id)),
           (D90IDP$VDG[,3]*qf(.95,1,dfd90idp)),(D2233$VDG[,3]*qf(.95,1,df2233))),
       max((DP$VDG[,2]*qf(.95,1,dfdp)),(AP$VDG[,2]*qf(.95,1,dfap)),
           (IP$VDG[,2]*qf(.95,1,dfip)),(ID$VDG[,2]*qf(.95,1,dfid)),(IDP$VDG[,2]*qf(.95,1,dfidp)),
           (DP30ID$VDG[,2]*qf(.95,1,dfdp30id)),(DP10ID$VDG[,2]*qf(.95,1,dfdp10id)),
           (D90IDP$VDG[,2]*qf(.95,1,dfd90idp)),(D2233$VDG[,2]*qf(.95,1,df2233))))

plot((DP$VDG[,2]*qf(.95,1,dfdp))~DP$VDG[,1],ylim=c(0,12),type="n",xlab="Distance",
     ylab=" ",col=4,lty=3,lwd=2.0)


lines((IP$VDG[,2]*qf(.95,1,dfip))~IP$VDG[,1],col=1,lty=1,lwd=2)
lines((IP$VDG[,3]*qf(.95,1,dfip))~IP$VDG[,1],col=1,lty=1,lwd=2)
lines((ID$VDG[,2]*qf(.95,1,dfid))~ID$VDG[,1],col=2,lty=1,lwd=2)      #CCD
lines((ID$VDG[,3]*qf(.95,1,dfid))~ID$VDG[,1],col=2,lty=1,lwd=2)
lines((IDP$VDG[,2]*qf(.95,1,dfidp))~IDP$VDG[,1],col=6,lty=1,lwd=2)
lines((IDP$VDG[,3]*qf(.95,1,dfidp))~IDP$VDG[,1],col=6,lty=1,lwd=2)
lines((DP30ID$VDG[,2]*qf(.95,1,dfdp30id))~DP30ID$VDG[,1],col=3,lty=1,lwd=2)
lines((DP30ID$VDG[,3]*qf(.95,1,dfdp30id))~DP30ID$VDG[,1],col=3,lty=1,lwd=2)
lines((DP10ID$VDG[,2]*qf(.95,1,dfdp10id))~DP10ID$VDG[,1],col=4,lty=1,lwd=2)
lines((DP10ID$VDG[,3]*qf(.95,1,dfdp10id))~DP10ID$VDG[,1],col=4,lty=1,lwd=2)
lines((D90IDP$VDG[,2]*qf(.95,1,dfd90idp))~D90IDP$VDG[,1],col=7,lty=1,lwd=2)
lines((D90IDP$VDG[,3]*qf(.95,1,dfd90idp))~D90IDP$VDG[,1],col=7,lty=1,lwd=2)
#lines((D2233$VDG[,2]*qf(.95,1,df2233))~D2233$VDG[,1],col=1,lty=1,lwd=2)
#points((D2233$VDG[,2]*qf(.95,1,df2233))~D2233$VDG[,1],col=1,pch=8,cex=.8)
#lines((D2233$VDG[,3]*qf(.95,1,df2233))~D2233$VDG[,1],col=1,lty=1,lwd=2)
#points((D2233$VDG[,3]*qf(.95,1,df2233))~D2233$VDG[,1],col=1,pch=8,cex=.8)

legend(0,12,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0.2,12,expression(I[D]/CCD),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(0.4,12,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0,11.2,expression(8:kappa[1]==0.3~kappa[7]==0.7),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0.4,11.2,expression(9:kappa[1]==0.1~kappa[7]==0.9),lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0,10.4,expression(10:kappa[0]==0.9~kappa[8]==0.1),lty=1,col=7,bty="n",cex=.8,lwd=2)
#legend(0.4, 10.2,expression(11:kappa[0]==~kappa[1]==0.2~kappa[4]==~kappa[7]==0.3),lty=1,col=1,bty="n",cex=.8,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var}\\left(\\hat{y}(\\mathbf{x})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}}$'),cex=.8,side=2,line=2,lwd=1)

############################################################################

###############################################################################
# VDG pe vs RELATIVE VOLUME (TO BE IN THE PAPER) 
X11()
par(mfrow=c(1,1),bty="n")
yvdg=c(min((DP$VDG[,3]*qf(.95,1,dfdp)),(AP$VDG[,3]*qf(.95,1,dfap)),
           (IP$VDG[,3]*qf(.95,1,dfip)),(ID$VDG[,3]*qf(.95,1,dfid)),(IDP$VDG[,3]*qf(.95,1,dfidp)),
           (DP30ID$VDG[,3]*qf(.95,1,dfdp30id)),(DP10ID$VDG[,3]*qf(.95,1,dfdp10id)),
           (D90IDP$VDG[,3]*qf(.95,1,dfd90idp)),(D2233$VDG[,3]*qf(.95,1,df2233))),
       max((DP$VDG[,2]*qf(.95,1,dfdp)),(AP$VDG[,2]*qf(.95,1,dfap)),
           (IP$VDG[,2]*qf(.95,1,dfip)),(ID$VDG[,2]*qf(.95,1,dfid)),(IDP$VDG[,2]*qf(.95,1,dfidp)),
           (DP30ID$VDG[,2]*qf(.95,1,dfdp30id)),(DP10ID$VDG[,2]*qf(.95,1,dfdp10id)),
           (D90IDP$VDG[,2]*qf(.95,1,dfd90idp)),(D2233$VDG[,2]*qf(.95,1,df2233))))

  plot((DP$VDG[,2]*qf(.95,1,dfdp))~RV,ylim=c(0,12),type="n",xlab="Relative Volume",
       ylab=" ",col=4,lty=3,lwd=2.0)

lines((IP$VDG[,2]*qf(.95,1,dfip))~RV,col=1,lty=1,lwd=2)
lines((IP$VDG[,3]*qf(.95,1,dfip))~RV,col=1,lty=1,lwd=2)
lines((ID$VDG[,2]*qf(.95,1,dfid))~RV,col=2,lty=1,lwd=2)      #CCD
lines((ID$VDG[,3]*qf(.95,1,dfid))~RV,col=2,lty=1,lwd=2)
lines((IDP$VDG[,2]*qf(.95,1,dfidp))~RV,col=6,lty=1,lwd=2)
lines((IDP$VDG[,3]*qf(.95,1,dfidp))~RV,col=6,lty=1,lwd=2)
lines((DP30ID$VDG[,2]*qf(.95,1,dfdp30id))~RV,col=3,lty=1,lwd=2)
lines((DP30ID$VDG[,3]*qf(.95,1,dfdp30id))~RV,col=3,lty=1,lwd=2)
lines((DP10ID$VDG[,2]*qf(.95,1,dfdp10id))~RV,col=4,lty=1,lwd=2)
lines((DP10ID$VDG[,3]*qf(.95,1,dfdp10id))~RV,col=4,lty=1,lwd=2)
lines((D90IDP$VDG[,2]*qf(.95,1,dfd90idp))~RV,col=7,lty=1,lwd=2)
lines((D90IDP$VDG[,3]*qf(.95,1,dfd90idp))~RV,col=7,lty=1,lwd=2)
#lines((D2233$VDG[,2]*qf(.95,1,df2233))~RV,col=1,lty=1,lwd=2)
#points((D2233$VDG[j,2]*qf(.95,1,df2233))~RV[j],col=1,pch=8,cex=.8)
#lines((D2233$VDG[,3]*qf(.95,1,df2233))~RV,col=1,lty=1,lwd=2)
#points((D2233$VDG[j,3]*qf(.95,1,df2233))~RV[j],col=1,pch=8,cex=.8)

legend(0,12,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0.2,12,expression(I[D]/CCD),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(0.4,12,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0,11.2,expression(8:kappa[1]==0.3~kappa[7]==0.7),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0.4,11.2,expression(9:kappa[1]==0.1~kappa[7]==0.9),lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0,10.4,expression(10:kappa[0]==0.9~kappa[8]==0.1),lty=1,col=7,bty="n",cex=.8,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var}\\left(\\hat{y}(\\mathbf{x})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}}$'),cex=.8,side=2,line=2,lwd=1)

##################################################################################################

#######################################################################################################
# FDS PE (TO BE IN THE PAPER)
X11()
par(mfrow=c(1,1))
{
i=c(seq(1,10000,by=500),9800,9950,9970,9980,9990,9995,9997,9998,9999,10000)
yvdg=c(min((DP$FDS[,2]*qf(.95,1,dfdp)),(AP$FDS[,2]*qf(.95,1,dfap)),
           (IP$FDS[,2]*qf(.95,1,dfip)),(ID$FDS[,2]*qf(.95,1,dfid)),(IDP$FDS[,2]*qf(.95,1,dfidp)),
           (DP30ID$FDS[,2]*qf(.95,1,dfdp30id)),(DP10ID$FDS[,2]*qf(.95,1,dfdp10id)),
           (D90IDP$FDS[,2]*qf(.95,1,dfd90idp)),(D2233$FDS[,2]*qf(.95,1,df2233))),
       max((DP$FDS[,2]*qf(.95,1,dfdp)),(AP$FDS[,2]*qf(.95,1,dfap)),
           (IP$FDS[,2]*qf(.95,1,dfip)),(ID$FDS[,2]*qf(.95,1,dfid)),(IDP$FDS[,2]*qf(.95,1,dfidp)),
           (DP30ID$FDS[,2]*qf(.95,1,dfdp30id)),(DP10ID$FDS[,2]*qf(.95,1,dfdp10id)),
           (D90IDP$FDS[,2]*qf(.95,1,dfd90idp)),(D2233$FDS[,2]*qf(.95,1,df2233))))

plot((DP$FDS[,2]*qf(.95,1,dfdp))~DP$FDS[,1],ylim=c(0,10),type="n",xlab="Fraction of design space",
     ylab=" ",col=4,lty=3,lwd=2.0)
lines((IP$FDS[,2]*qf(.95,1,dfip))~IP$FDS[,1],col=1,lty=1,lwd=2)
lines((ID$FDS[,2]*qf(.95,1,dfid))~ID$FDS[,1],col=2,lty=1,lwd=2)      #CCD
lines((IDP$FDS[,2]*qf(.95,1,dfidp))~IDP$FDS[,1],col=6,lty=1,lwd=2)
lines((DP30ID$FDS[,2]*qf(.95,1,dfdp30id))~DP30ID$FDS[,1],col=3,lty=1,lwd=2)
lines((DP10ID$FDS[,2]*qf(.95,1,dfdp10id))~DP10ID$FDS[,1],col=4,lty=1,lwd=2)
lines((D90IDP$FDS[,2]*qf(.95,1,dfd90idp))~D90IDP$FDS[,1],col=7,lty=1,lwd=2)


legend(0,10,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(0.2,10,expression(I[D]/CCD),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(0.4,10,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0,9.3,expression(8:kappa[1]==0.3~kappa[7]==0.7),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(0.4,9.3,expression(9:kappa[1]==0.1~kappa[7]==0.9),lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0,8.6,expression(10:kappa[0]==0.9~kappa[8]==0.1),lty=1,col=7,bty="n",cex=.8,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var}\\left(\\hat{y}(\\mathbf{x})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}}$'),cex=.8,side=2,line=2,lwd=1)
}

#######################################################################################
# DVDG RELATIVE VOLUME 
X11()
par(mfrow=c(1,1),bty="n")
{
yvdg=c(min((D$DVDG[,3]),(DP$DVDG[,3]),(A$DVDG[,3]),(AP$DVDG[,3]),(IP$DVDG[,3]),
           (ID$DVDG[,3]),(IDP$DVDG[,3]),
           (DP30ID$DVDG[,3]),(DP10ID$DVDG[,3]),(D90IDP$DVDG[,3]),(D2233$DVDG[,3])),
       max((DP$DVDG[,2]),(AP$DVDG[,2]),(IP$DVDG[,2]),(ID$DVDG[,2]),(IDP$DVDG[,2]),
           (DP30ID$DVDG[,2]),(DP10ID$DVDG[,2]),(D90IDP$DVDG[,2]),(D2233$DVDG[,2])))

plot((D$DVDG[,2])~RV,ylim=c(0,3),type="n",xlab="Relative Volume",
     ylab=" ",col=2,lty=2,lwd=2.0)
lines((IP$DVDG[,2])~RV,col=1,lty=1,lwd=2)
lines((IP$DVDG[,3])~RV,col=1,lty=1,lwd=2)
lines((ID$DVDG[,2])~RV,col=2,lty=1,lwd=2)      #CCD
lines((ID$DVDG[,3])~RV,col=2,lty=1,lwd=2)
lines((IDP$DVDG[,2])~RV,col=6,lty=1,lwd=2)
lines((IDP$DVDG[,3])~RV,col=6,lty=1,lwd=2)
lines((DP30ID$DVDG[,2])~RV,col=3,lty=1,lwd=2)
lines((DP30ID$DVDG[,3])~RV,col=3,lty=1,lwd=2)
lines((DP10ID$DVDG[,2])~RV,col=4,lty=1,lwd=2)
lines((DP10ID$DVDG[,3])~RV,col=4,lty=1,lwd=2)
lines((D90IDP$DVDG[,2])~RV,col=7,lty=1,lwd=2)
lines((D90IDP$DVDG[,3])~RV,col=7,lty=1,lwd=2)


legend(0,3,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(.2,3,expression(I[D]/CCD),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(.4,3,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0,2.83,expression(8:kappa[1]==0.3~kappa[7]==0.7),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(.4,2.83,expression(9:kappa[1]==0.1~kappa[7]==0.9),lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0,2.66,expression(10:kappa[0]==0.9~kappa[8]==0.1),lty=1,col=7,bty="n",cex=.8,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var}\\left(\\hat{y}(\\mathbf{x})-\\hat{y}(\\mathbf{x}_{0})\\right)}$'),cex=.8,side=2,line=2,lwd=1)
}

#######################################################################################
# DVDG pe
X11()
par(mfrow=c(1,1))
{
yvdg=c(min((DP$DVDG[,3]*qf(.95,1,dfdp)),(AP$DVDG[,3]*qf(.95,1,dfap)),
           (IP$DVDG[,3]*qf(.95,1,dfip)),(ID$DVDG[,3]*qf(.95,1,dfid)),(IDP$DVDG[,3]*qf(.95,1,dfidp)),
           (DP30ID$DVDG[,3]*qf(.95,1,dfdp30id)),(DP10ID$DVDG[,3]*qf(.95,1,dfdp10id)),
           (D90IDP$DVDG[,3]*qf(.95,1,dfd90idp)),(D2233$DVDG[,3]*qf(.95,1,df2233))),
       max((DP$DVDG[,2]*qf(.95,1,dfdp)),(AP$DVDG[,2]*qf(.95,1,dfap)),
           (IP$DVDG[,2]*qf(.95,1,dfip)),(ID$DVDG[,2]*qf(.95,1,dfid)),(IDP$DVDG[,2]*qf(.95,1,dfidp)),
           (DP30ID$DVDG[,2]*qf(.95,1,dfdp30id)),(DP10ID$DVDG[,2]*qf(.95,1,dfdp10id)),
           (D90IDP$DVDG[,2]*qf(.95,1,dfd90idp)),(D2233$DVDG[,2]*qf(.95,1,df2233))))
dist <- DP$DVDG[,1]

plot((DP$DVDG[,2]*qf(.95,1,dfdp))~dist,ylim=yvdg,type="l",xlab="Distance",
     ylab=" ",col=4,lty=3,lwd=2.0)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})-\\hat{y}(\\mathbf{x}_{0})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}}$'),cex=.8,side=2,line=2,lwd=1)
lines((DP$DVDG[,3]*qf(.95,1,dfdp))~DP$DVDG[,1],col=4,lty=3,lwd=2)
lines((AP$DVDG[,2]*qf(.95,1,dfap))~AP$DVDG[,1],col="cyan",lty=1,lwd=2)
points((AP$DVDG[,2]*qf(.95,1,dfap))~D$DVDG[,1],col="cyan",pch=23,cex=1.2,bg="white")
lines((AP$DVDG[,3]*qf(.95,1,dfap))~AP$DVDG[,1],col="cyan",lty=1,lwd=2)
points((AP$DVDG[,3]*qf(.95,1,dfap))~D$DVDG[,1],col="cyan", pch=23,cex=1.2,bg="white")
lines((IP$DVDG[,2]*qf(.95,1,dfip))~IP$DVDG[,1],col=1,lty=1,lwd=2)
points((IP$DVDG[,2]*qf(.95,1,dfip))~IP$DVDG[,1],col=1,pch=19,cex=.8)
lines((IP$DVDG[,3]*qf(.95,1,dfip))~IP$DVDG[,1],col=1,lty=1,lwd=2)
points((IP$DVDG[,3]*qf(.95,1,dfip))~IP$DVDG[,1],col=1,pch=19,cex=.8)
lines((ID$DVDG[,2]*qf(.95,1,dfid))~ID$DVDG[,1],col=2,lty=1,lwd=2)      #CCD
points((ID$DVDG[,2]*qf(.95,1,dfid))~ID$DVDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((ID$DVDG[,3]*qf(.95,1,dfid))~ID$DVDG[,1],col=2,lty=1,lwd=2)
points((ID$DVDG[,3]*qf(.95,1,dfid))~ID$DVDG[,1],col=2,pch=22,cex=.8,bg="white")
lines((IDP$DVDG[,2]*qf(.95,1,dfidp))~IDP$DVDG[,1],col=6,lty=1,lwd=2)
points((IDP$DVDG[,2]*qf(.95,1,dfidp))~IDP$DVDG[,1],col=6,pch=8,cex=1)
lines((IDP$DVDG[,3]*qf(.95,1,dfidp))~IDP$DVDG[,1],col=6,lty=1,lwd=2)
points((IDP$DVDG[,3]*qf(.95,1,dfidp))~IDP$DVDG[,1],col=6,pch=8,cex=1)
lines((DP30ID$DVDG[,2]*qf(.95,1,dfdp30id))~DP30ID$DVDG[,1],col=3,lty=1,lwd=2)
points((DP30ID$DVDG[,2]*qf(.95,1,dfdp30id))~DP30ID$DVDG[,1],col=3,pch=24,cex=.8,bg="white")
lines((DP30ID$DVDG[,3]*qf(.95,1,dfdp30id))~DP30ID$DVDG[,1],col=3,lty=1,lwd=2)
points((DP30ID$DVDG[,3]*qf(.95,1,dfdp30id))~DP30ID$DVDG[,1],col=3,pch=24,cex=.8,bg="white")
lines((DP10ID$DVDG[,2]*qf(.95,1,dfdp10id))~DP10ID$DVDG[,1],col=4,lty=4,lwd=2)
lines((DP10ID$DVDG[,3]*qf(.95,1,dfdp10id))~DP10ID$DVDG[,1],col=4,lty=4,lwd=2)
lines((D90IDP$DVDG[,2]*qf(.95,1,dfd90idp))~D90IDP$DVDG[,1],col=7,lty=1,lwd=2)
lines((D90IDP$DVDG[,3]*qf(.95,1,dfd90idp))~D90IDP$DVDG[,1],col=7,lty=1,lwd=2)
lines((D2233$DVDG[,2]*qf(.95,1,df2233))~D2233$DVDG[,1],col=1,lty=1,lwd=2)
points((D2233$DVDG[,2]*qf(.95,1,df2233))~D2233$DVDG[,1],col=1,pch=8,cex=.8)
lines((D2233$DVDG[,3]*qf(.95,1,df2233))~D2233$DVDG[,1],col=1,lty=1,lwd=2)
points((D2233$DVDG[,3]*qf(.95,1,df2233))~D2233$DVDG[,1],col=1,pch=8,cex=.8)

legend(0,15,expression((DP)[S]),lty=3,col=4,bty="n",cex=.8,lwd=2)
legend(.2,15," ",col="cyan",pch=23,cex=.8,bg="white",bty="n")
legend(.2,15,expression((AP)[S]),lty=1,col="cyan",bty="n",cex=.8,lwd=2)
legend(0,14.1," ",col=1,pch=19,cex=.8,bty="n")
legend(0,14.1,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(.2,14.1," ",col=2,pch=22,cex=.8,bg="white",bty="n")
legend(.2,14.1,expression(I[D]/CCD),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(.4,14.1," ",col=6,pch=8,cex=.8,bty="n")
legend(.4,14.1,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0,13.2," ",col=3,pch=24,cex=.8,bg="white",bty="n")
legend(0,13.2,expression(8:kappa[1]==0.3~kappa[7]==0.7),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(.4,13.2,expression(9:kappa[1]==0.1~kappa[7]==0.9),lty=4,col=4,bty="n",cex=.8,lwd=2)
legend(0,13.2,expression(10:kappa[0]==0.9~kappa[8]==0.1),lty=1,col=7,bty="n",cex=.8,lwd=2)
legend(0,12.3," ",col=1,pch=8,cex=.8,bty="n")
legend(0, 12.3,expression(11:kappa[0]==~kappa[1]==0.2~kappa[4]==~kappa[7]==0.3),lty=1,col=1,bty="n",cex=.8,lwd=2)
}
#######################################################################################
# DVDG pe RELATIVE VOLUME (TO BE IN THE PAPER)
X11()
par(mfrow=c(1,1),bty="n")
{
yvdg=c(min((DP$DVDG[,3]*qf(.95,1,dfdp)),(AP$DVDG[,3]*qf(.95,1,dfap)),
           (IP$DVDG[,3]*qf(.95,1,dfip)),(ID$DVDG[,3]*qf(.95,1,dfid)),(IDP$DVDG[,3]*qf(.95,1,dfidp)),
           (DP30ID$DVDG[,3]*qf(.95,1,dfdp30id)),(DP10ID$DVDG[,3]*qf(.95,1,dfdp10id)),
           (D90IDP$DVDG[,3]*qf(.95,1,dfd90idp)),(D2233$DVDG[,3]*qf(.95,1,df2233))),
       max((DP$DVDG[,2]*qf(.95,1,dfdp)),(AP$DVDG[,2]*qf(.95,1,dfap)),
           (IP$DVDG[,2]*qf(.95,1,dfip)),(ID$DVDG[,2]*qf(.95,1,dfid)),(IDP$DVDG[,2]*qf(.95,1,dfidp)),
           (DP30ID$DVDG[,2]*qf(.95,1,dfdp30id)),(DP10ID$DVDG[,2]*qf(.95,1,dfdp10id)),
           (D90IDP$DVDG[,2]*qf(.95,1,dfd90idp)),(D2233$DVDG[,2]*qf(.95,1,df2233))))
plot((DP$DVDG[,2]*qf(.95,1,dfdp))~RV,ylim=c(0,14),type="n",xlab="Relative Volume",
     ylab=" ",col=4,lty=3,lwd=2.0)

lines((IP$DVDG[,2]*qf(.95,1,dfip))~RV,col=1,lty=1,lwd=2)
lines((IP$DVDG[,3]*qf(.95,1,dfip))~RV,col=1,lty=1,lwd=2)
lines((ID$DVDG[,2]*qf(.95,1,dfid))~RV,col=2,lty=1,lwd=2)      #CCD
lines((ID$DVDG[,3]*qf(.95,1,dfid))~RV,col=2,lty=1,lwd=2)
lines((IDP$DVDG[,2]*qf(.95,1,dfidp))~RV,col=6,lty=1,lwd=2)
lines((IDP$DVDG[,3]*qf(.95,1,dfidp))~RV,col=6,lty=1,lwd=2)
lines((DP30ID$DVDG[,2]*qf(.95,1,dfdp30id))~RV,col=3,lty=1,lwd=2)
lines((DP30ID$DVDG[,3]*qf(.95,1,dfdp30id))~RV,col=3,lty=1,lwd=2)
lines((DP10ID$DVDG[,2]*qf(.95,1,dfdp10id))~RV,col=4,lty=1,lwd=2)
lines((DP10ID$DVDG[,3]*qf(.95,1,dfdp10id))~RV,col=4,lty=1,lwd=2)
lines((D90IDP$DVDG[,2]*qf(.95,1,dfd90idp))~RV,col=7,lty=1,lwd=2)
lines((D90IDP$DVDG[,3]*qf(.95,1,dfd90idp))~RV,col=7,lty=1,lwd=2)

legend(0,14,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(.2,14,expression(I[D]/CCD),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(.4,14,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0,13.2,expression(8:kappa[1]==0.3~kappa[7]==0.7),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(.4,13.2,expression(9:kappa[1]==0.1~kappa[7]==0.9),lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0,12.3,expression(10:kappa[0]==0.9~kappa[8]==0.1),lty=1,col=7,bty="n",cex=.8,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})-\\hat{y}(\\mathbf{x}_{0})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}$'),cex=.8,side=2,line=2,lwd=1)
}
#######################################################################################

###################################################################################################
# DFDS PE DF (TO BE IN THE PAPER)
X11()
par(mfrow=c(1,1),bty="n")
{
yvdg=c(min((DP$FDS[,3]*qf(.95,1,dfdp)),(AP$FDS[,3]*qf(.95,1,dfap)),
           (IP$FDS[,3]*qf(.95,1,dfip)),(ID$FDS[,3]*qf(.95,1,dfid)),(IDP$FDS[,3]*qf(.95,1,dfidp)),
           (DP30ID$FDS[,3]*qf(.95,1,dfdp30id)),(DP10ID$FDS[,3]*qf(.95,1,dfdp10id)),
           (D90IDP$FDS[,3]*qf(.95,1,dfd90idp)),(D2233$FDS[,3]*qf(.95,1,df2233))),
       max((DP$FDS[,3]*qf(.95,1,dfdp)),(AP$FDS[,3]*qf(.95,1,dfap)),
           (IP$FDS[,3]*qf(.95,1,dfip)),(ID$FDS[,3]*qf(.95,1,dfid)),(IDP$FDS[,3]*qf(.95,1,dfidp)),
           (DP30ID$FDS[,3]*qf(.95,1,dfdp30id)),(DP10ID$FDS[,3]*qf(.95,1,dfdp10id)),
           (D90IDP$FDS[,3]*qf(.95,1,dfd90idp)),(D2233$FDS[,3]*qf(.95,1,df2233))))

plot((DP$FDS[,3]*qf(.95,1,dfdp))~DP$FDS[,1],ylim=c(0,14),type="n",xlab="Fraction of design space",
     ylab=" ",col=4,lty=3,lwd=2.0)
lines((IP$FDS[,3]*qf(.95,1,dfip))~IP$FDS[,1],col=1,lty=1,lwd=2)
lines((ID$FDS[,3]*qf(.95,1,dfid))~ID$FDS[,1],col=2,lty=1,lwd=2)
lines((IDP$FDS[,3]*qf(.95,1,dfidp))~IDP$FDS[,1],col=6,lty=1,lwd=2)
lines((DP30ID$FDS[,3]*qf(.95,1,dfdp30id))~DP30ID$FDS[,1],col=3,lty=1,lwd=2)
lines((DP10ID$FDS[,3]*qf(.95,1,dfdp10id))~DP10ID$FDS[,1],col=4,lty=1,lwd=2)
lines((D90IDP$FDS[,3]*qf(.95,1,dfd90idp))~D90IDP$FDS[,1],col=7,lty=1,lwd=2)

legend(0,14,expression((IP)),lty=1,col=1,bty="n",cex=.8,lwd=2)
legend(.2,14,expression(I[D]/CCD),lty=1,col=2,bty="n",cex=.8,lwd=2)
legend(.4,14,expression((I[D]*P)),lty=1,col=6,bty="n",cex=.8,lwd=2)
legend(0,13.2,expression(8:kappa[1]==0.3~kappa[7]==0.7),lty=1,col=3,bty="n",cex=.8,lwd=2)
legend(.4,13.2,expression(9:kappa[1]==0.1~kappa[7]==0.9),lty=1,col=4,bty="n",cex=.8,lwd=2)
legend(0,12.3,expression(10:kappa[0]==0.9~kappa[8]==0.1),lty=1,col=7,bty="n",cex=.8,lwd=2)
mtext(TeX('$\\sigma^{-2}\\textit{Var} \\left(\\hat{y}(\\mathbf{x})-\\hat{y}(\\mathbf{x}_{0})\\right)\\times \\textit{F}_{1,\\textit{d},0.95}$'),cex=.8,side=2,line=2,lwd=1)

}

