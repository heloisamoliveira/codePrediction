#################################################################################################
######### PROGRAM FOR SEARCHING EXACT OPTIMUM UNBLOCKED DESIGN FOR INFERENCE  ################### 
#################################################################################################
source("functions_for_prediction_crit.R")

# NUMBER OF FACTORS
K <- 3

# LEVELS OF EACH FACTOR                          
Levels <- list(1:3,1:3,1:3) 

# DESIGN REGION: FOR CUBIC ENTER 'Y'  or 'N' FOR SPHERICAL     
Cubic <- 'Y'

# NUMBER OF DESIGN RUNS                       
N <- 26

# INDICATORS OF TERMS IN THE SECOND ORDER MODEL
# ENTRIES IN THE ORDER: LINEAR, QUADRATIC, INTERACTIONS  
Terms <- matrix(c(rep(1,K),rep(1,K),rep(1,(K*(K-1)/2))),nr=1)
Npar <- sum(Terms)

# NUMBER OF 'TRIES' FOR THE EXCHANGE ALGORITHM (RECOMMENDED TO START WITH 100 
# AND INCREASE IT IF THE BEST DESIGN DOES NOT SHOW UP OFTEN)
Ntries <- 100

# COMPOUND CRITERIA
# ENTER THE PRIORITY WEIGHT FOR EACH PROPERTY
# THE ENTRIES BELOW ILUSTRATE A COMPOUND OF "D" AND "I_DP" BOTH WITH THE SAME WEIGHT
Kap0 <- 1                 # weight for D 
Kap1 <- 0                 # weight for (DP)
Kap2 <- 0                 # weight for CIs (AP or LP)
Kap3 <- 0                 # weight for point estimation (A or L)
Kap4 <- 0                 # weight for DF efficiency (lof)
Kap5 <- 0                 # weight for point response predictions (I)
Kap6 <- 0                 # weight for CIs response predictions (IP)
Kap7 <- 0                 # weight for point prediction of differences (I_D)
Kap8 <- 1                # weight for CIs prediction of differences (I_DP)
Kappa <- c(Kap0,Kap1,Kap2,Kap3,Kap4,Kap5,Kap6,Kap7,Kap8)

# CONFIDENCE COEFFICIENT FOR DP
# (1-alpha1)
prob1 <- 0.95  

# CONFIDENCE COEFFICIENT FOR AP
# (1-alpha2)
prob2 <- 0.95

# CONFIDENCE COEFFICIENT FOR IP and I_DP
# (1-alpha3)
prob3 <- 0.95

# BONFERRONI'S CORRECTION OF ALPHA FOR AP CRITERION?                                  
MC <- 'N'      # 'Y'=yes or 'N'=no  
prob2 <- ifelse(MC=='Y',prob2^(1/Npar),prob2)

# WEIGHTS FOR PARAMETERS (W DIAGONAL)
# LEAVE AS IT IS FOR WEIGHTS SUGGESTED IN Gilmour & Trinca (2012)
# NOT USED FOR D, DP AND LOF CRITERIA
W <- matrix(rep(1,Npar),nr=1)                              
for(i in (K+1):(2*K))
{if(Terms[i]==1)
  W[sum(Terms[1:i])] <- W[sum(Terms[1:i])]/4
}
W <- matrix(W/sum(W),nr=1)

##############################################################################################
##### Warning: The implementation maximizes the criterion function which is writen ###########
##### in the form of reciprocals for A, AP, I, IP, I_D and I_DP quantities       #############
##############################################################################################

# TO RUN THE SEARCH
design <- SearchTreat()   

# TO SAVE AS R OBJECT
save(design,file="design1.RData") 

# TO SALVE THE RESULTING DESIGN AS txt FILE
write(t(design$XBest),file="design1.txt",ncol=K) 

###############################################################################################
################################## PRINTING THE RESULTS #######################################
###############################################################################################

print(design$time)
print(list("Weights actually used for A type criterion"=design$W))
print(list("Values of Kappa to criterion"=design$Kappa))
print("Criterion values for each Try")
stem(design$critall)
print(list("Criterion values, all parts"=design$criteria.opt))
print(list("Variances of each term for the Optimum Design"=design$V))
print("Optimal Design")
print(design$XBest)




