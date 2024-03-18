MC <- 1000
m <- 80
n <- 100
N <- 180

stat <- function(X,Y,m,n) {
  Z <- c(X,Y)
  T1 <- c()
  T2 <- c()
  for(i in 1:m) {
    T1[i] = -1/m 
  }
  for(i in (m+1):(m+n)) {
    T1[i] = 1/n
  }
  V1 <- cbind(Z,T1)
  V1 = V1[order(V1[,1],decreasing = FALSE),]
  T1 = V1[,2]
  SP1 <- c()
  SP1[1] = T1[1]
  for(i in 2:(m+n)) {
    SP1[i] = SP1[i-1]+T1[i]
  }
  KSY <- max(SP1)
  
  for(i in 1:m) {
    T2[i] = 1/m 
  }
  for(i in (m+1):(m+n)) {
    T2[i] = -1/n
  }
  V2 <- cbind(Z,T2)
  V2 = V2[order(V2[,1],decreasing = FALSE),]
  T2 = V2[,2]
  SP2 <- c()
  SP2[1] = T2[1]
  for(i in 2:(m+n)) {
    SP2[i] = SP2[i-1]+T2[i]
  }
  KSX <- max(SP2)
  KS  <- min(KSY,KSX)
  return(KS)
}

KSMC <- c()
KSMC1 <- c()
KSMC2 <- c()

for(j in 1:MC) {
  ZMC <- runif(N,0,1)
  TMC1 <- c()
  for(i in 1:N) {
    if(i <= m) {TMC1[i] = -1/m}
    else{TMC1[i] = 1/n}
  }
  VMC1 <- cbind(ZMC,TMC1)
  VMC1 = VMC1[order(VMC1[,1], decreasing = FALSE),]
  TMC1 = VMC1[,2]
  SPMC1 <- c()
  SPMC1[1] = TMC1[1]
  for(i in 2:N) {
    SPMC1[i] = SPMC1[i-1] + TMC1[i]
  }
  KSMC1[j] <- max(SPMC1)
  
  TMC2 <- c()
  for(i in 1:N) {
    if(i <= m) {TMC2[i] = 1/m}
    else{TMC2[i] = -1/n}
  }
  VMC2 <- cbind(ZMC,TMC2)
  VMC2 = VMC2[order(VMC2[,1], decreasing = FALSE),]
  TMC2 = VMC2[,2]
  SPMC2 <- c()
  SPMC2[1] = TMC2[1]
  for(i in 2:N) {
    SPMC2[i] = SPMC2[i-1] + TMC2[i]
  }
  KSMC2[j] <- max(SPMC2)
  KSMC[j] <- min(KSMC1[j],KSMC2[j])
}

plot(ecdf(KSMC))

p <- 75
alpha <- 0.05
s <- c()
M <- c()
Moc <- c()
pKS <- c()
KS <- c()
m <- 80
n <- 100

for(t in 1:p) {
  s[t] <- (t-1)*0.05+0.25
  for(i in 1:MC) {
    X = rnorm(m, mean = 0, sd = 1)
    Y = rnorm(n, mean = 0, sd = s[t])
    pKS <- 0
    KS[i] <- stat(X,Y,m,n)
    
    for(j in 1:MC) {
      if(KSMC[j] > KS[i]) {pKS = pKS + 1}
    }
    pKS[i] = pKS/MC
    
    M[i] <- 0
    if(pKS[i] < alpha) {M[i] <- 1}
  }
  Moc[t] <- mean(M)
  
}
plot(s,Moc)
abline(h=0.05)



