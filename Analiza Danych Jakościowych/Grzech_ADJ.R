wine <- read.csv("/Users/karolina/Documents/STUDIA/ADJ/winequalityN.csv")
wine1 <- read.csv("/Users/karolina/Documents/STUDIA/ADJ/winequalityN.csv")
#wine1 bedzie sluzylo do porownan
# tworzymy zmienną opdowiadajaca dobrej jakosci
# uznajemy ze jest to powyzej 6
wine <- na.omit(wine)
wine1 <- na.omit(wine)
wine$good <- factor(ifelse(wine$quality >= 6, "1","0"))
wine1$good <- factor(ifelse(wine$quality >= 6, "1","0"))
#musimy usunac quality
wine1 <- wine1[,-13]
wine <- wine[,-13]
wine_zmienne <- wine[,-13]
wine_zmienne <- wine_zmienne[,-1]
cor(subset(wine_zmienne))

summary(as.factor(wine$type))

summary(wine$free.sulfur.dioxide) #50 max
summary(wine$total.sulfur.dioxide) #260 max

wine <- subset(wine, free.sulfur.dioxide <= 50)
wine1 <- subset(wine1, free.sulfur.dioxide <= 50)
wine <- subset(wine, total.sulfur.dioxide <= 300)
wine1 <- subset(wine1, total.sulfur.dioxide <= 300)

summary(wine$pH)
wine1$pH <- factor(ifelse(wine$pH >= 3.8, "High",ifelse(wine$pH < 3.8 & wine$pH >= 3.3, "Recommended",
          ifelse(wine$pH < 3.3 & wine$pH >= 3.0, "Medium", "Low"))))

summary(wine$fixed.acidity)
wine1$fixed.acidity <- factor(ifelse(wine$fixed.acidity >= 9, "High",
          ifelse(wine$fixed.acidity < 9 & wine$fixed.acidity >= 6, "Medium","Low")))

summary(wine$volatile.acidity)
wine1$volatile.acidity <- factor(ifelse(wine$volatile.acidity >= 0.7, "High","Low"))

summary(wine$citric.acid)
wine1$citric.acid <- factor(ifelse(wine$citric.acid >= 1.0, "5",
                    ifelse(wine$citric.acid < 1.0 & wine$citric.acid >= 0.7, "4",
                    ifelse(wine$citric.acid < 0.7 & wine$citric.acid >= 0.55, "3",
                    ifelse(wine$citric.acid < 0.55 & wine$citric.acid >= 0.3, "2","1")))))

summary(wine$residual.sugar)
wine1$residual.sugar <- factor(ifelse(wine$residual.sugar >= 45, "Sweet",
            ifelse(wine$residual.sugar < 45 & wine$residual.sugar >= 12, "Semi-sweet",
            ifelse(wine$residual.sugar < 12 & wine$residual.sugar >= 4, "Semi-dry", "Dry"))))

summary(wine$chlorides)
wine1$chlorides <- factor(ifelse(wine$chlorides >= 0.1, "6",
                  ifelse(wine$chlorides < 0.1 & wine$chlorides >= 0.069, "5",
                  ifelse(wine$chlorides < 0.069 & wine$chlorides >= 0.056, "4",
                  ifelse(wine$chlorides < 0.056 & wine$chlorides >= 0.048, "3",
                  ifelse(wine$chlorides < 0.048 & wine$chlorides >= 0.03, "2","1"))))))

wine1$free.sulfur.dioxide <- factor(ifelse(wine$free.sulfur.dioxide >= 30, "High",
ifelse(wine$free.sulfur.dioxide < 30 & wine$free.sulfur.dioxide >= 15, "Medium","Low")))

wine1$total.sulfur.dioxide <- factor(ifelse(wine$total.sulfur.dioxide >= 200, "High",
ifelse(wine$total.sulfur.dioxide < 200 & wine$total.sulfur.dioxide >= 100, "Medium","Low")))

summary(wine$density)
wine1$density <- factor(ifelse(wine$density >= 0.99, "Normal","Sparkling"))

summary(wine$sulphates)
wine1$sulphates <- factor(ifelse(wine$sulphates >= 1.5, "6",
          ifelse(wine$sulphates < 1.5 & wine$sulphates >= 1.0, "5",
          ifelse(wine$sulphates < 1.0 & wine$sulphates >= 0.5, "4",
          ifelse(wine$sulphates < 0.5 & wine$sulphates >= 0.4, "3",
          ifelse(wine$sulphates < 0.4 & wine$sulphates >= 0.3, "2","1"))))))

summary(wine$alcohol)
wine1$alcohol <- factor(ifelse(wine$alcohol >= 12, "High",
          ifelse(wine$alcohol < 12 & wine$alcohol >= 10, "Medium","Low")))


########################################################################
odds.ratio <-
  function(x, pad.zeros=FALSE, conf.level=0.95) {
    if (pad.zeros) {
      if (any(x==0)) x <- x + 0.5
    }
    theta <- x[1,1] * x[2,2] / ( x[2,1] * x[1,2] )
    ASE <- sqrt(sum(1/x))
    CI <- exp(log(theta)
              + c(-1,1) * qnorm(0.5*(1+conf.level)) *ASE )
    p.value <- pnorm(abs(log(theta)/ASE), lower.tail=FALSE)
    list(estimator=theta,
         p.value=p.value,
         conf.interval=CI,
         conf.level=conf.level)
  }
GK.tau <- function(dat)
{ N <- sum(dat);

dat.rows <- nrow(dat);
dat.cols <- ncol(dat);
max.col <- sum.col <- L.col <- matrix(,dat.cols);
max.row <- sum.row <- L.row <- matrix(,dat.rows);
for(i in 1:dat.cols)
{ sum.col[i] <- sum(dat[,i]); max.col[i] <- max(dat[,i]); }
for(i in 1:dat.rows)
{ sum.row[i] <- sum(dat[i,]); max.row[i] <- max(dat[i,]); }

max.row.margin <- max(apply(dat,1,sum));   max.col.margin <- max(apply(dat,2,sum));

# Goodman-Kruskal tau (raws=indep.vars, cols=dep.vars)
n.err.unconditional <- N^2;
for(i in 1:dat.rows)
  n.err.unconditional <- n.err.unconditional-N*sum(dat[i,]^2/sum.row[i]);   
n.err.conditional <- N^2-sum(sum.col^2);   
tau <- 1-(n.err.unconditional/n.err.conditional);

v <- n.err.unconditional/(N^2);
d <- n.err.conditional/(N^2);
f <- d*(v+1)-2*v;

var.tau.CR <- 0;
for(i in 1:dat.rows)
  for(j in 1:dat.cols)
    var.tau.CR <- var.tau.CR + dat[i,j]*(-2*v*(sum.col[j]/N)+d*((2*dat[i,j]/sum.row[i])-sum((dat[i,]/sum.row[i])^2))-f)^2/(N^2*d^4);
ASE <- sqrt(var.tau.CR);

U.tau.CR <- (N-1)*(dat.cols-1)*tau; 
# Chi-squared approximation for H0 according to Margolin & Light JASA 1974, 755-764, 
# see also Liebetrau 1983   
p.value <- pchisq(U.tau.CR,df=(dat.rows-1)*(dat.cols-1),lower=FALSE); 

data.frame(tau, p.value, ASE);  
}

#############################################################################
summary(wine$good)
# 2372 bad, 4091 good
sum(wine$type == "white" & wine$good == "1")
sum(wine$type == "white" & wine$good == "0")
sum(wine$type == "red" & wine$good == "1")
sum(wine$type == "red" & wine$good == "0")

#dla white
binconf(2743,4005,method="asymptotic",alpha=0.05)
binconf(2743+2,4005+4,method="asymptotic",alpha=0.05)
binconf(2743,4005,method="wilson")
binconf(2743,4005,method="exact")
#dla red
binconf(844,1577,method="asymptotic",alpha=0.05)
binconf(844+2,1577+4,method="asymptotic",alpha=0.05)
binconf(844,1577,method="wilson")
binconf(844,1577,method="exact")
wiersze <- c("White","Red")
kolumny <- c("1","0")
dane <- c(2743,844,1262,733)
ocena <- cbind(expand.grid(type=wiersze, ocena=kolumny),liczn=dane)
tab1 <- xtabs(liczn~type+ocena,data=ocena)
tab1
odds.ratio(tab1)
chisq.test(tab1)
GK.tau(tab1)
#################################################################

sum(wine1$alcohol == "High" & wine1$good == "1")
sum(wine1$alcohol == "High" & wine1$good == "0")
sum(wine1$alcohol == "Medium" & wine1$good == "1")
sum(wine1$alcohol == "Medium" & wine1$good == "0")
sum(wine1$alcohol == "Low" & wine1$good == "1")
sum(wine1$alcohol == "Low" & wine1$good == "0")

#dla High
binconf(874,932,method="asymptotic",alpha=0.05)
binconf(874+2,932+4,method="asymptotic",alpha=0.05)
binconf(874,932,method="wilson")
binconf(874,932,method="exact")

#dla Medium
binconf(1868,2625,method="asymptotic",alpha=0.05)
binconf(1868+2,2625+4,method="asymptotic",alpha=0.05)
binconf(1868,2625,method="wilson")
binconf(1868,2625,method="exact")

#dla Low
binconf(845,2025,method="asymptotic",alpha=0.05)
binconf(845+2,2025+4,method="asymptotic",alpha=0.05)
binconf(845,2025,method="wilson")
binconf(845,2025,method="exact")

wiersze <- c("High","Medium","Low")
kolumny <- c("1","0")
dane <- c(874,1868,845,58,757,1180)
ocena <- cbind(expand.grid(alcohol=wiersze, ocena=kolumny),liczn=dane)
tab2 <- xtabs(liczn~alcohol+ocena,data=ocena)
tab2

chisq.test(tab2)
GK.tau(tab2)

#################################################################

sum(wine1$fixed.acidity == "High" & wine1$good == "1")
sum(wine1$fixed.acidity == "High" & wine1$good == "0")
sum(wine1$fixed.acidity == "Medium" & wine1$good == "1")
sum(wine1$fixed.acidity == "Medium" & wine1$good == "0")
sum(wine1$fixed.acidity == "Low" & wine1$good == "1")
sum(wine1$fixed.acidity == "Low" & wine1$good == "0")

wiersze <- c("High","Medium","Low")
kolumny <- c("1","0")
dane <- c(317,2860,410,228,1622,145)
ocena <- cbind(expand.grid(fixed.acidity=wiersze, ocena=kolumny),liczn=dane)
tab3 <- xtabs(liczn~fixed.acidity+ocena,data=ocena)
tab3

chisq.test(tab3)
GK.tau(tab3)

#################################################################

sum(wine1$volatile.acidity == "High" & wine1$good == "1")
sum(wine1$volatile.acidity == "High" & wine1$good == "0")
sum(wine1$volatile.acidity == "Low" & wine1$good == "1")
sum(wine1$volatile.acidity == "Low" & wine1$good == "0")

wiersze <- c("High","Low")
kolumny <- c("1","0")
dane <- c(80,3507,168,1827)
ocena <- cbind(expand.grid(volatile.acidity=wiersze, ocena=kolumny),liczn=dane)
tab4 <- xtabs(liczn~volatile.acidity+ocena,data=ocena)
tab4

odds.ratio(tab4)
chisq.test(tab4)
GK.tau(tab4)

#################################################################

sum(wine1$citric.acid == "5" & wine1$good == "1")
sum(wine1$citric.acid == "5" & wine1$good == "0")
sum(wine1$citric.acid == "4" & wine1$good == "1")
sum(wine1$citric.acid == "4" & wine1$good == "0")
sum(wine1$citric.acid == "3" & wine1$good == "1")
sum(wine1$citric.acid == "3" & wine1$good == "0")
sum(wine1$citric.acid == "2" & wine1$good == "1")
sum(wine1$citric.acid == "2" & wine1$good == "0")
sum(wine1$citric.acid == "1" & wine1$good == "1")
sum(wine1$citric.acid == "1" & wine1$good == "0")

wiersze <- c("5","4","3","2","1")
kolumny <- c("1","0")
dane <- c(5,44,117,1982,1439,1,28,89,831,1046)
ocena <- cbind(expand.grid(citric.acid=wiersze, ocena=kolumny),liczn=dane)
tab5 <- xtabs(liczn~citric.acid+ocena,data=ocena)
tab5

GK.tau(tab5)
#################################################################

sum(wine1$residual.sugar == "Sweet" & wine1$good == "1")
sum(wine1$residual.sugar == "Sweet" & wine1$good == "0")
sum(wine1$residual.sugar == "Semi-sweet" & wine1$good == "1")
sum(wine1$residual.sugar == "Semi-sweet" & wine1$good == "0")
sum(wine1$residual.sugar == "Semi-dry" & wine1$good == "1")
sum(wine1$residual.sugar == "Semi-dry" & wine1$good == "0")
sum(wine1$residual.sugar == "Dry" & wine1$good == "1")
sum(wine1$residual.sugar == "Dry" & wine1$good == "0")

wiersze <- c("Semi-sweet","Semi-dry","Dry")
kolumny <- c("1","0")
dane <- c(359,1113,2114,217,543,1235)
ocena <- cbind(expand.grid(residual.sugar=wiersze, ocena=kolumny),liczn=dane)
tab6 <- xtabs(liczn~residual.sugar+ocena,data=ocena)
tab6

chisq.test(tab6)
GK.tau(tab6)
fisher.test(tab6,simulate.p.value = T,B=10000)

#################################################################

sum(wine1$chlorides == "6" & wine1$good == "1")
sum(wine1$chlorides == "6" & wine1$good == "0")
sum(wine1$chlorides == "5" & wine1$good == "1")
sum(wine1$chlorides == "5" & wine1$good == "0")
sum(wine1$chlorides == "4" & wine1$good == "1")
sum(wine1$chlorides == "4" & wine1$good == "0")
sum(wine1$chlorides == "3" & wine1$good == "1")
sum(wine1$chlorides == "3" & wine1$good == "0")
sum(wine1$chlorides == "2" & wine1$good == "1")
sum(wine1$chlorides == "2" & wine1$good == "0")
sum(wine1$chlorides == "1" & wine1$good == "1")
sum(wine1$chlorides == "1" & wine1$good == "0")

wiersze <- c("6","5","4","3","2","1")
kolumny <- c("1","0")
dane <- c(134,542,354,456,1733,368,169,568,236,362,600,60)
ocena <- cbind(expand.grid(chlorides=wiersze, ocena=kolumny),liczn=dane)
tab7 <- xtabs(liczn~chlorides+ocena,data=ocena)
tab7

GK.tau(tab7)

#################################################################

sum(wine1$free.sulfur.dioxide == "High" & wine1$good == "1")
sum(wine1$free.sulfur.dioxide == "High" & wine1$good == "0")
sum(wine1$free.sulfur.dioxide == "Medium" & wine1$good == "1")
sum(wine1$free.sulfur.dioxide == "Medium" & wine1$good == "0")
sum(wine1$free.sulfur.dioxide == "Low" & wine1$good == "1")
sum(wine1$free.sulfur.dioxide == "Low" & wine1$good == "0")

wiersze <- c("High","Medium","Low")
kolumny <- c("1","0")
dane <- c(1558,1367,662,646,728,621)
ocena <- cbind(expand.grid(free.sulfur.dioxide=wiersze, ocena=kolumny),liczn=dane)
tab8 <- xtabs(liczn~free.sulfur.dioxide+ocena,data=ocena)
tab8

chisq.test(tab8)
GK.tau(tab8)

#################################################################

sum(wine1$total.sulfur.dioxide == "High" & wine1$good == "1")
sum(wine1$total.sulfur.dioxide == "High" & wine1$good == "0")
sum(wine1$total.sulfur.dioxide == "Medium" & wine1$good == "1")
sum(wine1$total.sulfur.dioxide == "Medium" & wine1$good == "0")
sum(wine1$total.sulfur.dioxide == "Low" & wine1$good == "1")
sum(wine1$total.sulfur.dioxide == "Low" & wine1$good == "0")

wiersze <- c("High","Medium","Low")
kolumny <- c("1","0")
dane <- c(96,2021,1470,80,1051,864)
ocena <- cbind(expand.grid(alcohol=wiersze, ocena=kolumny),liczn=dane)
tab9 <- xtabs(liczn~alcohol+ocena,data=ocena)
tab9

chisq.test(tab9)
GK.tau(tab9)

#################################################################

sum(wine1$density == "Normal" & wine1$good == "1")
sum(wine1$density == "Normal" & wine1$good == "0")
sum(wine1$density == "Sparkling" & wine1$good == "1")
sum(wine1$density == "Sparkling" & wine1$good == "0")

wiersze <- c("Normal","Sparkling")
kolumny <- c("1","0")
dane <- c(3279,308,1981,14)
ocena <- cbind(expand.grid(density=wiersze, ocena=kolumny),liczn=dane)
tab10 <- xtabs(liczn~density+ocena,data=ocena)
tab10

odds.ratio(tab10)
chisq.test(tab10)
GK.tau(tab10)

#################################################################

sum(wine1$pH == "High" & wine1$good == "1")
sum(wine1$pH == "High" & wine1$good == "0")
sum(wine1$pH == "Recommended" & wine1$good == "1")
sum(wine1$pH == "Recommended" & wine1$good == "0")
sum(wine1$pH == "Medium" & wine1$good == "1")
sum(wine1$pH == "Medium" & wine1$good == "0")
sum(wine1$pH == "Low" & wine1$good == "1")
sum(wine1$pH == "Low" & wine1$good == "0")

#dla High
binconf(8,9,method="asymptotic",alpha=0.05)
binconf(8+2,9+4,method="asymptotic",alpha=0.05)
binconf(8,9,method="wilson")
binconf(8,9,method="exact")

#dla Recommended
binconf(1194,1815,method="asymptotic",alpha=0.05)
binconf(1194+2,1815+4,method="asymptotic",alpha=0.05)
binconf(1194,1815,method="wilson")
binconf(1194,1815,method="exact")

#dla Medium
binconf(2116,3351,method="asymptotic",alpha=0.05)
binconf(2116+2,3351+4,method="asymptotic",alpha=0.05)
binconf(2116,3351,method="wilson")
binconf(2116,3351,method="exact")

#dla Low
binconf(3507,5334,method="asymptotic",alpha=0.05)
binconf(3507+2,5334+4,method="asymptotic",alpha=0.05)
binconf(3507,5334,method="wilson")
binconf(3507,5334,method="exact")

wiersze <- c("High","Recommended","Medium","Low")
kolumny <- c("1","0")
dane <- c(8,1194,2116,269,1,621,1235,138)
ocena <- cbind(expand.grid(pH=wiersze, ocena=kolumny),liczn=dane)
tab11 <- xtabs(liczn~pH+ocena,data=ocena)
tab11

GK.tau(tab11)
#nie mamy podstaw do odrzucenia H0
fisher.test(tab11,simulate.p.value = T,B=10000)

#################################################################

sum(wine1$sulphates == "6" & wine1$good == "1")
sum(wine1$sulphates == "6" & wine1$good == "0")
sum(wine1$sulphates == "5" & wine1$good == "1")
sum(wine1$sulphates == "5" & wine1$good == "0")
sum(wine1$sulphates == "4" & wine1$good == "1")
sum(wine1$sulphates == "4" & wine1$good == "0")
sum(wine1$sulphates == "3" & wine1$good == "1")
sum(wine1$sulphates == "3" & wine1$good == "0")
sum(wine1$sulphates == "2" & wine1$good == "1")
sum(wine1$sulphates == "2" & wine1$good == "0")
sum(wine1$sulphates == "1" & wine1$good == "1")
sum(wine1$sulphates == "1" & wine1$good == "0")

wiersze <- c("6","5","4","3","2","1")
kolumny <- c("1","0")
dane <- c(3,28,1981,982,560,33,5,27,1087,596,265,15)
ocena <- cbind(expand.grid(sulphates=wiersze, ocena=kolumny),liczn=dane)
tab12 <- xtabs(liczn~sulphates+ocena,data=ocena)
tab12

GK.tau(tab12)

###############################################################

full.model <- glm(good ~ ., family = binomial, data = wine)
summary(full.model)
anova(full.model, test = "Chisq")

#sprobujmy usunac zmienne o najwiekszej p-wartosci
model1 <- glm(good ~ type+volatile.acidity+citric.acid+residual.sugar+free.sulfur.dioxide+
          total.sulfur.dioxide+sulphates+alcohol+chlorides+density+pH,
          family = binomial, data = wine)
summary(model1)
anova(model1, test = "Chisq")

model2 <- glm(good ~ type+volatile.acidity+citric.acid+residual.sugar+free.sulfur.dioxide+
          total.sulfur.dioxide+sulphates+alcohol+chlorides+density,
          family = binomial, data = wine)
summary(model2)
anova(model2, test = "Chisq")

model3 <- glm(good ~ type+volatile.acidity+residual.sugar+free.sulfur.dioxide+
          total.sulfur.dioxide+sulphates+alcohol+chlorides+density,
          family = binomial, data = wine)
summary(model3)
anova(model3, test = "Chisq")

model4 <- glm(good ~ type+volatile.acidity+residual.sugar+free.sulfur.dioxide+
          total.sulfur.dioxide+sulphates+alcohol+density,
          family = binomial, data = wine)
summary(model4)
anova(model4, test = "Chisq")

model3 <- glm(good ~ type+volatile.acidity+residual.sugar+free.sulfur.dioxide+
          total.sulfur.dioxide+sulphates+alcohol+chlorides+density,
          family = binomial, data = wine)
summary(model3)
anova(model3, test = "Chisq")

model4 <- glm(good ~ type+volatile.acidity+residual.sugar+free.sulfur.dioxide+
          total.sulfur.dioxide+sulphates+alcohol+density,
          family = binomial, data = wine)
summary(model4)
anova(model4, test = "Chisq")

model5 <- glm(good ~ type+volatile.acidity+residual.sugar+free.sulfur.dioxide+
          total.sulfur.dioxide+sulphates+alcohol,
          family = binomial, data = wine)
summary(model5)
anova(model5, test = "Chisq")

#sprobujmy type
model6 <- glm(good ~ volatile.acidity+residual.sugar+free.sulfur.dioxide+
          total.sulfur.dioxide+sulphates+alcohol,
          family = binomial, data = wine)
summary(model6)
anova(model6, test = "Chisq")

# sprawdzmy jeszcze model z interakcjami
model7 <- glm(good ~ volatile.acidity*residual.sugar*free.sulfur.dioxide*
          total.sulfur.dioxide*sulphates*alcohol,
          family = binomial, data = wine)
summary(model7)
anova(model6,model7,test="Chisq")

model8 <- glm(good ~ residual.sugar + volatile.acidity:residual.sugar + volatile.acidity:sulphates +
          volatile.acidity:residual.sugar:sulphates + volatile.acidity:residual.sugar:alcohol +
          volatile.acidity:residual.sugar:free.sulfur.dioxide:sulphates +
          volatile.acidity:residual.sugar:total.sulfur.dioxide:sulphates +
          volatile.acidity:residual.sugar:sulphates:alcohol, family = binomial, data = wine)
summary(model8)

model9 <- glm(good ~ residual.sugar + volatile.acidity:residual.sugar + volatile.acidity:sulphates +
          volatile.acidity:residual.sugar:alcohol +
          volatile.acidity:residual.sugar:free.sulfur.dioxide:sulphates +
          volatile.acidity:residual.sugar:total.sulfur.dioxide:sulphates +
          volatile.acidity:residual.sugar:sulphates:alcohol, family = binomial, data = wine)
summary(model9)

model0 <- glm(good ~ 1, family = binomial, data = wine)

#sprawdzamy dewiancje oraz AIC
c(model6$deviance,model7$deviance)
c(model6$aic,model7$aic)

c(model0$deviance, full.model$deviance, model1$deviance, model2$deviance, model3$deviance,
model4$deviance, model5$deviance, model6$deviance, model7$deviance, model8$deviance,model9$deviance)
which.min(c(model0$deviance, full.model$deviance, model1$deviance, model2$deviance, model3$deviance,
model4$deviance, model5$deviance, model6$deviance, model7$deviance, model8$deviance,model9$deviance))

c(model0$aic, full.model$aic, model1$aic, model2$aic, model3$aic,
model4$aic, model5$aic, model6$aic, model7$aic, model8$aic,model9$aic)
which.min(c(model0$aic, full.model$aic, model1$aic, model2$aic, model3$aic,
model4$aic, model5$aic, model6$aic, model7$aic, model8$aic,model9$aic))

#model7 ma najlepsze parametry
model7 <- glm(good ~ volatile.acidity*residual.sugar*free.sulfur.dioxide*
          total.sulfur.dioxide*sulphates*alcohol,
          family = binomial, data = wine)
summary(model7)

library(lmtest)
lrtest(model0,model7)
anova(model0,model7,test="Chisq")

#sprawdzmy dokladnosc modelu
set.seed(46)
#tworzymy zbior treningowy i testowy
prop <- 0.7
# Wygenerowanie indeksów dla zbioru treningowego
train_idx <- sample(1:nrow(wine), prop * nrow(wine))
# Podział danych na zbiór treningowy i testowy
train <- wine[train_idx, ]
test <- wine[-train_idx, ]
 
modell <- glm(good ~ volatile.acidity*residual.sugar*free.sulfur.dioxide*
                total.sulfur.dioxide*sulphates*alcohol,
               family = binomial, data = train)
        
wine.predict <- predict(modell, newdata=test, type="response")
labels2 <- ifelse(wine.predict > 0.5, 1, 0)
tabela <- data.frame(Actual = test$good, Predicted = labels2)
conf_matrix2 <- table(Przewidziane = tabela$Predicted, Aktualne = tabela$Actual)
conf_matrix2

accuracy2 <- sum(diag(conf_matrix2)) / sum(conf_matrix2)
print(paste("Dokładność modelu:", round(accuracy2, 4)))
