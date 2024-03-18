#importujemy dane do ramki
data <- read.csv("heart_failure_clinical_records_dataset.csv", stringsAsFactors = FALSE)
#sprawdzamy strukturÄ™ zaimportowanego zbioru
str(data)
#wyÅ›wietlamy podstawowe statystyki
summary(data)
#sprawdzamy null
colSums(is.na(data))

table(data$DEATH_EVENT)

data$DEATH_EVENT <-factor(data$DEATH_EVENT)

normalize <- function(x) {
  return ((x-min(x)) / (max(x)-min(x)))
}

#normalizacja danych
data_n <- as.data.frame(lapply(data[1:12], normalize))
summary(data_n)

#sprawdzenie, czy normalizacja zadzia³a³a
summary(data_n$creatinine_phosphokinase)

data1 <- data_n[sample(1:nrow(data_n)),]
summary(data1)

data_train <- data1[1:229,]
data_test <- data1[230:299,]

data_train_labels <- data[1:229,13]
data_test_labels <- data[230:299,13]

#k=15
library(class)
data_test_pred <- knn(train = data_train, test = data_test, cl = data_train_labels, k = 15)
data_test_pred

library(gmodels)
CrossTable(x = data_test_labels, y = data_test_pred, prop.chisq = FALSE)

#k=7
library(class)
data_test_pred <- knn(train = data_train, test = data_test, cl = data_train_labels, k = 5)

library(gmodels)
CrossTable(x = data_test_labels, y = data_test_pred, prop.chisq = FALSE)

#k=23
library(class)
data_test_pred <- knn(train = data_train, test = data_test, cl = data_train_labels, k = 21)

library(gmodels)
CrossTable(x = data_test_labels, y = data_test_pred, prop.chisq = FALSE)

#k=1
library(class)
data_test_pred <- knn(train = data_train, test = data_test, cl = data_train_labels, k = 1)

library(gmodels)
CrossTable(x = data_test_labels, y = data_test_pred, prop.chisq = FALSE)

