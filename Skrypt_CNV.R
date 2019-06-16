#ppkt.a

CNV <- read.table(file = "CNV.txt", header = TRUE)
NO_CNV <- read.table(file = "no_CNV.txt", header = TRUE, fill=TRUE)

#########################
attach(CNV)
wektorCNV <- c()

# oblicz œrednia wartoœæ dla ka¿dego ekotypu/probki
for(i in seq_along(CNV$PrÃ³bka)){
  wektorCNV[i]<- rowMeans(CNV[i,-(1:1)])
}
# wykresy
boxplot(wektorCNV, ylim = c(0,11), ylab = "œrednia waroœæ ekotypu", main = "œrednia wartoœæ próbek CNV")
plot(wektorCNV, type = "o", col = "pink", xlab = "próbka", ylab = "œrednia waroœæ ekotypu", main = "œrednia wartoœæ próbek CNV")

# stworz posortowana liste
tmp <- sort(wektorCNV, decreasing = TRUE)
head(tmp)
tail(tmp)

max <- which.max(wektorCNV)
wektorCNV[max]
#############################
attach(NO_CNV)
wektorNO_CNV <- c()

#oblicz œrednia wartoœæ dla ka¿dego ekotypu/probki
for(i in seq_along(NO_CNV$PrÃ³bka)){
  wektorNO_CNV[i]<- rowMeans(NO_CNV[i,-(1:1)])
}
# wykresy
boxplot(wektorNO_CNV, ylim = c(), ylab = "œrednia waroœæ ekotypu", main = "œrednia wartoœæ próbek NO_CNV")
plot(wektorNO_CNV, type = "o", col = "cyan", xlab = "próbka", ylab = "œrednia waroœæ ekotypu", main = "œrednia wartoœæ próbek NO_CNV")

# posortowana lista
tmp <- sort(wektorNO_CNV, decreasing = TRUE)
head(tmp)
tail(tmp)

max <- which.max(wektorNO_CNV)
wektorNO_CNV[max]

##############################
plot(wektorCNV, ylim = c(1.5,11), type = "o", col = "pink", xlab = "próbka", ylab = "œrednia waroœæ ekotypu", main = "œrednia wartoœæ próbek CNV")
lines(wektorNO_CNV, type = "o", col = "cyan", xlab = "próbka", ylab = "œrednia waroœæ ekotypu", main = "œrednia wartoœæ próbek NO_CNV")

wektorCNV[318]  #[1] 7.439302 œredni wynik kopii
wektorNO_CNV[318]  #[1] 2.025752

#############################
# CNV
counter <- 0
CNVL0318 <- c() # wektor do zapisu wartoœci z wiersza L0318
for (i in seq_along(CNV[318,-(1:1)])){
  if(CNV[318, i]==2){
    counter <- counter +1
  }
  CNVL0318[i]<-CNV[318,i+1] # przy okazji dodaje sobie Wartosci na wektor
}
print(counter) #==0 zawsze bo nie ma nigdzie dokladnie 2, sa blisko ale nie 2

plot(CNVL0318, type = "o", ylim = c(1,5) , col = "violet", xlab = "index genu", ylab = "wartoœæ CNV", main = "Wartoœæi cnv L0318 z CNV")

plot(wektorCNV, type = "o", col = "pink", xlab = "próbka", ylab = "œrednia waroœæ ekotypu", main = "œrednia wartoœæ próbek CNV")
lines(wektorCNV[318], x = 318, type = "o", col = "blue")

odstajace <- c()
for (i in seq_along(CNVL0318)){
  if(CNVL0318[i]>5){
    odstajace[length(odstajace)+1] <- names(CNV[i+1])  
  }
}
print(odstajace)
#####################################
# NO_CNV
counter <- 0
NO_CNVL0318 <- c() # wektor do zapisu wartoœci z wiersza L0318
for (i in seq_along(NO_CNV[318,-(1:1)])){
  if(NO_CNV[318, i]==2){
    counter <- counter +1
  }
  NO_CNVL0318[i]<-NO_CNV[318,i+1] # przy okazji dodaje sobie Wartosci na wektor
}
print(counter) #==0 zawsze bo nie ma nigdzie dokladnie 2, sa blisko ale nie 2

plot(NO_CNVL0318, type = "o", col = "cyan", xlab = "index genu", ylab = "wartoœæ CNV", main = "Wartoœæi cnv L0318 z no_CNV")

plot(wektorNO_CNV, type = "o", col = "pink", xlab = "próbka", ylab = "œrednia waroœæ ekotypu", main = "œrednia wartoœæ próbek no_CNV")
lines(wektorNO_CNV[318], x = 318, type = "o", col = "blue")

odstajace <- c()
for (i in seq_along(NO_CNVL0318)){
  if(NO_CNVL0318[i]>2.7){
    odstajace[length(odstajace)+1] <- names(NO_CNV[i]) 
  }
}
print(odstajace)

#########################################
#ppkt. b

# CNV
selekcja <- c()

for (i in seq_along(CNVL0318)){
  if(CNVL0318[i]<2.5 && CNVL0318[i]>1.5){
    selekcja[length(selekcja)+1] <- names(CNV[i+1]) 
  }
}
print(length(CNVL0318))
print(length(CNVL0318) - length(selekcja))  # roznica pomiedzy calym wierszem 318 
                                            # a tymi wyselekcjonowanymi
odstajace <- c()
odst <- c(1)
for (i in seq_along(CNVL0318)){ # przypadki odstajace
  if(CNVL0318[i]>2.5 || CNVL0318[i]<1.5){
    odstajace[length(odstajace)+1] <- names(CNV[i+1]) 
    odst[length(odst)+1] <- i 
  }
}
kopia_CNV <- CNV[,-odst]
sredniokopie <- colMeans(kopia_CNV)
plot(sredniokopie, type = "o", col = "pink" , main = "srednie wartoœci cnv po selekcji L0318 dla CNV", ylab = "srednia wartoœæ cnv", xlab = "index kolejnych genów")

print(odstajace)

# NO_CNV
selekcja2 <- c()

for (i in seq_along(NO_CNVL0318)){
  if(NO_CNVL0318[i]<2.5 && NO_CNVL0318[i]>1.5){
    selekcja2[length(selekcja2)+1] <- names(NO_CNV[i+1])
  }
}
print(length(NO_CNVL0318))
print(length(NO_CNVL0318) - length(selekcja2))  # roznica pomiedzy calym wierszem 318 
                                                # a tymi wyselekcjonowanymi

odstajace2 <- c()
odst2 <- c(1)
for (i in seq_along(NO_CNVL0318)){ # przypadki odstajace
  if(NO_CNVL0318[i]>2.5 || NO_CNVL0318[i]<1.5){
    odstajace2[length(odstajace2)+1] <- names(NO_CNV[i+1])
    odst2[length(odst2)+1] <- i 
  }
}
kopia_NO_CNV <- NO_CNV[,-odst2]
sredniokopie2 <- colMeans(kopia_NO_CNV)
plot(sredniokopie2, type = "o", col = "cyan" , main = "srednie wartoœci cnv po selekcji L0318 dla no_CNV", ylab = "srednia wartoœæ cnv", xlab = "index kolejnych genów")

print(odstajace2)

plot(sredniokopie, ylim = c(0,6)  ,type = "o", col = "pink" , main = "srednie wartoœci cnv po selekcji L0318 dla CNV", ylab = "srednia wartoœæ cnv", xlab = "index kolejnych genów")
lines(sredniokopie2, type = "o", col = "cyan" , main = "srednie wartoœci cnv po selekcji L0318 dla no_CNV", ylab = "srednia wartoœæ cnv", xlab = "index kolejnych genów")


########################################################
# ppkt. c

# CNV
cp <- sort(sredniokopie, decreasing = TRUE)
head(cp)
which.max(sredniokopie)

hist(sredniokopie, xlim = c(0,10), ylim = c(0,600), breaks = 1000, ylab = "czêstotliwoœæ" , xlab = "zmiennoœæ", main = "Zmiennoœæ liczby kopii CNV", col = c("pink", "cyan", "magenta", "red", "green"))


# NO_CNV
cp <- sort(sredniokopie2, decreasing = TRUE)
head(cp)
which.max(sredniokopie2)

hist(sredniokopie2, xlim = c(0,8),ylim = c(0,600), breaks = 6,  ylab = "czêstotliwoœæ" , xlab = "zmiennoœæ", main = "Zmienoœæ liczby kopii NO_CNV", col = c("pink", "cyan", "magenta", "red", "green"))







