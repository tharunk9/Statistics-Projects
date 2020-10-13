#Install and load necessary R packages
install.packages("MASS")
install.packages("fitdistrplus")
install.packages("goftest")
install.packages("extraDistr")

library(MASS)
library(fitdistrplus)
library(goftest)
library(extraDistr)

#----------------------------------------------------------------------------------------#

#Read in Google stock data
Gg <- read.csv(file="Googl.csv", header=TRUE, sep=",")

#Save adjusted closing price as a variable and calculate the log-return
G <- Gg$Adj.Close
Gr<- diff(log(G), lag = 1)

#Calculate the 0.05 sample quantile
quantile(Gr, 0.05)

#Calculate the conditional mean of the log-return values < 0.05 sample quantile
mean(Gr[Gr<quantile(Gr, 0.05)])

#Use the plotdist and descdist functions
plotdist(Gr, breaks = 25)
descdist(Gr)

#Use the fitdist function to fit "Laplace" distribution to the log-return values
Gfit <- fitdist(Gr, dist = "laplace", start = list(mu = 0, sigma = 1))

#Test the fit
ad.test(Gr, "plaplace", Gfit$estimate["mu"], Gfit$estimate["sigma"])

#Create a random sample of n observations 
#from the Laplace distribution with the parameters estimated from the fitdist function
n = 10000
Gsample <- rlaplace(n, mu = Gfit$estimate["mu"], sigma = Gfit$estimate["sigma"])
hist(Gsample, breaks = 50, main = "Histogram of a Random Sample from Laplace Distribution", freq = FALSE)

#Calculate the 0.05 sample quantile
Gc = quantile(Gsample,0.05)
Gc


#Calculate the conditional mean of the random sample < 0.05 sample quantile
f <- function(r){
  r*dlaplace(r, mu = Gfit$estimate["mu"], sigma = Gfit$estimate["sigma"]) / plaplace(Gc, mu = Gfit$estimate["mu"], sigma = Gfit$estimate["sigma"])
}
Gresult <- integrate(f, -Inf, Gc)
Gresult$value


#----------------------------------------------------------------------------------------#

#Read in Amazon stock data
AMZN <- read.csv(file="AMZN.csv", header=TRUE, sep=",")

#Save adjusted closing price as a variable and calculate the log-return
A <- AMZN$Adj.Close
Ar<- diff(log(A), lag = 1)

#Calculate the 0.05 sample quantile
quantile(Ar, 0.05)

#Calculate the conditional mean of the log-return values < 0.05 sample quantile
mean(Ar[Ar<quantile(Ar, 0.05)])

#Use the plotdist and descdist functions
plotdist(Ar, breaks = 25)
descdist(Ar)

#Use the fitdist function to fit "Huber" distribution to the log-return values
Afit <- fitdist(Ar, dist = "huber", start = list(mu = 0  , sigma = 1, epsilon = 1))

#Test the fit
ad.test(Ar, "phuber", Afit$estimate["mu"], Afit$estimate["sigma"], Afit$estimate["epsilon"])

#Create a random sample of n observations 
#from the Huber distribution with the parameters estimated from the fitdist function
n = 10000
Asample <- rhuber(n, mu = Afit$estimate["mu"], sigma = Afit$estimate["sigma"],epsilon=Afit$estimate["epsilon"])
hist(Asample, breaks = 50, main = "Histogram of a Random Sample from Huber Distribution", freq = FALSE)

#Calculate the 0.05 sample quantile
Ac <- quantile(Asample,0.05)
Ac

#Calculate the conditional mean of the random sample < 0.05 sample quantile
h <- function(r){
  r*dhuber(r, mu = Afit$estimate["mu"], sigma = Afit$estimate["sigma"],epsilon=Afit$estimate["epsilon"]) / phuber(Ac, mu = Afit$estimate["mu"], sigma = Afit$estimate["sigma"],epsilon=Afit$estimate["epsilon"])
}
Aresult <- integrate(h, -Inf, Ac)
Aresult$value


#----------------------------------------------------------------------------------------#

#Calculate the combined prices of Google and Amazon for each wE[0,1]
w = 0
w1 <- 1/48
counter <- 49

weightedprice <- function(w){ w*G + (1-w)*A}

C <- data.frame(V1 = weightedprice(w))
for (i in 1:counter){
  C[,i+1] <- weightedprice(w)
  w <- w+w1
}

#Calculate the log-return of the combined prices for each wE[0,1]
Cr = C[-1,]
counter = 50
for (i in 1:counter){
  Cr[,i] <- diff(log(C[,i]), lag = 1)
}


#Calculate the 5% quantiles for each wE[0,1]
Cq <-data.frame(V1 = matrix(vector(), 50, 1))
counter <- 50
alpha <- 0.05
for (i in 1:counter){
  Cq[i,] <- quantile(Cr[,i], alpha)
  
}

#Plot the 5% quantiles on grid of equally spaced 50 values in [0, 1]
Cq$w <- seq(0, 1, 1/49)
plot(Cq$w, Cq$V1,
     main = "0.05 Sample Quantile Q0.05(w)",
     xlab="w",
     ylab="Q0.05(w)")
grid(50,1)

Cq1 <- Cq$V1

#Calculate the conditional mean of the random sample < 0.05 sample quantile
Ce <-data.frame(V1 = matrix(vector(), 50, 1))
counter <- 50
for (i in 1:counter){
  Ce[i,] <- mean(Cr[Cr[,i] < Cq1[i], i] )
  
}

#Plot the expected values on grid of equally spaced 50 values in [0, 1]
Ce$w <- seq(0, 1, 1/49)
plot(Ce$w, Ce$V1,
     main = "Sample Estimates of ES(w)",
     xlab="w",
     ylab="ES(w)")
grid(50,1)