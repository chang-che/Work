##########Q1
####Method 1
mydat <-read.delim('http://www.ams.sunysb.edu/~pfkuan/Teaching/AMS597/Data/Qn1Data.txt',header=T)
BIC.all <- NULL
X <- mydat[,-1]

BIC.all <- BIC(glm(mydat$y~1,family='binomial'))
names(BIC.all) <- c('Intercept')
for(i in 1:6){
  tmp <- combn(6,i)
  for(j in 1:dim(tmp)[2]){
    subX <- X[,tmp[,j]]
    subdat <- data.frame(cbind(mydat$y,subX))
    colnames(subdat)[1] <- 'y'
    BIC.all <- c(BIC.all,BIC(glm(y~.,data=subdat,family='binomial')))
    names(BIC.all)[length(BIC.all)] <- paste('x',tmp[,j],collapse='',sep='')
  }
}
BIC.all[which(BIC.all==min(BIC.all))]
#> BIC.all[which(BIC.all==min(BIC.all))]
#      x3 
#266.9233
# y~x3 is the best model using BIC

####Method 2
library(leaps);#install.packages("leaps");
mydat <-read.delim('http://www.ams.sunysb.edu/~pfkuan/Teaching/AMS597/Data/Qn1Data.txt',header=T)
leaps1<-regsubsets(y~x1+x2+x3+x4+x5+x6,data=mydat,nbest=choose(6,6/2))
summary(leaps1)
which.min(summary(leaps1)$bic);

##########Q2
### This is a Weibull(2,1) distribution
## F(x)= 1-exp(-x^2)
u <- runif(1000)
x <- sqrt(-log(u)) ## sqrt(-log(1-u))
## check your results with built in Weibull
qqplot(x,rweibull(1000,2,1));abline(0,1)

##########Q3
u1 <- runif(1500)
u2 <- runif(1500)
v <- sqrt(-2*log(u1))*cos(2*pi*u2)
w <- sqrt(-2*log(u1))*sin(2*pi*u2)
z1 <- matrix(v[1:300],ncol=3)
z2 <- matrix(v[301:800],ncol=5)
z3 <- matrix(v[801:1500],ncol=7)
w1 <- rowSums(z1^2)
w2 <- rowSums(z2^2)
w3 <- rowSums(z3^2)
t1 <- w[1:100]/sqrt(w1/3)
t2 <- w[101:200]/sqrt(w2/5)
t3 <- w[201:300]/sqrt(w3/7)
s <- sample(c(1:3),100,prob=c(0.3,0.35,0.35),replace=TRUE)
x<-rep(NA,100)
x[s==1] <- t1[1:sum(s==1)]
x[s==2] <- t2[1:sum(s==2)]
x[s==3] <- t3[1:sum(s==3)]
x
##check your results 
a<-rep(NA,100);b=c(3,5,7);for(i in 1:3){a[s==i]<-rt(length(which(s==i)),b[i])};qqplot(x,a); abline(0,1,col='red')

