#install from our github branch
library(devtools)
install_github("gsbDBI/bayesTree",ref="modSA",force=TRUE)
library(BayesTree)
#code to test bayesTree
#simulate data (example from Friedman MARS paper)
f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
sigma = 1.0 #y = f(x) + sigma*z , z~N(0,1)
n = 100 #number of observations
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(x)
y=Ey+sigma*rnorm(n)
lmFit = lm(y~.,data.frame(x,y)) #compare lm fit to BART later
##run BART
#set dummy train weights for now
trainw<-y*0+0.5
weights_flag<-0
xtest<-matrix(0.0,0,0)
set.seed(99)
bartFit = bart(x,y,ndpost=200,xtest,trainw,weights_flag) #default is ndpost=1000, this is to run example fast.
par(mar=c(1,1,1,1))
plot(bartFit) # plot bart fit
##compare BART fit to linear matter and truth = Ey
fitmat = cbind(y,Ey,lmFit$fitted,bartFit$yhat.train.mean)
colnames(fitmat) = c('y','Ey','lm','bart')