#Does Metropolis-Hasting algorithm work well on 
#the simulated data from the first order Markov chain model?
#Name: Eunji Kim(i6073884) with: Wen Xu(i6054652)
#Date: 25/01/2018

set.seed(123) # setting random seed 12345 to be able to replicate results

#input: x is old(scalar)
#       y is new(scalar)
#output: qx is the density of (y|x) of first order markov chain model
fn_MarkovChain <- function(x,y){
  qx            <- 0
  qx[y==1 & x==1] =  0.8 #true p11 is 0.8
  qx[y==2 & x==1] =  0.2 #true p12 is 0.2
  qx[y==1 & x==2] =  0.2 #true p21 is 0.2
  qx[y==2 & x==2] =  0.8 #true p22 is 0.8
  return(qx)
}  

#to draw 1 observation from the above Markov chain model we defined
#input:  x is  last draw(scalar)
#output: y is new draw(scalar)
fn_draw1obs <- function(x){
  prob          <- fn_MarkovChain(x,x) # probability(new draw=last draw)=0.7
  tmp          <- runif(1)            
  if(tmp<=prob){ y=x }
  else{ y=3-x }
  return(y) #if tmp <= 0.8, new draw is same as the old draw. 
            #Otherwise new draw is different from the old draw
}

#then, draw N observation from the above Markov chain model
#we first initialize the first observation on our own.
#input:  initial value, sample size N
#output: N observation in one vector
fn_drawNobs <- function(Initial,N){
  MCdraws=rep(0,N)      
  MCdraws[1]=Initial
  for(i in 2:N){
    MCdraws[i]=fn_draw1obs(MCdraws[i-1]) }
  return(MCdraws)
}
MC_draws=fn_drawNobs(1,1000) #we simulate 1000 observation with initial value 1

#log likelihood of single pair
#input:  ynew is new draw(scalar)
#        ypast: old draw(scalar)
#        p11,p21
#output: log likelihood 
oneloglikeli <-function(ynew, ypast, p11, p21){
  if(ypast==1){
    l=(p11^(2-ynew))*(1-p11)^(ynew-1)}
  if(ypast!=1 && ypast==2){
    l=(p21^(2-ynew))*(1-p21)^(ynew-1)}
  l=log(l)
  return (l)
}

#log likelihood of n observations
#input:  chain is N draws(vector)
#        p11,p12
#output: log likelihood
loglikeli <-function(chain, p11, p21){
  N=length(chain)
  ini=1
  L=oneloglikeli(chain[1],ini,p11,p21)
  
  for(i in 2:N){L=L+oneloglikeli(chain[i],chain[i-1],p11,p21)}
  return (L)
}

#the density of the flat prior
#input:  p11,p12
#output: density of prior
fn_flatprior <- function(p11,p21){
  p=0
  if(p11>=0 && p11<=1 && p21>=0 && p21<=1)
    p =1
  return (p)
}

#log of posterior density
#input:  n draws(1xn vector)
#output: log of posterior
log_posterior <-function(Mchain, p11, p21){
  L=log(fn_flatprior(p11, p21)) + loglikeli(Mchain, p11, p21)
  return (L)
}

#the candidate is the same as flat prior
#draw one candidate(1 pair of p11 and p21) from uniform distribution
#note that we dont use this in our program
draw1Candi_unif <- function(){
  p <-c(0,0)
  p[1]<- runif(1)            
  p[2]<- runif(1)
  return(p)
}

#draw one candidate(1 pair of p11 and p21) from normal distribution
#output: p is vector of 2 elements (p11, p21)
draw1Candi_normal <-function(meancandi){
  p=c(2,2) #c(2,2) doesn't have any meaning. just make us to run the while loop that
           #initial value of p[1] and p[2] are not probabilities
  while(p[1]<0 || p[1]>1|| p[2]<0|| p[2]>1){
    p[1]=rnorm(1,meancandi[1],1)
    p[2]=rnorm(1,meancandi[2],1)
  }
  return (p)
}

# Function to get `acceptance probability'
# input:  logpx_old [scalar] log of target density values for last draw
#         qx_old [scalar] candidate density values for last draw
#         logpx_new [scalar] log oftarget density values for new draw
#         qx_new [scalar] candidate density values for new draw
# output: alpha  [scalar] acceptance probability
fn_acc_prob <- function(logpx_old,logpx_new,qx_old,qx_new){
  tmp   <- logpx_new-logpx_old+log(qx_old)-log(qx_new)
  logalpha <- min(tmp,log(1))
  alpha=exp(logalpha)
  return(alpha)
}

# MH steps:
n <- 10000  #number of draws

# initialize algorithm
p11drawsaccepted<-rep(0,10000) #storage of vectors for the accepted draws of p11
p21drawsaccepted<-rep(0,10000) #for p21
p11draws<-rep(0,10000) #storage of vectors for the candidate dtwas of p11
p21draws<-rep(0,10000) #for p21

#storage of vectors of alphas from MH algorithm
alphalist <- rep(0,10000)

counter <- 0 # counter for number of accepted draws

#approximate MLE p11 and p21
#find the starting value for p11 and p21
Diff=rep(0,999)

Diff=MC_draws[1:999]-MC_draws[2:1000]

ini21=sum(Diff==1)  #2 jump to 1
ini12=sum(Diff==-1) #1 jump to 2

inip12=ini12/sum(MC_draws==1) #approximate p12
inip11=1-inip12
inip21=ini21/sum(MC_draws==2) #approximate p21

old_draws=rep(0,2)
old_draws[1]=inip11 #starting value for p11
old_draws[2]=inip21 #for p21

inip11
inip21

#for loop
for(i in 1:n){
  # Step 1: draw from the candidate density  
  new_draws = draw1Candi_normal(c(inip11,inip21)) #get one random vector of (p11,p21)
  
  # Step 2: get acceptance rate 
  px_old <- log_posterior(MC_draws, old_draws[1], old_draws[2])
  px_new <- log_posterior(MC_draws, new_draws[1], new_draws[2])
  qx_old <- dnorm(old_draws[1],inip11,1,log=FALSE)*dnorm(old_draws[2],inip21,1,log=FALSE)
  qx_new <- dnorm(new_draws[1],inip11,1,log=FALSE)*dnorm(new_draws[2],inip21,1,log=FALSE)
  
  # get probability of alpha
  alpha  <- fn_acc_prob(px_old,px_new,qx_old,qx_new)
  ###record alpha
  alphalist[i]<-alpha
  
  # Step 3: accept draw with probability alpha
  if(runif(1)<=alpha){
    counter = counter + 1         # increase count of accepted draws
    old_draws   = new_draws       # accept new draw
  }
  
  # Step 4: store draws (repeated or accepted)
  p11drawsaccepted[i] = old_draws[1]
  p11draws[i]=new_draws[1]          #store all the draws p11 from candidate density
  p21drawsaccepted[i] = old_draws[2]
  p21draws[i]=new_draws[2]          #store all the draws p21 from candidate density
}

counter

##burn in 
p11afterburnin<-rep(0,8000) #make a new vector to store the draws after burn in
p21afterburnin<-rep(0,8000)

for(i in 1:8000){
  p11afterburnin[i]=p11drawsaccepted[2000+i]
  p21afterburnin[i]=p21drawsaccepted[2000+i]
}

#trim
p11trimed<-rep(0,2000)     #make a new vector to store the draws after trim
p21trimed<-rep(0,2000)

for(i in 1:2000){
  p11trimed[i]=p11afterburnin[4*i] #keep 1 for each 5 draws
  p21trimed[i]=p21afterburnin[4*i]
}

#to check autocorrelation before/after trimming
acf(p11afterburnin)
acf(p11trimed)
acf(p21afterburnin)
acf(p21trimed)

#########
mean(p11drawsaccepted)
mean(p21drawsaccepted)

mean(p11trimed) #true value 0.8
var(p11trimed)
quantile(p11trimed,seq(0,1,0.05))

mean(p21trimed) #true value 0.2
var(p21trimed)
quantile(p21trimed,seq(0,1,0.05))

p12trimed=1-p11trimed
mean(p12trimed) #true value 0.2
var(p12trimed)
quantile(p12trimed,seq(0,1,0.05))

p22trimed=1-p21trimed
mean(p22trimed) #true value 0.8
var(p22trimed)
quantile(p22trimed,seq(0,1,0.05))

#calculate w1 and w2
w1=p11trimed/p21trimed
w2=p12trimed/p22trimed
mean(w1)
var(w1)
quantile(w1,seq(0,1,0.05)) #true value 0.8/0.2


########
plot(1:1000,MC_draws,type="l")
plot(1:10000, p11drawsaccepted)
plot(1:8000, p11afterburnin)
plot(1:2000, p11trimed)
plot(1:10000, p21drawsaccepted)
plot(1:8000, p21afterburnin)
plot(1:2000, p21trimed)
plot(seq(0,1,0.05),quantile(p11trimed,seq(0,1,0.05)))
plot(seq(0,1,0.05),quantile(p21trimed,seq(0,1,0.05)))
plot(seq(0,1,0.05),quantile(p12trimed,seq(0,1,0.05)))
plot(seq(0,1,0.05),quantile(p22trimed,seq(0,1,0.05)))
plot(1:2000, w1)
plot(1:2000, w2)

# PRINT SUMMARY 
cat('Metropolis Hasting algorithm',fill=TRUE)

cat('True value of p11:', '0.8',fill=TRUE)
cat('Mean of MH p11',mean(p11trimed),fill=TRUE)
cat('Variance of MH p11',var(p11trimed),fill=TRUE)

cat('True value of p12:', '0.2',fill=TRUE)
cat('Mean of MH p12',mean(p12trimed),fill=TRUE)
cat('Variance of MH p12',var(p12trimed),fill=TRUE)

cat('True value of p21:', '0.2',fill=TRUE)
cat('Mean of MH p21',mean(p21trimed),fill=TRUE)
cat('Variance of MH p21',var(p21trimed),fill=TRUE)

cat('True value of p22:', '0.8',fill=TRUE)
cat('Mean of MH p22',mean(p22trimed),fill=TRUE)
cat('Variance of MH p22',var(p22trimed),fill=TRUE)

cat('True value of w1 (p11/p21):', 0.8/0.2,fill=TRUE)
cat('Mean of MH w1',mean(w1),fill=TRUE)
cat('Variance of MH w1',var(w1),fill=TRUE)
cat('Credible interval of MH w1',quantile(w1,c(0.05,0.95)),fill=TRUE)

cat('True value of w2 (p12/p22):', 0.2/0.8,fill=TRUE)
cat('Mean of MH w2',mean(w2),fill=TRUE)
cat('Variance of MH w2',var(w2),fill=TRUE) 
cat('Credible interval of MH w2',quantile(w2,c(0.05,0.95)),fill=TRUE)

cat('Total simulation from MH:',n,fill=TRUE)
cat('The number of draws accepted in the application of MH algorithm',counter,fill=TRUE)
cat('acceptance rate:',counter/n,fill=TRUE)

