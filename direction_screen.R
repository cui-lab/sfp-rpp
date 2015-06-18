

direction.screen<-function(Y)
{

## projection direction screen, based on function SNP3_2nd
## detection methods include
## 1) Hampel identifier
## 2) Mahalanobis distance
## 3) Maximum within variance

n<-length(Y[,1])
p<-length(Y[1,])   #n by p

# - - (1 initial computation

alpha<-0.05 # percentile of tail

mean.Y<-apply(Y,2,mean) #p by 1
s<-cov(Y)               #p by p, sample cov.matrix
sinv<-solve(s)          #inverse of s, p by p

# -- (2) compute Hampel identifier

#pY=[];                 # prejected values of Y
#med_pY=[];             % medians of pY
#mad_pY=[];             % median of absulate deviation of pY
#ham_pY=[];             % Hampel identifier of Y                       
#pos_ham=[];            % sig. pos. by Hampel iden.



ch.v<-eigen(s,symmetric=T)$vectors # ch. vectors and roots of
ch.r<-eigen(s,symmetric=T)$values  # the sample cov. matrix
                       
 

 pY<-Y%*%ch.v[,1:3]    # take first three eigenvectors as projection directions
                       # Y projected on i^st ch.v.n by 3

 med.pY<-apply(pY,2,median)   # median of pY(i), 1 by 3
 mad.pY<-apply(abs(pY-matrix(rep(med.pY,n),ncol=3,byrow=T)),2,median)

                              # medians of abs. deviation of pY, 1 by 3

ham.pY<-matrix(NA,ncol=3,nrow=n)
for (i in 1:3)
  ham.pY[,i]<-(pY[,i]-med.pY[i])/mad.pY[i];
                       # Hampel identifier, n by 3




perc.h<-apply(ham.pY, 2,function(x)quantile(x,1-alpha)) 
                       #95 percentiles of Hampel

pos.ham<-c()                       # identifiers, 1 by 3
for (i in 1:3)
    {
    
    m.pos<-c(1:n)[ham.pY[,i]>perc.h[i]]
                       # locate significant pos.
    pos.ham<-c(pos.ham, m.pos)
    }
                               # take union of sig. pos.
pos.ham<-unique(pos.ham)       # from three directions
                       


# -- (3) compute Mahalanobis identifier
MD<-rep(NA,n)

for(i in 1:n)
MD[i]<-t(Y[i,]-mean.Y)%*%sinv%*%((Y[i,]-mean.Y))

                       # Mahalanobis distances
                       # 1 by n

perc.MD<-quantile(MD,1-alpha)
                       # 95 percentile of MD
pos.MD=c(1:n)[MD>perc.MD];
                       # 95 sig. pos. by MD 
                       
# -- (5) Maximum with-in variance method

var.Y<-apply(Y,1,var)  #with-in variance of Y
                       # 1 by n
perc.v<-quantile(var.Y,1-alpha);
                       # 95 percentile of with-in 
                       # variances 
pos.v<-c(1:n)[var.Y>perc.v]
                       # 95 sig. pos. by maximum var.

# -- (6) combine all the sig. positions

pos<-c(pos.ham,pos.MD,pos.v)
pos<-unique(pos)
                       # take union of three
                       # sig. positions

return(pos)

}
  
