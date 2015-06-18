


compute.uscores.XC<-function(workdata)
{
  #compute outlying scores u,u_factor1,and u_factor2
  #combine function SNP3_2nd and function SNP8_5th

  n<-length(workdata[,1])
  p<-length(workdata[1,])-1   #n by p

 #---(1) initial computation
 Y<-matrix(as.numeric(workdata[,-1]),ncol=11,byrow=F) 
 
 row.median.Y<-apply(Y,1,median)  # n by 1 row median of Y

 dir.pos<-direction.screen(Y)   #directions

 L<-length(dir.pos)


  # ---(2) project pursuit overall outlying scores


 proj.dire<-Y[dir.pos,]-matrix(1,ncol=1,nrow=L)%*%apply(Y,2,median)
					   # projection directions
                                           #L by p 
 
 proj.Y<-Y%*%t(proj.dire)

 
 
 
 MAD<-rep(NA,L)
   
 med.proj.Y<-c()

 for(i in 1:L)
  {
   medPROJY<-median(proj.Y[,i])
   med.proj.Y<-c(med.proj.Y,medPROJY)
   abs.devi<-abs(proj.Y[,i]-medPROJY)
   MAD[i]<-median(abs.devi)
 }
  
 u<-rep(NA,n)
 max.dire<-rep(NA,n)

 for(i in 1:n)
 {abs.devi<-abs(proj.Y[i,]-med.proj.Y)/MAD
  u[i]<-max(abs.devi)
  max.dire[i]<-c(1:L)[!is.na(match(abs.devi,u[i]))]
 }
 rm(proj.Y)
 gc()
  


        #u: n by 1
        #max.dire: n by 1 record the direction
        #on which maximum outlying score is obtained
        #end of u score computation 

rm(abs.devi)
gc()

# -- (3) project pursuit individual outlying scores
 
u.factor1<-matrix(0,ncol=p,nrow=n)
u.factor2<-matrix(0,ncol=p,nrow=n)

# --first, compute the second portion of u_factor

for (j in 1:p)
    {    
         { 
         if(j==1)
         proj.Y.j<-cbind(row.median.Y,Y[,(j+1):p])%*%t(proj.dire)
         else if(j==11)
          proj.Y.j<-cbind(Y[,1:(j-1)],row.median.Y)%*%t(proj.dire)
              else
              proj.Y.j<-cbind(Y[,1:(j-1)],row.median.Y,Y[,(j+1):p])%*%t(proj.dire)
                              # new projected values
                              # after replacing jth component
                              # by the rwo median
                              # n by L matrix  
     
      for( i in 1:L) 
      proj.Y.j[,i]<-abs(proj.Y.j[,i]-med.proj.Y[i])/MAD[i]

      u.factor1[,j]<-apply(proj.Y.j,1,max)
                              # n by 1
                              # first way to compute u.factor
                              # saved in u.factor1





        if(j==1)
       
        
         proj.Y.j<-cbind(row.median.Y,Y[,(j+1):p])
        else if(j==11)
         
         proj.Y.j<-cbind(Y[,1:(j-1)],row.median.Y)
              else 

              proj.Y.j<-cbind(Y[,1:(j-1)],row.median.Y,Y[,(j+1):p])

        tmp<-proj.Y.j[,1]
        

        for(k in 1:n)
        {tmp[k]<-proj.Y.j[k,]%*%proj.dire[max.dire[k],]
        u.factor2[k,j]<-abs(tmp[k]-med.proj.Y[max.dire[k]])/MAD[max.dire[k]]
        }
      } 
     cat(j)
     gc()
     }

   
# - - second, compute the u.facotor


u.factor1<- u%*%matrix(1,ncol=11,nrow=1)-u.factor1
u.factor2<- u%*%matrix(1,ncol=11,nrow=1)-u.factor2

a<-list(workdata[,1],u,u.factor1,u.factor2)

return(a)
}