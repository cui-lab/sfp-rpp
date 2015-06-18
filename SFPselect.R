SFPselection<-function(threshold,uname,uscore,w2score,name1,name2,fnpm,affinity,probesetnames)
{



cutoff<-quantile(uscore,1-threshold)
SFPsetid<-c(1:length(uscore))[uscore>cutoff]

SFPset<-uname[SFPsetid]
SFPset<-cbind(SFPset,uscore[SFPsetid])

nonSFPset<-uname[-SFPsetid]

SFPprobes<-w2score[SFPsetid,]
nonSFPprobes<-w2score[-SFPsetid,]


reorder<-order(as.numeric(SFPset[,2]),decreasing=TRUE)
SFPset<-SFPset[reorder,]
SFPprobes<-SFPprobes[reorder,]



INDpercent<-abs(SFPprobes)/apply(abs(SFPprobes),1,sum)
maxINDpercent<-apply(INDpercent,1,max)


SFPprobeCall<-ifelse(INDpercent>0.4,1,0)



SFPprobeCallName<-c()
SFPprobeCallNameScore<-c()
SFPprobesID<-as.list(rep(NA,length(SFPsetid)))
for(i in 1:length(SFPsetid))
   {if(sum(SFPprobeCall[i,])==0)
     { j<-c(1:11)[INDpercent[i,]==maxINDpercent[i]]
       SFPprobeCallName<-c(SFPprobeCallName,paste(SFPset[i],j,sep=""))
       SFPprobeCallNameScore<-c(SFPprobeCallNameScore,SFPprobes[i,j])
       SFPprobesID[[i]]<-c(SFPprobesID[[i]],j)
     }
     else
     {
      for(j in 1:11)
      if(SFPprobeCall[i,j]==1)
      {SFPprobeCallName<-c(SFPprobeCallName, paste(SFPset[i],j,sep=""))
      SFPprobeCallNameScore<-c(SFPprobeCallNameScore,SFPprobes[i,j])
      SFPprobesID[[i]]<-c(SFPprobesID[[i]],j)
      }
     }
   }

write.table(SFPset,paste(name1,name2,threshold,"SFPprobesets.csv",sep="_"),sep=",",row.names=FALSE,col.names=FALSE)
b<-cbind(SFPprobeCallName,SFPprobeCallNameScore)
write.table(b,paste(name1,name2,threshold,"SFPprobes.csv",sep="_"),sep=",",row.names=FALSE,col.names=FALSE)

save(SFPset,b,file=paste(name1,name2,threshold,"result.RData",sep="_"))

pdf(paste("plot_",threshold,"_",name1,"vs",name2,".pdf"))    
  
par(mfrow = c(2,2))
for (id in 1:length(SFPsetid))
{index<-match(SFPset[id,1],uname)
indscore<-w2score[index,]


#(2) estimate affinity for each chip in group A and group B

NumA<-2

NumB<-2


#fnPMlog<-as.list(rep(NA,NumA+NumB))

#for (i in 1:(NumA+NumB))
#      {fnPMlog[[i]]<-affinity[[cols[which(!is.na(cols))][i]]]
#       fnPMlog[[i]]<-fnPMlog[[i]][match(uname[index],probesetnames),]
#       }


SFPmatrix<-matrix(NA,ncol=11,nrow=NumA+NumB)

for(i in 1:(NumA+NumB)) 
SFPmatrix[i,]<-affinity[[i]][match(uname[index],probesetnames),]


SFPmatrix.lograw<-matrix(NA,ncol=11,nrow=NumA+NumB)

for(i in 1:(NumA+NumB)) 
SFPmatrix.lograw[i,]<-log(pm(fnpm,probesetnames)[probeNames(fnpm,probesetnames)==uname[index],i],2)


#par(mar = c(2.5, 3.5, 1, 0.5))

#par(mgp = c(1.5, 0.5, 0))
#par(oma = c(0, 0, , 0))

matplot(c(1:11),seq(0,15.5,1.5), type= "n", ylab ="Log2(transformed Intensity)", xlab = "Probe Number",main=uname[index],axes=F,frame=T)
matlines(t(SFPmatrix.lograw), col=c(1,1,1,2,2,2))
legend(8,2,c(name1,name2),col=c(1,2),lty=1,cex=0.5)
axis(1,1:11,c(0,1,2,3,4,5,6,7,8,9,10))
axis(2,seq(0,15.5,1.5),seq(0,15.5,1.5))


matplot(c(1:11),seq(-13,13,2.6), type= "n", ylab ="Affinity", xlab = "Probe Number",main=uname[index],axes=F,frame=T)
matlines(t(SFPmatrix),col=c(1,1,1,2,2,2))
legend(8,8,c(name1,name2),col=c(1,2),lty=1,cex=0.5)
axis(1,1:11,c(0,1,2,3,4,5,6,7,8,9,10))
axis(2,seq(-13,13,2.6),seq(-13,13,2.6))


#(4) compute working data, obtained by choosing all "present" call probesets and estimating affinity values
#####estimate affinity for group A and group B using median

SFPmatrixDiff<-SFPmatrix[1,]  
   		for(j in 1:11)
                {tmp1<-c()
                 for(k in 1:NumA)
                 tmp1<-cbind(tmp1,SFPmatrix[k,j])
                  tmp2<-c()
                 for(k in (NumA+1):(NumA+NumB))
                 tmp2<-cbind(tmp2,SFPmatrix[k,j])

                 SFPmatrixDiff[j]<-median(tmp1)-median(tmp2)
          
                
                }


   plot(c(1:11),seq(-10.5,10.5,2), type= "n", xlab="", ylab ="Affinity difference", axes=F, frame=T)
   lines(c(1:11), SFPmatrixDiff)
   axis(1,1:11,c(0,1,2,3,4,5,6,7,8,9,10))
    axis(2,seq(-10.5,10.5,2),seq(-10.5,10.5,2))



plot(c(1:11),indscore, type= "h", col=1, xlab="", ylab ="Indi.outlying score", axes=F, frame=T)
axis(1,1:11,c(0,1,2,3,4,5,6,7,8,9,10))
axis(2,seq(-1,30,1),seq(-1,30,1))

cat(id)

}

dev.off()

  }

######################################################################################################



