
datascreen.XC<-function(fnPM,NumA,NumB,CALL)
{
#probe sets preparation, including
#(1) select probe sets with "present" calls in all chips by MAS5
#(2) use XC model to compute the affinity level for each probe set
#(3) working data for SFP detection is obtained by taking the difference of the average of affinity value (evaluated by sample median)

## (1) select "present"-called-in-all-chips probe sets

# read in probe set names from Affy TXT CALL file
# including some Affy control probe set names

# subset probeset with present calls
MAS.Name<-CALL[,1]
CALLP<-apply(ifelse(CALL[,-1]=="P",1,0),1,sum)
presentName<-MAS.Name[CALLP==(NumA+NumB)]


#(2) estimate affinity for each chip in group A and group B


for (i in 1:(NumA+NumB))

 fnPM[[i]]<-fnPM[[i]][CALLP==(NumA+NumB),]


#(4) compute working data, obtained by choosing all "present" call probesets and estimating affinity values
#####estimate affinity for group A and group B using median
N<-length(fnPM[[1]][,1])
affiA<-fnPM[[1]]
affiB<-fnPM[[1]]

         
   		for(j in 1:11)
                {tmp<-c()
                 for(k in 1:NumA)
                 tmp<-cbind(tmp,fnPM[[k]][,j])
                  affiA[,j]<-apply(tmp,1,median)
                 tmp<-c()
                 for(k in (NumA+1):(NumA+NumB))
                 tmp<-cbind(tmp,fnPM[[k]][,j])
                  affiB[,j]<-apply(tmp,1,median)
                
                }
 

 	 

workdata<-affiA-affiB
workdata<-cbind(presentName,workdata)
rm(affiA,affiB,MAS.Name,CALLP,presentName)
gc()
return(workdata)
}