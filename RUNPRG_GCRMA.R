


preprocessing<-function(celfiles,sample1,sample2)

{
fn.batch<-ReadAffy(filenames=celfiles)


#look at all probesets with 11 probes each
ProbeNames <- probeNames(fn.batch)
wheatprobeset<-names(which(table(ProbeNames)==11))


ai <- compute.affinities.local(fn.batch, Array=NULL)

#ai2 <- compute.affinities(cdfName(fn.batch))


fn.batch.trans.GCRMA<- bg.adjust.gcrma(fn.batch, affinity.info = ai, type="fullmodel", fast = FALSE)
fn.batch.trans.GCRMA.exprs<-rma(fn.batch.trans.GCRMA, subset=wheatprobeset,background=FALSE)

fn.batch.trans.GCRMA.pm<-pm(fn.batch.trans.GCRMA,wheatprobeset)
fn.batch.trans.GCRMA.names<-probeNames(fn.batch.trans.GCRMA,wheatprobeset)
fn.batch.trans.GCRMA.probeset<-exprs(fn.batch.trans.GCRMA.exprs)
fn.batch.trans.GCRMA.probesetName<-featureNames(fn.batch.trans.GCRMA.exprs)

#  find probesets that are very likely to have absent target transcripts
#  probesets that are very likely to have an absent target transcript
#  are defined as having MAS5 present/absent P values for all replicates > 0.50 

raw.calls <- mas5calls.AffyBatch(fn.batch)

f1 <- kOverA(length(sampleNames(fn.batch)), 0.5)
ff1 <- filterfun(f1)

ab.50 <- genefilter(assayData(raw.calls)[["se.exprs"]], ff1)
gn <- rownames(exprs(raw.calls))
ab.50 <- gn[ab.50]

#  replace MM probe values with a PM threshold value

for(i in 1:length(sampleNames(fn.batch))){
mm(fn.batch.trans.GCRMA)[,i] <- mean(pm(fn.batch.trans.GCRMA, ab.50)[,i], trim=0.02)}


#  calculate present/absent calls
#  alpha values set the "present" and "marginal" P value cutoffs
#  alpha1 is set to 0.03 and the expected false discover rate at this P value is roughly 5% to 8%
#  alpha1 is set to 0.12 and the expected false discover rate at this P value is roughly 10% to 13%

thres.calls <- mas5calls(fn.batch.trans.GCRMA, tau=0.015, alpha1 = 0.024, alpha2 = 0.111 )


#Affy's text files reporting different 'present','Absent', or 'Marginal' calls.
N<-length(exprs(thres.calls)[1,])
n<-length(fn.batch.trans.GCRMA.probesetName)
fnCALL<-matrix(NA,ncol=(N+1),nrow=n) #matrix containing probe names and Affy's calls
fnCALL[,1]<-fn.batch.trans.GCRMA.probesetName

fnCALL[,2:(N+1)]<-exprs(thres.calls)[match(wheatprobeset,featureNames(fn.batch)),]

            
fnPMlog<-as.list(rep(NA,4))
for (i in 1:4)
{fnPMlog[[i]]<-matrix(log(as.numeric(fn.batch.trans.GCRMA.pm[,i]),2),ncol=11,byrow=T)
 medianA<-fn.batch.trans.GCRMA.probeset[,i]
 fnPMlog[[i]]<-fnPMlog[[i]]-medianA
} 


workdata<-datascreen.XC(fnPMlog, NumA=2,NumB=2,fnCALL)

outlyingoutput<-compute.uscores.XC(workdata)
uname<-outlyingoutput[[1]]
uscore<-outlyingoutput[[2]]
w1score<-outlyingoutput[[3]]
w2score<-outlyingoutput[[4]]
probesetnames<-wheatprobeset
affinity<-fnPMlog
SFPselection(0.15,uname,uscore,w2score,file1,file2,fn.batch,affinity,probesetnames)
          
SFPselection(0.10,uname,uscore,w2score,file1,file2,fn.batch,affinity,probesetnames)
          
SFPselection(0.05,uname,uscore,w2score,file1,file2,fn.batch,affinity,probesetnames)
          
}






