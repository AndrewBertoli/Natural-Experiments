# This R function estimates the p-value for an experiment where every unit has an equal probability of treatment
# assignment and the size treatment group is fixed. Users simply enter the outcomes for the treatment and control groups.
# If the treatment and control outcomes are the same length, users can specify whether the data is paired or unpaired.
# This function is meant to be used in conjunction with BalancePlot()

PermutationTest=function(Treatment, Control, paired=FALSE, simulations=100000, na.rm=FALSE, output="full"){

form_note<-0   
if(class(Treatment)=="formula"){
form_note<-1 
Treat=as.character(Treatment)[3]
Left=paste("+",as.character(Treatment)[2])
Left=strsplit(Left," ")[[1]]
plus=which(Left=="+")
minus=which(Left=="-")
variables_index=(1:length(Left))[-c(plus,minus)]
outcome=0
for(i in variables_index){
outcome=Control[,Left[i]]*(2*as.numeric((i-1)%in%plus)-1)+outcome}
Condition1<-sort(unique(Control[,Treat]))[1]
Condition2<-sort(unique(Control[,Treat]))[2]  
Treatment=outcome[Control[,Treat]==Condition2]
Control=outcome[Control[,Treat]==Condition1]}

if(na.rm==FALSE){if(any(is.na(c(Treatment,Control)))==TRUE){return("NAs Detected. Must set na.rm=TRUE")}}



if(paired==TRUE){
differences=Treatment-Control
differences=differences[!is.na(differences)]
new.t.stats=rep(0,simulations)
for(i in 1:simulations){
assignment=sample(c(-1,1),length(differences), replace=TRUE)
new.t.stats[i]=mean(assignment*differences)
}
pvalue=length(which(abs(new.t.stats)>=abs(mean(differences))))/simulations
}

if(paired==FALSE){
Treatment=Treatment[!is.na(Treatment)]
Control=Control[!is.na(Control)]
new.t.stats=rep(0,simulations)
for(i in 1:simulations){
assignment=sample(c(rep(0,length(Control)),rep(1,length(Treatment))),length(c(Treatment,Control)), replace=FALSE)
new.t=c(Treatment,Control)[assignment==1]
new.c=c(Treatment,Control)[assignment==0]
new.t.stats[i]=mean(new.t)-mean(new.c)
}
pvalue=length(which(abs(new.t.stats)>=abs(mean(Treatment)-mean(Control))))/simulations
}

est=mean(Treatment)-mean(Control)
se=sd(new.t.stats)
if(output=="p") return(pvalue)
if(output=="short") return(c(est,se,pvalue))  
n=length(c(Treatment,Control))

if(output=="full"&form_note==0) return(cbind(paste("Group 1 Mean=",mean(Treatment),sep=""),paste("Group 2 Mean=",mean(Control),sep=""),paste("Estimate=",est,sep=""),paste("p-value=",pvalue,sep=""),paste("SE=",se,sep=""),paste("n=",n,sep="")))
if(output=="full"&form_note==1) return(cbind(paste("Group 1: ",Treat,"=",Condition2,sep=""),paste("Group 2: ",Treat,"=",Condition1,sep=""),paste("Group 1 Mean=",mean(Treatment),sep=""),paste("Group 2 Mean=",mean(Control),sep=""),paste("Estimate=",est,sep=""),paste("p-value=",pvalue,sep=""),paste("SE=",se,sep=""),paste("n=",n,sep="")))}
