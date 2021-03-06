# This R function estimates the p-value for an experiment where every unit has an equal probability of treatment
# assignment and the size treatment group is fixed. Users simply enter the outcomes for the treatment and control groups.
# If the treatment and control outcomes are the same length, users can specify whether the data is paired or unpaired.
# This function is meant to be used in conjunction with BalancePlot()

PermutationTest=function(Treatment, Control, Paired=FALSE, Simulations=500000, na.rm=FALSE, Output="p"){

if(class(Treatment)=="formula"){
Treat=as.character(Treatment)[3]
Left=paste("+",as.character(Treatment)[2])
Left=strsplit(Left," ")[[1]]
plus=which(Left=="+")
minus=which(Left=="-")
variables_index=(1:length(Left))[-c(plus,minus)]
outcome=0
for(i in variables_index){
outcome=Control[,Left[i]]*(2*as.numeric((i-1)%in%plus)-1)+outcome}
Treatment=outcome[Control[,Treat]==1]
Control=outcome[Control[,Treat]==0]}

if(na.rm==FALSE){if(any(is.na(c(Treatment,Control)))==TRUE){return("NAs Detected. Must set na.rm=TRUE")}}



if(Paired==TRUE){
differences=Treatment-Control
differences=differences[!is.na(differences)]
new.t.stats=rep(0,Simulations)
for(i in 1:Simulations){
assignment=sample(c(-1,1),length(differences), replace=TRUE)
new.t.stats[i]=mean(assignment*differences)
}
pvalue=length(which(abs(new.t.stats)>=abs(mean(differences))))/Simulations
}

if(Paired==FALSE){
treatmeant=Treatment[!is.na(Treatment)]
Control=Control[!is.na(Control)]
new.t.stats=rep(0,Simulations)
for(i in 1:Simulations){
assignment=sample(c(rep(0,length(Control)),rep(1,length(Treatment))),length(c(Treatment,Control)), replace=FALSE)
new.t=c(Treatment,Control)[assignment==1]
new.c=c(Treatment,Control)[assignment==0]
new.t.stats[i]=mean(new.t)-mean(new.c)
}
pvalue=length(which(abs(new.t.stats)>=abs(mean(Treatment)-mean(Control))))/Simulations
}

est=mean(Treatment)-mean(Control)
se=sd(new.t.stats)
if(Output=="p") return(pvalue)

if(Output=="full") return(cbind(paste("Estimate=",est,colapse=""),paste("Two-tailed p-value=",pvalue,colapse=""),paste("SE=",se,colapse="")))}
