# The following function will compute the p-values for different Gamma levels for paired or unpaired data. 
# The argument "Gamma" is the Gamma level, defined as the maximum difference in treatment odds between any 
# two units in the sample. In general, setting Gamma=j will return the same p-value as Gamma=1/j. When "Type" 
# is set at "mean", the function will do a permutation test with the mean as the test statistic. When "Type" is
# set at "rank", the function will do the Wilcoxon Signed-Rank Test or Rank Sum Tests, which are less sensitive 
# to outliers. The p-values are two-tailed and computed using Monte Carlo simulation, so the function will take 
# some time to run. This method differs from the functions in the "rbound" package, which rely on the normal 
# approximation and do not allow the user to use the mean as the test statistic.

require(Rlab)

Sens=function(Treat.Outcomes,Control.Outcomes,Gamma=1,Type="Mean",Paired=TRUE,na.rm=FALSE,Simulations=1000000){

if(na.rm==FALSE){if(any(is.na(c(Treat.Outcomes,Control.Outcomes)))==TRUE){return("NAs Detected. Must set na.rm=TRUE")}}

if(na.rm==TRUE){
Treat.Outcomes <- Treat.Outcomes[!is.na(Treat.Outcomes)]
Control.Outcomes <- Control.Outcomes[!is.na(Control.Outcomes)]
}

if(Paired==TRUE){

differences=Treat.Outcomes-Control.Outcomes

differences=differences[differences!=0]    

if(Type=="Mean"){
real.t.stat=mean(differences)
fake.t.stats=rep(0,Simulations)
for(i in 1:Simulations){
treatment.assignment.t.higher=sample(c(1,-1),size=sum(differences>0),replace=TRUE,prob=c(Gamma/(Gamma+1),1/(Gamma+1)))    
treatment.assignment.c.higher=sample(c(1,-1),size=sum(differences<0),replace=TRUE,prob=c(1/(Gamma+1),Gamma/(Gamma+1)))    
first.group=differences[differences>0]*treatment.assignment.t.higher
second.group=differences[differences<0]*treatment.assignment.c.higher
new.differences=c(first.group,second.group)
fake.t.stats[i]=mean(new.differences)
}
p.val=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)
return(signif(p.val,2))}

if(Type=="Rank"){
ranks=rank(abs(differences))
real.t.stat=sum(ranks[differences>0])-sum(ranks)/2    
fake.t.stats=rep(0,Simulations)
for(i in 1:Simulations){
treatment.assignment.t.higher=sample(c(1,-1),size=sum(differences>0),replace=TRUE,prob=c(Gamma/(Gamma+1),1/(Gamma+1)))    
treatment.assignment.c.higher=sample(c(1,-1),size=sum(differences<0),replace=TRUE,prob=c(1/(Gamma+1),Gamma/(Gamma+1)))    
first.group=differences[differences>0]*treatment.assignment.t.higher
second.group=differences[differences<0]*treatment.assignment.c.higher
new.differences=c(first.group,second.group)
new.ranks=rank(abs(new.differences))
fake.t.stats[i]=sum(new.ranks[new.differences>0])-sum(new.ranks)/2
}
p.val=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)
return(signif(p.val,2))}

}

if(Paired==FALSE){

if(Type=="Mean"){
real.t.stat=mean(Treat.Outcomes)-mean(Control.Outcomes)
outcomes=c(Treat.Outcomes,Control.Outcomes)[order(c(Treat.Outcomes,Control.Outcomes))]
p.vals=rep(0,length(outcomes))
if(length(p.vals)>50){
for(k in 1:length(outcomes)){
fake.t.stats=rep(0,1000)
for(i in 1:1000){
treatment.assignment=rbern(length(outcomes),c(rep(Gamma/(Gamma+1),k),rep(1/(Gamma+1),length(outcomes)-k)))
fake.t.stats[i]=mean(outcomes[treatment.assignment==1])-mean(outcomes[treatment.assignment==0])
}
p.vals[k]=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)}
candidate.ks=tail(sort(p.vals),50)}

if(length(p.vals)<=50) candidate.ks=1:length(outcomes)

for(k in 1:length(candidate.ks)){
fake.t.stats=rep(0,2000)
for(i in 1:2000){
treatment.assignment=rbern(length(outcomes),c(rep(Gamma/(Gamma+1),candidate.ks[k]),rep(1/(Gamma+1),length(outcomes)-candidate.ks[k])))
fake.t.stats[i]=mean(outcomes[treatment.assignment==1])-mean(outcomes[treatment.assignment==0])
}
p.vals[k]=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)}

candidate.ks=tail(sort(p.vals),10)
p.vals=rep(0,10)
for(k in 1:10){
fake.t.stats=rep(0,20000)
for(i in 1:20000){
treatment.assignment=rbern(length(outcomes),c(rep(Gamma/(Gamma+1),candidate.ks[k]),rep(1/(Gamma+1),length(outcomes)-candidate.ks[k])))
fake.t.stats[i]=mean(outcomes[treatment.assignment==1])-mean(outcomes[treatment.assignment==0])
}
p.vals[k]=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)}
k=candidate.ks[which.max(p.vals)]
fake.t.stats=rep(0,Simulations)
for(i in 1:Simulations){
treatment.assignment=rbern(length(outcomes),c(rep(Gamma/(Gamma+1),k),rep(1/(Gamma+1),length(outcomes)-k)))
fake.t.stats[i]=mean(outcomes[treatment.assignment==1])-mean(outcomes[treatment.assignment==0])
}
p.value=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)
return(signif(p.value,2))}

if(Type=="Rank"){

ranks=rank(c(Treat.Outcomes,Control.Outcomes))
real.t.stat=mean(ranks[1:length(Treat.Outcomes)])-mean(ranks[(length(Treat.Outcomes)+ 1):length(c(Treat.Outcomes,Control.Outcomes))])
ranks=ranks[order(ranks)]

p.vals=rep(0,length(ranks))
if(length(p.vals)>50){
for(k in 1:length(ranks)){
fake.t.stats=rep(0,1000)
for(i in 1:1000){
treatment.assignment=rbern(length(ranks),c(rep(Gamma/(Gamma+1),k),rep(1/(Gamma+1),length(ranks)-k)))
fake.t.stats[i]=mean(ranks[treatment.assignment==1])-mean(ranks[treatment.assignment==0])
}
p.vals[k]=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)}
candidate.ks=tail(sort(p.vals),50)}

if(length(p.vals)<=50) candidate.ks=1:length(ranks)

for(k in 1:length(candidate.ks)){
fake.t.stats=rep(0,2000)
for(i in 1:2000){
treatment.assignment=rbern(length(ranks),c(rep(Gamma/(Gamma+1),candidate.ks[k]),rep(1/(Gamma+1),length(ranks)-candidate.ks[k])))
fake.t.stats[i]=mean(ranks[treatment.assignment==1])-mean(ranks[treatment.assignment==0])}    
p.vals[k]=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)}
candidate.ks=tail(sort(p.vals),10)
p.vals=rep(0,10)
for(k in 1:10){
fake.t.stats=rep(0,20000)
for(i in 1:20000){    
treatment.assignment=rbern(length(ranks),c(rep(Gamma/(Gamma+1),candidate.ks[k]),rep(1/(Gamma+1),length(ranks)-candidate.ks[k])))
fake.t.stats[i]=mean(ranks[treatment.assignment==1])-mean(ranks[treatment.assignment==0])}    
p.vals[k]=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)}
k=candidate.ks[which.max(p.vals)]

fake.t.stats=rep(0,Simulations)
for(i in 1:Simulations){
treatment.assignment=rbern(length(ranks),c(rep(Gamma/(Gamma+1),k),rep(1/(Gamma+1),length(ranks)-k)))
fake.t.stats[i]=mean(ranks[treatment.assignment==1])-mean(ranks[treatment.assignment==0])    }
p.value=length(which(abs(fake.t.stats)>=abs(real.t.stat)))/length(fake.t.stats)

return(signif(p.value,2))}}}

# Example

Sens(t$AGGAfter-t$AGGBefore,c$AGGAfter-c$AGGBefore, Gamma=1.5, Paired=TRUE, Type="Mean")
