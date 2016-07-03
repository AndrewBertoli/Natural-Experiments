External_Validity=function(Sample,Population,Covs,Names,Title="External Validity Graph",Sample_Color="cornflowerblue",Population_Color="lightgray",XLab="",YLab="Value",ln=c()){
	
Index=length(Covs):1	

Covs=rev(Covs)
Names=rev(Names)

Data1=Sample[,Covs]	
Data2=Population[,Covs]	



for(i in ln){Data1[,ln]=log(Data1[,ln]+1);Data2[,ln]=log(Data2[,ln]+1)}


if(nrow(Sample)>nrow(Population)){print("Warning: Popoulation size is smaller than sample size.")
Fill=matrix(NA,nrow=nrow(Data1)-nrow(Data2),ncol=ncol(Data2))

colnames(Fill)=colnames(Data2)

Data1=rbind(Data2,Fill)

}

if(nrow(Sample)<nrow(Population)){Fill=matrix(NA,nrow=nrow(Data2)-nrow(Data1),ncol=ncol(Data1))
	
colnames(Fill)=colnames(Data1)

Data1=rbind(Data1,Fill)

}	
	
Merged=cbind(Data1[,1],Data2[,1])
for(i in 2:length(Covs)){Merged=cbind(Merged,Data1[,i],Data2[,i])}

Names2=rep(c("Sample","Population"),length(Covs))
	
colnames(Merged)=paste(Names2,rep(Names,each=2))

Merged=melt(Merged)

colnames(Merged)[2:3]=c("Variable","Value")	

Merged=Merged[-which(is.na(Merged$Value)==TRUE),]

Merged$Variable=factor(Merged$Variable,levels=paste(Names2,rep(Names,each=2)),ordered=TRUE)	
	
p2 = ggplot(Merged, aes(Variable,Value)) + geom_boxplot(fill=rep(c(Sample_Color,Population_Color),length(Covs))) + coord_flip() + ylab(YLab)+ xlab(XLab)  + theme_bw() +theme(axis.title=element_text(size=16)) + ggtitle(Title) +theme(plot.title = element_text(lineheight=1.8,size=rel(1.5),face="bold"))	

return(p2)	
	
}
