# This function creates regression discontinuity plots with bootstrapped standard errors. 
# It provides users with a number of options, including whether to plot the raw data 
# or binned means, which smoother and kernel to use, and where to set the bandwidth. 
# For bandwidth choice, I recommend using the banwidth selection function rdbwselect() from the rdrobust package. 

RDPlot=function (X, Y, C = 0, xlab = "Forcing Variable",
ylab = "Outcome", Main = "Regression Discontinuity Plot",Title.Size=20,
Plot.Means = TRUE, Plot.Raw.Data = FALSE,Raw.Data.Point.Size=4,Raw.Data.Point.Density=0.6, Mean.Colors = c("blue",
"red"), Raw.Data.Colors = c("lightblue", "pink"),Line.Colors=c("black","black"), Point.Size = 4,
Shade.Color = "gray87", Window = "None", Plot.p.value = TRUE,
Tick.Marks = seq(-.2,.2,by=0.05),Labels=c("-20%","-15%","-10%","-5%","0%","5%","10%","15%","20%"),
xlim = range(Tick.Marks), ylim = "Automatic", Bandwidth = 2, NBoots = "Automatic", Breaks = "Automatic",
Parametric = FALSE,Smoother="Kernel", Kernel=NULL, Poly.Order=NULL, Cluster = NULL, Jitter = FALSE, 
Confidence.Level=0.95, Type="Two-Sided", Loc.p.value = "BR",Absolute.Y.Min=-10000000,Absolute.Y.Max=10000000,...)
{
 if (length(Cluster) != 0 & Parametric == TRUE) {
 if (length(names(table(as.vector(Cluster)))[table(as.vector(Cluster)) ==
 max(table(as.vector(Cluster)))]) != length(unique(Cluster))) {
 return("Clusters Must Be of Equal Size for Parametric Bootstrapping")
 }
 }
 if (length(Breaks) == 1) {
 Breaks = c(seq(min(X) - 1e-05, C, by = (C - min(X) +
 1e-05)/10), seq(max(X)/10, max(X), by = max(X)/10))
 }
 bins = cut(X, breaks = Breaks)
 dat = aggregate(Y ~ bins, FUN = mean)
 midpoints = rep(0, length(Breaks) - 1)
 for (i in 1:length(midpoints)) {
 midpoints[i] = (Breaks[i] + Breaks[i + 1])/2
 }
 midpoints = midpoints[as.numeric(dat[[1]])]
 negative.midpoints = midpoints[midpoints < C]
 positive.midpoints = midpoints[midpoints > C]
 pointsize = rep(0, length(Breaks) - 1)
 for (i in 1:length(pointsize)) {
 pointsize[i] = length(which(X > Breaks[i] & X <= Breaks[i +
 1]))
 }
 pointsize = Point.Size * sqrt(pointsize[as.numeric(dat[[1]])])
 if (length(ylim) == 1 & Plot.Raw.Data == TRUE) {
 ylim = range(Y)
 }
 if (length(ylim) == 1 & Plot.Raw.Data == FALSE) {
 ylim = range(dat[, 2])
 }


if(NBoots=="Automatic"){
if(length(X)>200){NBoots=100}
if(length(X)<=200){NBoots=200}
if(length(X)<150){NBoots=500}
if(length(X)<100){NBoots=1000}
if(length(X)<70){NBoots=2000}
if(length(X)<50){NBoots=3000}
if(length(X)<30){NBoots=50000}
if(length(X)<20){NBoots=10000}


}


 bootstrapmatrix1 = matrix(0, nrow = 100 * NBoots, ncol = 2)
 if (Parametric == FALSE) {
 for (i in 1:NBoots) {
 if (length(Cluster) == 0) {
 randompoints = sample(seq(1, sum(X < C)), sum(X <
 C), replace = TRUE)
 new.X = X[X < C][randompoints]
 new.Y = Y[X < C][randompoints]
 }
 if (length(Cluster) != 0) {
 randomClusters = sample(unique(Cluster[X < C]),
 length(unique(Cluster[X < C])), replace = TRUE)
 new.X = c()
 for (j in 1:length(randomClusters)) {
 new.X = c(new.X, X[X < C][Cluster[X < C] ==
 randomClusters[j]])
 }
 new.Y = c()
 for (j in 1:length(randomClusters)) {
 new.Y = c(new.Y, Y[X < C][Cluster[X < C] ==
 randomClusters[j]])
 }
 }



weight.function=function(v,p){
if(length(Kernel)==0){Kernel="Box"}
if(Kernel=="Box"){return(rep(1,length(v)))}
if(Kernel=="Triangular"){return(Bandwidth-abs(v-p))}
if(Kernel=="Epanechnikov"){return(3/4(1-((v-p)/(Bandwidth))^2))}
if(Kernel=="Quartic"){return(15/16*((1-((v-p)/(Bandwidth))^2)^2))}
if(Kernel=="Triweight"){return(35/32*((1-((v-p)/(Bandwidth))^2)^3))}
if(Kernel=="Tricube"){return(70/81*((1-(abs(v-p)/(Bandwidth))^3)^3))}}


if(Smoother=="Nearest Neighbor"){

Smoother="Kernel"
Kernel="Box"

}

if(Smoother=="Kernel"){

smoother=function(x,y,c,bw){
if(c<mean(x)){range=c(c,max(X))}
if(c>mean(x)){range=c(min(X),c)}
xs=seq(range[1],range[2],by=(range[2]-range[1])/99)
ys=rep(0,length(xs))
for(k in 1:length(xs)){
obsx=x[x<(xs[k]+bw)&x>(xs[k]-bw)]
obsy=y[x<(xs[k]+bw)&x>(xs[k]-bw)]
if(length(obsx)==0){ys[k]=Y[which.min(abs(X-xs[k]))]}
if(length(obsx)>=1){
weights=weight.function(obsx,xs[k])
ys[k]=sum(obsy*weights)/sum(weights)
}}
combined=data.frame(cbind(xs,ys))
colnames(combined)=c("x","y")
return(combined)
}}


if(Smoother=="Local Linear"){

smoother=function(x,y,c,bw){
if(c<mean(x)){range=c(c,max(X))}
if(c>mean(x)){range=c(min(X),c)}
xs=seq(range[1],range[2],by=(range[2]-range[1])/99)
ys=rep(0,length(xs))
for(k in 1:length(xs)){
obsx=x[x<(xs[k]+bw)&x>(xs[k]-bw)]
design=cbind(1,obsx)
obsy=y[x<(xs[k]+bw)&x>(xs[k]-bw)]
if(length(unique(obsx))<=1){ys[k]=Y[which.min(abs(X-xs[k]))]}
if(length(unique(obsx))>1){
weights=weight.function(design[,2],xs[k])
beta=summary(lm(obsy~obsx,weights=weights))$coefficients[1:2,1]
ys[k]=beta[1]+beta[2]*xs[k]}
}
combined=data.frame(cbind(xs,ys))
colnames(combined)=c("x","y")
return(combined)
}}




if(Smoother=="Local Polynomial"){

if(length(Poly.Order)==0){Poly.Order=4}

smoother=function(x,y,c,bw){
if(c<mean(x)){range=c(c,max(X))}
if(c>mean(x)){range=c(min(X),c)}
xs=seq(range[1],range[2],by=(range[2]-range[1])/99)
ys=rep(0,length(xs))
for(k in 1:length(xs)){
obsx=x[x<(xs[k]+bw)&x>(xs[k]-bw)]
if(length(unique(obsx))<=(Poly.Order+1)){ys[k]=Y[which.min(abs(X-xs[k]))]}
if(length(unique(obsx))>(Poly.Order+1)){
design=cbind(1,obsx,obsx^2,obsx^3,obsx^4,obsx^5,obsx^6,obsx^7,obsx^8,obsx^9)
design=design[,1:(Poly.Order+1)]
obsy=y[x<(xs[k]+bw)&x>(xs[k]-bw)]
weights=weight.function(design[,2],xs[k])
if(Poly.Order==1) beta=summary(lm(obsy~obsx,weights=weights))$coefficients[1:2,1]
if(Poly.Order==2) beta=summary(lm(obsy~obsx+I(obsx^2),weights=weights))$coefficients[1:3,1]
if(Poly.Order==3) beta=summary(lm(obsy~obsx+I(obsx^2)+I(obsx^3),weights=weights))$coefficients[1:4,1]
if(Poly.Order==4) beta=summary(lm(obsy~obsx+I(obsx^2)+I(obsx^3)+I(obsx^4),weights=weights))$coefficients[1:5,1]
if(Poly.Order==5) beta=summary(lm(obsy~obsx+I(obsx^2)+I(obsx^3)+I(obsx^4)+I(obsx^5),weights=weights))$coefficients[1:6,1]
vect=c(1,xs[k],xs[k]^2,xs[k]^3,xs[k]^4,xs[k]^5,xs[k]^6,xs[k]^7,xs[k]^8,xs[k]^9)[1:(Poly.Order+1)]
ys[k]=sum(beta*vect)}
}
combined=data.frame(cbind(xs,ys))
colnames(combined)=c("x","y")
return(combined)
}}


 regression = smoother(new.X, new.Y,c=C,bw=Bandwidth)







 bootstrapmatrix1[(1 + (100 * (i - 1))):(100 * i),
 1:2] = cbind(regression$x, regression$y)

 }
 }
 if (Parametric == TRUE) {
 get.resid = function(x, y, x.close, y.hat) {
 return(y - y.hat[which.min(abs(x - x.close))])
 }



 regression = smoother(X[X < C], Y[X < C],c=C,bw=Bandwidth)


 residuals = rep(0, length(X[X < C]))
 for (i in 1:sum(X < C)) {
 residuals[i] = get.resid(X[X < C][i], Y[X < C][i],
 regression$x, regression$y)
 }
 for (i in 1:NBoots) {
 if (length(Cluster) == 0) {
 randompoints = sample(seq(1, sum(X < C)), sum(X <
 C), replace = TRUE)
 new.X = X[X < C]
 new.Y = rep(0, sum(X < C))
 for (k in 1:sum(X < C)) {
 new.Y[k] = regression$y[which.min(abs(X[X <
 C][k] - regression$x))] + residuals[randompoints[k]]
 }
 }
 if (length(Cluster) != 0) {
 if (length(names(table(as.vector(Cluster)))[table(as.vector(Cluster)) ==
 max(table(as.vector(Cluster)))]) == length(unique(Cluster))) {
 new.X = X[X < C]
 randomClusters = sample(unique(Cluster[X <
 C]), length(unique(Cluster[X < C])), replace = TRUE)
 new.residuals = c()
 for (j in 1:length(randomClusters)) {
 new.residuals = c(new.residuals, residuals[Cluster[X <
 C] == randomClusters[j]])
 }
 new.residuals = sample(new.residuals, length(new.residuals))
 new.Y = rep(0, sum(X < C))
 for (k in 1:sum(X < C)) {
 new.Y[k] = regression$y[which.min(abs(X[X <
 C][k] - regression$x))] + new.residuals[k]
 }
 }
 }






 regression2 = smoother(new.X, new.Y,c=C,bw=Bandwidth)




 bootstrapmatrix1[(1 + (100 * (i - 1))):(100 * i),
 1:2] = cbind(regression2$x, regression2$y)
 }
 }
 lowerreg = cbind(bootstrapmatrix1[1:100, ])
 upperreg = cbind(bootstrapmatrix1[1:100, ])
 for (i in 1:100) {
 pointsforthatx = bootstrapmatrix1[bootstrapmatrix1[,
 1] == bootstrapmatrix1[i, 1], ]
 lowerreg[i, 1:2] = c(pointsforthatx[1, 1], quantile(pointsforthatx[,
 2], prob = (1-Confidence.Level)/2, na.rm = TRUE))
 upperreg[i, 1:2] = c(pointsforthatx[1, 1], quantile(pointsforthatx[,
 2], prob = 0.5+Confidence.Level/2, na.rm = TRUE))
 }
 lowerreg = lowerreg[which(!is.na(lowerreg[, 2])), ]
 upperreg = upperreg[which(!is.na(upperreg[, 2])), ]

lowerreg[,2][lowerreg[,2]<Absolute.Y.Min]=Absolute.Y.Min
upperreg[,2][upperreg[,2]>Absolute.Y.Max]=Absolute.Y.Max







 bootstrapmatrix2 = matrix(0, nrow = 100 * NBoots, ncol = 2)
 if (Parametric == FALSE) {
 for (i in 1:NBoots) {
 if (length(Cluster) == 0) {
 randompoints = sample(seq(1, sum(X > C)), sum(X >
 C), replace = TRUE)
 new.X = X[X > C][randompoints]
 new.Y = Y[X > C][randompoints]
 }
 if (length(Cluster) != 0) {
 randomClusters = sample(unique(Cluster[X > C]),
 length(unique(Cluster[X > C])), replace = TRUE)
 new.X = c()
 for (j in 1:length(randomClusters)) {
 new.X = c(new.X, X[X > C][Cluster[X > C] ==
 randomClusters[j]])
 }
 new.Y = c()
 for (j in 1:length(randomClusters)) {
 new.Y = c(new.Y, Y[X > C][Cluster[X > C] ==
 randomClusters[j]])
 }
 }




 regression = smoother(new.X, new.Y, c=C, bw=Bandwidth)



 bootstrapmatrix2[(1 + (100 * (i - 1))):(100 * i),
 1:2] = cbind(regression$x, regression$y)
 }
 }
 if (Parametric == TRUE) {
 get.resid = function(x, y, x.close, y.hat) {
 return(y - y.hat[which.min(abs(x - x.close))])
 }

 regression2 = smoother(X[X > C], Y[X > C],c=C,bw=Bandwidth)


 residuals = rep(0, length(X[X > C]))
 for (i in 1:sum(X > C)) {
 residuals[i] = get.resid(X[X > C][i], Y[X > C][i],
 regression$x, regression$y)
 }
 for (i in 1:NBoots) {
 if (length(Cluster) == 0) {
 randompoints = sample(seq(1, sum(X > C)), sum(X >
 C), replace = TRUE)
 new.X = X[X > C]
 new.Y = rep(0, sum(X > C))
 for (k in 1:sum(X > C)) {
 new.Y[k] = regression$y[which.min(abs(X[X >
 C][k] - regression$x))] + residuals[randompoints[k]]
 }
 }
 if (length(Cluster) != 0) {
 new.X = X[X < C]
 randomClusters = sample(unique(Cluster[X > C]),
 length(unique(Cluster[X > C])), replace = TRUE)
 new.residuals = c()
 for (j in 1:length(randomClusters)) {
 new.residuals = c(new.residuals, residuals[Cluster[X >
 C] == randomClusters[j]])
 }
 new.residuals = sample(new.residuals, length(new.residuals))
 new.Y = rep(0, sum(X > C))
 for (k in 1:sum(X > C)) {
 new.Y[k] = regression$y[which.min(abs(X[X >
 C][k] - regression$x))] + new.residuals[k]
 }
 }

 regression2 = smoother(new.X, new.Y,c=C,bw=Bandwidth)



 bootstrapmatrix2[(1 + (100 * (i - 1))):(100 * i),
 1:2] = cbind(regression2$x, regression2$y)
 }
 }
 lowerreg2 = cbind(bootstrapmatrix2[1:100, ])
 upperreg2 = cbind(bootstrapmatrix2[1:100, ])
 for (i in 1:100) {
 pointsforthatx = bootstrapmatrix2[bootstrapmatrix2[,
 1] == bootstrapmatrix2[i, 1], ]
 lowerreg2[i, 1:2] = c(pointsforthatx[1, 1], quantile(pointsforthatx[,
 2], prob = (1-Confidence.Level)/2, na.rm = TRUE))
 upperreg2[i, 1:2] = c(pointsforthatx[1, 1], quantile(pointsforthatx[,
 2], prob = 0.5+Confidence.Level/2, na.rm = TRUE))
 }
 lowerreg2 = lowerreg2[which(!is.na(lowerreg2[, 2])), ]
 upperreg2 = upperreg2[which(!is.na(upperreg2[, 2])), ]

lowerreg2[,2][lowerreg2[,2]<Absolute.Y.Min]=Absolute.Y.Min
upperreg2[,2][upperreg2[,2]>Absolute.Y.Max]=Absolute.Y.Max


y1=smoother(X[X < C], Y[X < C],c=C,bw=Bandwidth)[,2]
y1[y1<Absolute.Y.Min]=Absolute.Y.Min
y1[y1>Absolute.Y.Max]=Absolute.Y.Max

y2=smoother(X[X > C], Y[X > C],c=C,bw=Bandwidth)[,2]
y2[y2<Absolute.Y.Min]=Absolute.Y.Min
y2[y2>Absolute.Y.Max]=Absolute.Y.Max


 if (Plot.Raw.Data == TRUE & Jitter == FALSE) {
	
plot=ggplot()+geom_point(aes_string(x=X[X < C],y=Y[X < C]),colour=Raw.Data.Colors[1],alpha=Raw.Data.Point.Density,size=Raw.Data.Point.Size)+xlab(xlab)+ylab(ylab)+ggtitle(Main)+theme(legend.position="none",plot.title = element_text(size=Title.Size))+geom_point(aes_string(x=X[X > C],y=Y[X > C]),colour=Raw.Data.Colors[2],size=Raw.Data.Point.Size,alpha=Raw.Data.Point.Density)+

geom_ribbon(aes_string(x=upperreg[,1],ymin=lowerreg[,2],ymax=upperreg[,2]),colour="gray",alpha=0.2)+# geom_line(aes_string(x=lowerreg[,1],y=lowerreg[,2]),linetype="dashed")+geom_line(aes_string(x=upperreg[,1],y=upperreg[,2]),linetype="dashed")+

geom_ribbon(aes_string(x=upperreg2[,1],ymin=lowerreg2[,2],ymax=upperreg2[,2]),colour="gray",alpha=0.2)+# geom_line(aes_string(x=lowerreg2[,1],y=lowerreg2[,2]),linetype="dashed")+geom_line(aes_string(x=upperreg2[,1],y=upperreg2[,2]),linetype="dashed")+

geom_line(aes_string(x=smoother(X[X < C], Y[X < C],c=C,bw=Bandwidth)[,1],y=y1),colour=Line.Colors[1])+geom_line(aes_string(x=smoother(X[X > C], Y[X > C],c=C,bw=Bandwidth)[,1],y=y2),colour=Line.Colors[2])+geom_vline(xintercept=C)+geom_hline(yintercept=0) +

scale_x_continuous(breaks=Tick.Marks,labels=Labels) + coord_cartesian(xlim = xlim, ylim = ylim)

 }
 if (Plot.Raw.Data == TRUE & Jitter == TRUE) {
 JX = jitter(X[X < C])
 JX[JX > C] = 2 * C - JX[JX > C]
 X[X < C] = JX
 JX = jitter(X[X > C])
 JX[JX < C] = 2 * C - JX[JX < C]
 X[X > C] = JX

plot=ggplot()+geom_point(aes(X[X < C],jitter(Y[X < C])),colour=Raw.Data.Colors[1],alpha=Raw.Data.Point.Density,size=Raw.Data.Point.Size)+xlab(xlab)+ylab(ylab)+ggtitle(Main)+theme(legend.position="none",plot.title = element_text(size=Title.Size))+geom_point(aes(X[X > C],jitter(Y[X > C])),colour=Raw.Data.Colors[2],size=Raw.Data.Point.Size,alpha=Raw.Data.Point.Density)+

geom_ribbon(aes_string(x=upperreg[,1],ymin=lowerreg[,2],ymax=upperreg[,2]),colour="gray",alpha=0.2)+# geom_line(aes_string(x=lowerreg[,1],y=lowerreg[,2]),linetype="dashed")+geom_line(aes_string(x=upperreg[,1],y=upperreg[,2]),linetype="dashed")+

geom_ribbon(aes_string(x=upperreg2[,1],ymin=lowerreg2[,2],ymax=upperreg2[,2]),colour="gray",alpha=0.2)+# geom_line(aes_string(x=lowerreg2[,1],y=lowerreg2[,2]),linetype="dashed")+geom_line(aes_string(x=upperreg2[,1],y=upperreg2[,2]),linetype="dashed")+

geom_line(aes_string(x=smoother(X[X < C], Y[X < C],c=C,bw=Bandwidth)[,1],y=y1),colour=Line.Colors[1])+geom_line(aes_string(x=smoother(X[X > C], Y[X > C],c=C,bw=Bandwidth)[,1],y=y2),colour=Line.Colors[2])+geom_vline(xintercept=C)+geom_hline(yintercept=0) +

scale_x_continuous(breaks=Tick.Marks,labels=Labels) + coord_cartesian(xlim = xlim, ylim = ylim)



 }
 if (Plot.Means == TRUE) {



plot=ggplot()+geom_point(aes(negative.midpoints,dat[, 2][midpoints < C]),size=pointsize[midpoints < C],colour=Mean.Colors[1])+xlab(xlab)+ylab(ylab)+ggtitle(Main)+theme(legend.position="none",plot.title = element_text(size=Title.Size))+geom_point(aes(positive.midpoints,dat[, 2][midpoints > C]),size=pointsize[midpoints > C],colour=Mean.Colors[2])+


geom_ribbon(aes_string(x=upperreg[,1],ymin=lowerreg[,2],ymax=upperreg[,2]),colour="gray",alpha=0.2)+# geom_line(aes_string(x=lowerreg[,1],y=lowerreg[,2]),linetype="dashed")+geom_line(aes_string(x=upperreg[,1],y=upperreg[,2]),linetype="dashed")+

geom_ribbon(aes_string(x=upperreg2[,1],ymin=lowerreg2[,2],ymax=upperreg2[,2]),colour="gray",alpha=0.2)+# geom_line(aes_string(x=lowerreg2[,1],y=lowerreg2[,2]),linetype="dashed")+geom_line(aes_string(x=upperreg2[,1],y=upperreg2[,2]),linetype="dashed")+

geom_line(aes_string(x=smoother(X[X < C], Y[X < C],c=C,bw=Bandwidth)[,1],y=y1),colour=Line.Colors[1])+geom_line(aes_string(x=smoother(X[X > C], Y[X > C],c=C,bw=Bandwidth)[,1],y=y2),colour=Line.Colors[2])+geom_vline(xintercept=C)+geom_hline(yintercept=0) +

scale_x_continuous(breaks=Tick.Marks,labels=Labels) + coord_cartesian(xlim = xlim, ylim = ylim)

}


 estimate=smoother(X[X > C], Y[X > C],c=C,bw=Bandwidth)[2][[1]][smoother(X[X > C], Y[X > C],c=C,bw=Bandwidth)[1][[1]]==C]-smoother(X[X < C], Y[X < C],c=C,bw=Bandwidth)[2][[1]][smoother(X[X < C], Y[X < C],c=C,bw=Bandwidth)[1][[1]]==C]

# if (length(Window) == 2) {
# abline(v = Window[1], lty = 2)
# abline(v = Window[2], lty = 2)
# }
# if(length(Labels)>0){axis(1, at = Tick.Marks, labels = Labels)}
   # if(length(Labels)==0){axis(1, at = Tick.Marks, labels = Tick.Marks)}
 p = 2 * length(which(bootstrapmatrix2[which(bootstrapmatrix2[,
 1] == C), 2] - bootstrapmatrix1[which(bootstrapmatrix1[,
 1] == C), 2] < 0))/NBoots
 if(p>1){p=2-p} 
    if(Type=="One-Sided"){p=p/2}
 p.text = paste("p=", signif(p, digits = 2), collapse = "")
 p.text = gsub("\\s", "", p.text)
 if (Plot.p.value == TRUE) {
 if (Loc.p.value == "BR" & length(ylim) == 1) {
 location = c(xlim[2] - (xlim[2] - xlim[1])/20, min(dat[,
 2]) + (max(dat[, 2]) - min(dat[, 2]))/25)
 }
 if (Loc.p.value == "BR" & length(ylim) == 2) {
 location = c(xlim[2] - (xlim[2] - xlim[1])/20, y = ylim[1])
 }
 if (Loc.p.value == "TR" & length(ylim) == 1) {
 location = c(xlim[2] - (xlim[2] - xlim[1])/20, max(dat[,
 2]) - (max(dat[, 2]) - min(dat[, 2]))/25)
 }
 if (Loc.p.value == "TR" & length(ylim) == 2) {
 location = c(xlim[2] - (xlim[2] - xlim[1])/20, y = ylim[2])
 }
 if (Loc.p.value == "BL" & length(ylim) == 1) {
 location = c(xlim[1] + (xlim[2] - xlim[1])/20, min(dat[,
 2]) + (max(dat[, 2]) - min(dat[, 2]))/25)
 }
 if (Loc.p.value == "BL" & length(ylim) == 2) {
 location = c(xlim[1] + (xlim[2] - xlim[1])/20, y = ylim[1])
 }
 if (Loc.p.value == "TL" & length(ylim) == 1) {
 location = c(xlim[1] + (xlim[2] - xlim[1])/20, max(dat[,
 2]) - (max(dat[, 2]) - min(dat[, 2]))/25)
 }
 if (Loc.p.value == "TL" & length(ylim) == 2) {
 location = c(xlim[1] + (xlim[2] - xlim[1])/20, y = ylim[2])
 }
 #text(location[1], location[2], p.text, cex = 0.75)
 }
 output=list(estimate,signif(p, digits = 2))
 names(output)=c("Estimate","p.value")

if(NBoots<10000){print(paste("To increase the precision of the p-value, increase the number of bootstrapped samples. Currently, NBoots =",NBoots))}

if(sum(X[X<(C+Bandwidth)&X>C])<100|sum(X[X<C&X>(X-Bandwidth)])<100){if(Smoother=="Local Linear"|Smoother=="Local Polynomial"){print("Note: Local linear regression and local polynomial regression can have high variablility in cases where there is not a large sample size around the cut-piont. Kernel Reggression may be more appropiate.")}}

if(sum(X==C)>0){print("Note: Some observations are on the cut-point. These units have been dropped from the test.")}


print(output)

plot

}
