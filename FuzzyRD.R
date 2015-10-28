# This functions creates two graphs that illustrate a Fuzzy RD. The first graph
# shows the relationship between the forcing variable and treatment, and the second
# graph shows the relationship between the forcing variable and the outcome. The function
# also calculates the estimated LATE for compliers and the p-value. It requires RDPlot
# to be loaded into the workspace.

FuzzyRD=function (X, T, Y, C = 0, xlim = range(X), ylim = list("Automatic",c(0,1)), 
xlab = "Forcing Variable", Main.Title="Fuzzy RD Plot", ITT.ylab = "Outcome", 
ITT.Title = "Discontinuity in Outcome", Treatment.ylab = "Treated", 
Treatment.Title = "Discontinuity in Treatment", Plot.Means = c(TRUE,TRUE), 
Plot.Raw.Data = c(FALSE,FALSE), Mean.Colors = c("blue", "red", "blue", "red"), 
Raw.Data.Colors = c("lightblue", "pink", "lightblue", "pink"), Point.Size = c(0.25,0.25), 
Shade.Color = c("gray87","gray87"), Window = "None", 
Tick.Marks = c(-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
Bandwidth = c(2,2), NBoots = 10000, Breaks = "Automatic",  Parametric = FALSE, 
Cluster = NULL, Jitter = c(FALSE,FALSE)){

mfrow=par()$mfrow
mar=par()$mar
oma=par()$oma
mgp=par()$mgp

op=par(mfrow=c(1,2), mar=c(3.9,3.9,2,2), oma=c(0,3,3,1),mgp=c(2.3,1,0))

Treatment.Status=RDPlot(X=X, Y=T, C=C, xlim=range(X), ylim=ylim[[2]], xlab=xlab, ylab=Treatment.ylab, 
Main=Treatment.Title, Plot.Means=Plot.Means[2], Plot.Raw.Data=Plot.Raw.Data[2], Mean.Colors=Mean.Colors[3:4], 
Raw.Data.Colors=Raw.Data.Colors[3:4], Point.Size=Point.Size[2], Shade.Color=Shade.Color[2], Window=Window, 
Plot.p.value=FALSE, Tick.Marks=Tick.Marks, Bandwidth=Bandwidth[2], NBoots=NBoots, Breaks=Breaks,
Parametric=Parametric, Cluster=Cluster, Jitter=Jitter[2])

ITT=RDPlot(X=X, Y=Y, C=C, xlim=range(X), ylim=ylim[[1]], xlab=xlab, ylab=ITT.ylab, Main=ITT.Title, 
Plot.Means=Plot.Means[1], Plot.Raw.Data=Plot.Raw.Data[1], Mean.Colors=Mean.Colors[1:2], 
Raw.Data.Colors=Raw.Data.Colors[1:2], Point.Size=Point.Size[1], Shade.Color=Shade.Color[1],
Window=Window, Plot.p.value=FALSE, Tick.Marks=Tick.Marks, Bandwidth=Bandwidth[1], NBoots=NBoots, 
Breaks=Breaks,Parametric=Parametric, Cluster=Cluster, Jitter=Jitter[1])

title(Main.Title,outer=TRUE,cex.main=2)

estimate=ITT$Estimate/Treatment.Status$Estimate
p=ITT$p.value

op=par(mfrow=mfrow, mar=mar, oma=oma, mgp=mgp)

output=list(estimate,signif(p, digits = 2))
names(output)=c("Estimate","p.value")
return(output)
}
