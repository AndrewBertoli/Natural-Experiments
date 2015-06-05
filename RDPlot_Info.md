### RDPlot()

RDPlot() creates a regression discontinuity graph with regression lines on either side of the cut-point. The user can choose the smoother, kernel, and bandwidth, along with whether to plot the raw data or the binned-means. Confidence intevals and p-values are estimated using parametic or non-parametric bootstrapping.

##### Usage

RDPlot=function (X, Y, C = 0, xlim = range(X), ylim = "Automatic", 
xlab = "Forcing Variable", ylab = "Outcome", Main = "Regression Discontinuity Plot", 
Plot.Means = TRUE, Plot.Raw.Data = FALSE, Mean.Colors = c("blue", "red"), 
Raw.Data.Colors = c("lightblue", "pink"), Point.Size = 0.25, 
Shade.Color = "gray87", Window = "None", Plot.p.value = TRUE, 
Tick.Marks = c(-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
Bandwidth = 2, NBoots = "Automatic", Breaks = "Automatic", Parametric = FALSE,
Smoother="Kernel", Kernel=NULL, Poly.Order=NULL, Cluster = NULL, Jitter = FALSE, 
Loc.p.value = "BR",...) 

##### Arguments

X:	A vector containing the values of the forcing variable for each unit.

Y:	A vector containing the outcomes for each unit.

C:	The cut-point.

xlim: A vector of length 2 giving the range of X to be graphed.

ylim: A vector of length 2 giving the range of Y to be graphed.

xlab: The label for the x-axis.

ylab: The label for the y-axis.

Main: The title of the graph.

Plot.Means: If TRUE, the means of groups of units with similar values of the forcing variable will be plotted. The groups will be determined based on Breaks.

Mean.Colors:  A vector of length 2 specifying the colors of the untreated and treated unit means.

Plot.Raw.Data:  If TRUE, the raw data will be plotted.

Raw.Data.Colors:  A vector of length 2 specifying the colors of the untreated and treated units.

Point.Size:	The size of the points representing the means. These points will be sized in proportion to how many units they represent.

Shade.Color:	The color of the confidence interval.

Window: A vector of length 2 specifying the x-values for the regression discontinuity window.

Plot.p.value:	If TRUE, the p-value will be plotted on the graph.

Tick.Marks:	The locatoins to put the tick marks on the x-axis.

Bandwidth:  The bandwidth to be used in constructing the regression lines.

NBoots:	The number of bootstrapped samples to draw.

Breaks:	The breaks for grouping the units when plotting the group means

Parametric:	If TRUE, the method used will be parametric bootstrapping.

Smoother:	The regression smoother. Options include "Nearest Neighbor", "Kernel", "Local Linear", and "Local Polynomial".

Kernel: The regression kernel. Options include "Box", "Tricube", "Triweight", "Triangular", "Epanechnikov", and "Quartic".

Poly.Order: The order of the polynomical if Smoother="Local Polynomial". This value can be any integer form 1 to 9.

Cluster: A vector with the same length as X that will specify a category for each unit that will be used in the bootstrapping stage.

Jitter:	If TRUE, the raw data points will be jittered.

Loc.p.value:	The location where the p-value will be plotted on the graph ("BR"-Bottom Right, "TR"-Top Right, "BL"-Bottom Left, "TL"-Top Left)

##### Value

Estimate:	The estimated treatment effect at the cut-point

p.value:	The two-tailed p-value for the estimated treatement effect



##### Example

data(WorldCup)

qualifiers=data[data$Treat==1,]
nonqualifiers=data[data$Treat==0,]

nonqualifiers=nonqualifiers[qualifiers$Score>=5,]
qualifiers=qualifiers[qualifiers$Score>=5,]
 
nonqualifiers[nonqualifiers$PointsFromCutpoint==0,]$PointsFromCutpoint=-0.0001
qualifiers[qualifiers$PointsFromCutpoint==0,]$PointsFromCutpoint=0.0001

sample=rbind(qualifiers,nonqualifiers)
	
op=par(mfrow=c(1,3), mar=c(4,3.9,3,0),oma=c(0,1,2,1))

RDPlot(X=sample$PointsFromCutpoint,Y=sample$AGG3YearsBefore,xlim=c(-3.7,3.7),
ylim=c(-1.5,1.5), Main="3 Years Before Qualification", xlab="Points Above/Below Cut-Point",
ylab="Militarized Interstate Disputes Initiated",Bandwidth=2.85621, Window=c(-2.5,2.5))

RDPlot(X=sample$PointsFromCutpoint,Y=sample$AGG3YearsAfter,xlim=c(-3.7,3.7),ylim=c(-1.5,1.5), 
Main="3 Years After Qualification",xlab="Points Above/Below Cut-Point", ylab="", Bandwidth=3.2316, 
Window=c(-2.5,2.5))

RDPlot(X=sample$PointsFromCutpoint,Y=sample$AGG3YearsAfter-sample$AGG3YearsBefore,
xlim=c(-3.7,3.7),ylim=c(-1.5,1.5), Main="Change in Aggression", xlab="Points Above/Below Cut-Point",
ylab="", Bandwidth=4.1937, Window=c(-2.5,2.5))

title("Figure 5. Change in Aggression for the World Cup",outer=TRUE,cex.main=2)	
