# Natural-Experiments
 This file contains a number of R functions for analyzing natural experiments and regression discontinuities. The functions are described below. If you have any questions, please email me at abertoli@berkeley.edu.

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


### Sens()

Sens() computes the p-values for different Gamma levels for paired or unpaired data. The argument "Gamma" is the Gamma level, defined as the maximum difference in treatment odds between any two units in the sample. So when Gamma=j (or 1/j), the maximum p-value is computed under the assumption that no unit was more than j-times more likely to be treated than any other unit in the sample. When "Type" is set at "mean", the function will use the mean as the test statistic. When "Type" is set at "rank", the function will do the Wilcoxon Signed-Rank Test or Rank Sum Tests, which are less sensitive to outliers.

##### Usage

Sens(Treatment.Outcomes, Control.Outcomes, Gamma = 1, Type = "Mean", Paired = TRUE, 
na.rm = FALSE, Simulations = 1000000)

##### Arguments

Treatment.Outcomes:	A vector containing the outcomes for the treated units.

Control.Outcomes: A vector containing the outcomes for the control units. If Paired=TRUE, then this vector must have the same length as Treatment.Outcomes.

Gamma:	The maximum difference in treatment odds between any two units in the sample.

Type:	The test statistic that should be used, either "Mean" or "Rank".

Paired:	If TRUE, the test will treat the units as paired.

na.rm:	If TRUE, NAs will be removed. Note that removing NAs may induce bias.

Simulations:	The number of Monte Carlo simulations to perform.

##### References

Rosenbaum, Paul R. Observational studies. Springer New York, 2002.

##### Example

data(WorldCup)

Sens(t$AGGAfter-t$AGGBefore,c$AGGAfter-c$AGGBefore, Gamma=1.5, 
Paired=TRUE, Type="Mean")

### BalancePlot

This function shows the balance between the treatment and control groups by plotting the distrubtion of p-values for all of the covariates. By presenting the data in this way, it is easier to assess whether the p-values appear to be distrubed uniformly between 0 and 1, as would be expected in an experiment.

##### Usage

BalancePlot=function(Data, Treat, Covariates, Names.To.Print, Shade.Color = "blanchedalmond",
Title = "Balance Plot", Title.Size = 1.2, na.rm = FALSE,
Built.In.Tests = c("T.Test"), Point.Color = "red", Other.Tests = NULL,
pch = NULL, Year.Covariates = NULL, Observational.Data = NULL,
Observational.Treat = NULL, Observational.Point.Color = "blue",
Sample.Name = "Sub-Sample", O.Name = "All Units", Legend = FALSE,
Paired = FALSE, Observational.Paired = FALSE)

##### Arguments

Data:	A dataframe where the units are represented by rows and the covariates are included in the columns.

Treat:	Either the name of the (0,1) treatment variable in the dataframe or a vector of 0's and 1's indicating the treatment assignments for the units.

Covariates:	A vector listing the covariates, as they are named in the dataframe.

Names.To.Print:	A vector listing the names to be printed for the covariates.

Shade.Color	
The color of the alternating background bars.

Title:	The title of the balance plot.

Title.Size:	The size of the title.

na.rm:	Logical argument indicating whether NAs in the covariates should be ignored. Best practice is usually to impute and match on missingness.

Built.In.Tests:	Can be set at T.Test and/or KS.Test, or neither.

Point.Color:	The color of the points that will be plotted.

Other.Tests:	A vector with the names of other tests to be included, such as Permutation.Test. The functions for these tests must be loaded into R's gloval enviorment. The functions should take a treatment and control group as the first two arguments and return a p-value.

pch:	A vector specifying the style of the points for the additional tests.

Year.Covariates:	The names of any covariates that are years. These covariates will be printed in a different format.

Observational.Data:	A separate dataframe for the non-experimental data.

Observational.Treat:	Either the name of the (0,1) treatment variable in Observational.Data or a vector of 0's and 1's indicating the treatment assignments for the observational units.

Observational.Point.Color:	The color of the points that will be plotted for the observational data.

Legend:	Logical indicating whether to include a legend.

Sample.Name:	If specified, the name of the sample data to be printed in the legend.

O.Name:	If specified, the name of the observational data to be printed in the legend.

Paired:	If TRUE, then the t-tests and KS-tests will be paired.

Observational.Paired:	If TRUE, then the t-tests and KS-tests will be paired for the observational data.

##### Example

data(WorldCup)
	
BalancePlot(Data=sample, Treat=sample$Treat, Title="Balance Between the Qualifiers and Non-qualifiers",
Covariates=c("Irst",'Milex','Milper','Tpop','Upop','BirthRate','DeathRate','InfantMortality','Energy',
'Imports','Exports','LandArea','CINC','Democracy','GreatPower','EngagedCivilWar','EndedCivilWar',
'EntranceYear','SexRatio','LifeExpectancy','MedianAge','Alliances','USAlly','SoccerMostPopular',
'PrevAppear','AGGYearBefore','AGG3YearsBefore','AGG5YearsBefore'), Names.To.Print=
c('Iron and Steel Production','Military Expenditures', 'Military Personnel', 'Total Population', 
'Urban Population', 'Birth Rate', 'Death Rate', 'Infant Mortality', 'Energy Production', 'Imports', 
'Exports', 'Land Area', 'Material Power Score', 'Level of Democracy', 'Great Power Status', 
'Engaged in Civil War', 'Resolved Civil War', 'Year of State Formation', 'Sex Ratio', 'Life Expectancy', 
'Median Age', "Number of Alliances", "U.S. Ally", "Soccer Most Popular Sport", 'Appearance at Previous World Cup',  
'MIDs Initiated in the Year Before', 'MIDs Initiated in the 3 Years Before', 'MIDs Initiated in the 5 Years Before'),
Shade.Color="cadetblue2", na.rm=FALSE, Built.In.Tests=c("T.Test"), Point.Color="black", 
Sample.Name="RD Sample", Year.Covariates=c("EntranceYear"), Observational.Data=data, 
Observational.Treat=data$Treat, Observational.Point.Color="gray90", Legend=TRUE, Paired=TRUE, 
Observational.Paired=TRUE)

### Permutation.Test()

Permutation.Test computes p-values using Monte Carlo simulation under the sharp null hypothesis of no treatment effect.

##### Usage

Permutation.Test(Treatment.Outcomes, Control.Outcomes, Paired = FALSE, 
Simulations = 100000, na.rm = FALSE, Output = "p")

##### Arguments

Treatment.Outcomes:	A vector of treatment outcomes.

Control.Outcomes:	A vector of control outcomes.

Paired:	If true, the test will be a paired permuation test.

Simulations:	The number of simulations to carry out.

na.rm:	If TRUE, NAs will be removed. Dropping NAs may induce bias.

Output:	If Output="p", just the p-value will be returned. If Output="Full", both the estimate and the p-value will be returned.



##### Example

data(WorldCup)
	
Permutation.Test(Treatment.Outcomes=t$AGGAfter-t$AGGBefore, Control.Outcomes=c$AGGAfter-c$AGGBefore, 
Paired = TRUE, Simulations = 100000, na.rm = FALSE, Output = "Full")	
