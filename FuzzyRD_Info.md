### FuzzyRD

FuzzyRD() creates two graphs to illustrate the treatment effect at the cut-point. The first plots the outcome 
as a function of the forcing variable, which provides an intention to treat estimate.  The second graph plots 
treatment status against the outcome to illustrate the extent to which treatment assignemnt changed at the 
cut-point. FuzzyRD also returns an estimated treatment effect and a p-value.

##### Usage

FuzzyRD(X, T, Y, C = 0, xlim = range(X), ylim = list("Automatic", c(0, 1)), 
xlab = "Forcing Variable", Main.Title = "Fuzzy RD Plot", ITT.ylab = "Outcome", 
ITT.Title = "Discontinuity in Outcome", Treatment.ylab = "Treated", 
Treatment.Title = "Discontinuity in Treatment", Plot.Means = c(TRUE, TRUE),  
Plot.Raw.Data = c(FALSE, FALSE), Mean.Colors = c("blue", "red", "blue", "red"),  
Raw.Data.Colors = c("lightblue", "pink", "lightblue", "pink"), 
Point.Size = c(0.25, 0.25), Shade.Color = c("gray87", "gray87"), Window = "None", 
Tick.Marks = c(-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
Bandwidth = c(2, 2), NBoots = 10000, Breaks = "Automatic", Parametric = FALSE, 
Cluster = NULL, Jitter = c(FALSE, FALSE))

##### Arguments

X:  A vector containing the values of the forcing variable for each unit.

T:  A vector containing the treatment status for each unit.

Y:  A vector containing the outcomes for each unit.

C:  The cut-point.

xlim: A vector of length 2 giving the range of X to be graphed. 

ylim: A vector of length 2 giving the range of Y to be graphed. 

xlab: The label for the x-axis.

Main: The title of the graph.

ITT.ylab: The label for the y-axis on the intention to treat graph. 

ITT.Title:  The title for the the intention to treat graph. 

Treatment.ylab: The label for the y-axis on the the treatment status graph. 

Treatment.Title:  The title for the the treatment status graph. 

Plot.Means: A logical vector of length two specifying whether to plot the means of groups of units with similar values of the forcing variable. The groups will be determined based on Breaks.

Plot.Raw.Data:  A logical vector of length two specifying whether the raw data should be plotted in the two graphs

Mean.Colors:  A vector of length 4 specifying the colors of the untreated and treated unit means in the two graphs.

Raw.Data.Colors:  A vector of length 4 specifying the colors of the untreated and treated units in the two graphs.

Point.Size: A vector of length two specifying the size of the points representing the means. These points will be sized in proportion to how many units they represent. 

Shade.Color:  A vector of length two specifying the colors of the confidence intervals for the two graphs.

Window: A vector of length 2 specifying the x-values for the regression discontinuity window.

Tick.Marks: The locatoins to put the tick marks on the x-axis.

Bandwidth:  A vector of length 2 specifying the bandwidth to be used when constructing the regression lines in the two graphs.

NBoots: The number of bootstrapped samples to draw.

Breaks: The breaks for grouping the units when plotting the group means

Parametric  If TRUE, the method used will be parametric bootstrapping.

Cluster:  A vector with the same length as X that will specify a category for each unit that will be used in the bootstrapping stage.

Jitter: If logical vector of length 2. If either value is TRUE, the raw data points for that graph will be jittered.


##### Value

Estimate: The estimated effect of the treatment on the outcome.

p.value:  The p-value for this estimate, which is the same p-value as for the ITT estimate.
