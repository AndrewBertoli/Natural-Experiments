
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

Other.Tests:	A vector with the names of other tests to be included, such as PermutationTest. The functions for these tests must be loaded into R's gloval enviorment. The functions should take a treatment and control group as the first two arguments and return a p-value.

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
