
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
