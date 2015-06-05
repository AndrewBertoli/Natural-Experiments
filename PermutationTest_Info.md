
### PermutationTest()

PermutationTest() computes p-values using Monte Carlo simulation under the sharp null hypothesis of no treatment effect.

##### Usage

PermutationTest(Treatment.Outcomes, Control.Outcomes, Paired = FALSE, 
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
	
PermutationTest(Treatment.Outcomes=t$AGGAfter-t$AGGBefore, Control.Outcomes=c$AGGAfter-c$AGGBefore, 
Paired = TRUE, Simulations = 100000, na.rm = FALSE, Output = "Full")	
