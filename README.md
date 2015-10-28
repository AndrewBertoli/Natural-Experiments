# Natural-Experiments
 This file contains a number of R functions for analyzing natural experiments and regression discontinuities. The functions are described below. If you have any questions, please email me at abertoli@berkeley.edu.

### RDPlot()

RDPlot() creates a regression discontinuity graph with regression lines on either side of the cut-point. The user can choose the smoother, kernel, and bandwidth, along with whether to plot the raw data or the binned-means. Confidence intevals and p-values are estimated using parametic or non-parametric bootstrapping.


### ![alt tag](https://cloud.githubusercontent.com/assets/7791421/7993424/e9b03a48-0aba-11e5-8ea3-78962dbf99e3.jpg)

### FuzzyRD()

FuzzyRD() creates two regression discontinuity graphs that illustrate the treatment effect for a fuzzy RD design. The first graph shows the relationship between the forcing variable and the treatment, and the second graph shows the relationship between the forcing variable and the outcome. The function also computes an estimated treatment effect and p-value for the fuzzy RD.

### Sens()

Sens() computes the p-values for different Gamma levels for paired or unpaired data. The argument "Gamma" is the Gamma level, defined as the maximum difference in treatment odds between any two units in the sample. So when Gamma=j (or 1/j), the maximum p-value is computed under the assumption that no unit was more than j-times more likely to be treated than any other unit in the sample. When "Type" is set at "mean", the function will use the mean as the test statistic. When "Type" is set at "rank", the function will do the Wilcoxon Signed-Rank Test or Rank Sum Tests, which are less sensitive to outliers.

### BalancePlot()

This function shows the balance between the treatment and control groups by plotting the distrubtion of p-values for all of the covariates. By presenting the data in this way, it is easier to assess whether the p-values appear to be distrubed uniformly between 0 and 1, as would be expected in an experiment.

### Permutation.Test()

Permutation.Test computes p-values using Monte Carlo simulation under the sharp null hypothesis of no treatment effect.
