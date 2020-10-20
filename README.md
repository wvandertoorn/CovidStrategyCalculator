# CovidStrategyCalculator
A standalone GUI to assess SARS-CoV2 testing- and  quarantine  strategies for incoming travelers, contact person management and de-isolation.

The CovidStrategyCalculator (CSC) calculates the residual risk (probability of being infectious upon release from quarantine or isolation), risk reduction and test efficacy for arbitrary testing and  quarantine strategies. The user determines whether the time of exposure/infection (or symptom onset) is known, the intended duration of isolation/quarantine and the time points for SARS-CoV2 PCR tests. Underneath, we implemented the analytical solution of a stochastic transit compartment model of the infection time course that captures published temporal changes and variability in test sensitivities, incubation- and infectious periods, as well as times to symptom onset. The infection time course is modeled using five compartments: a predetection phase (not infectious, negative PCR), a presymptomatic phase (infectious, positive PCR), a symptomatic phase (infectious, positive PCR), a postsymptomatic phase (not infectious, positive PCR) and a postdetection phase (not infectious, negative PCR).

**Disclaimer: CSC is in active development and is pending peer review. We strongly advice to re-build the tool from source before usage to ensure you have the latest version, and to use caution and sanity check your results.**

## Table of contents
* [Usage](#1-usage)
* [Results](#2-results)
* [Building from source](#3-building-from-source)

## 1. Usage
The model configuration can be viewed and changed in the upper left corner of the main window.
The configuration is split over three tabs: Strategy, Parameters and Prevalence estimator.

CSC's main task is to calculate the residual risk of a strategy.
The residual risk is defined as the probability of releasing infectious individuals, or individuals who are yet to become infectious, from quarantine or isolation.

### Mode
We provide three scenarios, or 'modes', in which the residual risk of a strategy can be assessed. The mode of simulation can be selected in the drop down menu in the first row of the Strategy tab. Each mode has a different use-case:

* 'exposure': This mode is meant to assess strategies for contact person management. The simulation starts at the date of exposure and assumes successful transmission, i.e. there is a 100% probability of becoming infectious at the start of the simulation. Exposure mode assumes symptomatic screening. This means that individuals who develop symptoms are assumed to go into isolation. They will not be released at the end of the quarantine and thus do not pose risk.
* 'symptom onset': This mode is meant to assess de-isolation strategies. The simulation starts at the date of symptom onset, i.e. there is a 100% probability of being infectious (and symptomatic) at the start of the simulation.
* 'prevalence estimation': This mode is meant to assess strategies for incoming travelers. This simulation is two-fold. A preliminary simulation is run based on the incidence reports of the last five weeks. Based on these reports, the prevalence of individuals that are- or will become infectious, and the initial states of the main simulation are determined. Prevalence estimation mode assumes symptomatic screening. Incoming travelers with symptoms are assumed to be denied entry or to go into isolation. Travelers who develop symptoms during the duration of the quarantine are assumed to go into isolation as well. In both cases, these individuals do no longer pose risk.

### Strategy
The intended duration of isolation/quarantine and the time points for SARS-CoV2 PCR tests can be set in the Strategy tab. By selecting a checkbox at day x, you indicate that a test should be performed on day x. You can plan multiple tests.  
Additionally, a time delay since the start of the simulation can be enforced. For example, by setting the 'time passed since exposure' to 3 days, the simulation starts at day -3.

### Parameters
In the Parameters tab, the user can set the parameters of model, such as the residence times per phase of the infection and testing statistics such as the average sensitivity and specificity.
These parameters can be tweaked to represent your data and are the responsibly of the user. We have merely provided realistic default values.

For the residence time parameters, you can provide a mean estimate and extreme values. The extreme values are used to calculate the uncertainty range of results.   

The parameter for the percentage of asymptomatic infections (last row) is only used in the modes 'exposure' and 'prevalence estimation'. This parameter determines which percentage of individuals in the symptomatic phase is NOT filtered out by the symptomatic screening. By changing this parameter to 100%, the assumption of symptomatic screening is dropped.

### Prevalence estimator
This tab is part of the mode 'prevalence estimation' but can also be used on its own to estimate the prevalence of individuals that are- or will become infectious, and the percent of the population in each phase of infection for a given region. The estimation is based on the incidence reports of the last 5 weeks. The reported cases are transformed to probabilities per day and distributed over the first four compartment of the model. For each of the last 35 days, a simulation is run to propagate the states said number of days. The final states of these simulations are summed and reported.

## 2. Results
The simulation calculates the residual risk and the fold risk reduction for the provided strategy. The results are reported in the result log in the upper right corner. The pre-procedure risk is the risk at the start of the simulation and the post-procedure risk is the risk at the last day of the quarantine/isolation.
Per result three values are reported: a mean estimate and the range of uncertainty based on the provided extreme values.

In the lower panel a chart and table are shown.
The chart can be used to determine the best time to place a test. The red curve represents the probability of becoming- or already being infectious, i.e. the probability of being in the incubation or symptomatic phase of the disease. No symptomatic screening is assumed. The black curve depicts the probability of infection and that infection being detectable.
Below the graph the test efficacy per day is stated. Again, three values are reported: a mean estimate and the range of uncertainty based on the provided extreme values.


## 3. Building from source
The CSC application can be compiled from source using the Qt5 framework. CSC requires the Eigen 3.3.7 library. When building from source, the path to the Eigen source code (root folder) has to be set in the
`src/CovidTestCalculator.pro` file:

```{c}
INCLUDEPATH += path/to/eigen
```
### Versions
This application was developed using:
* Eigen 3.3.7 (`https://gitlab.com/libeigen/eigen/-/releases#3.3.7`)
* Qt 5.9.5 (`https://download.qt.io/archive/qt/5.12/5.12.9/`)

Other versions of these two libraries might work, but have not been tested.
