# CovidStrategyCalculator
A standalone GUI to assess Covid-19 testing and quarantine strategies for contact person management.

## Table of contents
* [Usage](#1-usage)
* [Results](#2-results)
* [Building from source](#3-building-from-source)

## 1. Usage
The user can set the variables of the simulation in the left-upper panel.

In the tab "Input" the strategy is selected. The user has to provide three variables:
* The time-delay since infection. The time of infection is modeled as the time of exposure with an assumed one-probability of disease transmission.
* The duration of the quarantine in days.  
* The days on which a PCR test is to be performed. To ensure meaningful results, the time it takes to get test results back should be taken into consideration here. For example, with a delay of one day between testing and receiving results the (last) test should be scheduled at the second to last day of the quarantine.


In the tab "Parameters" the internal model parameters can be tweaked, these include:
* The duration of the incubation period, which is defined as the period before symptom onset.  
* The duration of the pre-detectable period, which is defined as the period during incubation in which infected individuals cannot yet be identified by a PCR-test. The duration of the pre-detable period is specified as a percentage of the total incubation period.
* The duration of the symptomatic period.
* The percentage of infections in which no symptoms are developed (asymptomatic disease). Asymptomatic infections are assumed to follow the same disease time course as symptomatic ones.
* The duration of the post-symptomatic period, which is defined as the period after symptoms have passed but during which the disease is still detectable (positive PCR-test).
* The (average) sensitivity of the PCR-test, which is defined as the percentage of cases in which a test returns positive, given the subject is infected.
* The (average) specificity of the PCR-test, which is defined as the percentage of cases in which a test returns negative, given the subject is not infected.

The user can set mean, lower-extreme and upper-extreme values for the duration of the different disease phases.
The application runs three separate simulations based on these parameters: an expected scenario based on the provided mean value, a best-case based on the lower extreme value and a worst-case scenario based on the upper extreme value.

## 2. Results
The simulation calculates the risk reduction for the provided strategy. The risk reduction is reported both as percental reduction and as factor:

<img src="https://render.githubusercontent.com/render/math?math=%5Ctext%7Brisk%20reduction%20%5B%5C%25%5D%7D%20%3D%20%5Cfrac%7B%5Ctext%7Binitial%20risk%7D%20-%20%5Ctext%7Bend%20risk%7D%7D%7B%5Ctext%7Binitial%20risk%7D%7D%20%5Ccdot%20100%5C%25">

<img src="https://render.githubusercontent.com/render/math?math=%5Ctext%7Brisk%20reduction%20%5Bfactor%5D%7D%20%3D%20%5Cfrac%7B%5Ctext%7Binitial%20risk%7D%7D%7B%5Ctext%7Bend%20risk%7D%7D">

We define risk as the probability that an infected individual without symptoms is released into society. The initial risk is the risk at day zero of the simulation (probability of 1) and the end risk is the risk at the last day of the quarantine.


The results are reported in the upper-right panel. Per risk reduction three values are stated, the first value corresponds to the expected scenario based on the provided mean duration times, the second and third reflect the best- and worst-case scenario respectively based on the provided extreme values.

In the lower panel two time trajectories are shown. The red curve represents the population that is- or will become infectious, i.e. individuals in the incubation or symptomatic phase. The black curve depicts the fraction of this population that can be detected through testing.
Below the graph the size of the detectable population is stated as a percentage of the total population (combined population of incubation, symptomatic, post-symptomatic and post-detection). Again, the first value reflects the mean scenario, the second and third value reflect the best- and worst-case scenario respectively based on the provided extreme values.  


## 3. Building from source
The CovidStrategyCalculator application can be compiled from source using the Qt5 framework. CovidStrategyCalculator requires the Eigen 3.3.7 library. When building from source, the path to the Eigen source code (root folder) has to be set in the `src/CovidTestCalculator.pro` file:

```{c}
INCLUDEPATH += path/to/eigen
```
### Versions
This application was developed using:
* Eigen 3.3.7 (`https://gitlab.com/libeigen/eigen/-/releases#3.3.7`)
* Qt 5.9.5 (`https://download.qt.io/archive/qt/5.12/5.12.9/`)

Other versions of these two libraries might well work, but have not been tested.
