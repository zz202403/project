# Project 
HW7
A GitHub repository for collaborative health data analysis.

## Objective
  
The focus will be on using public data to explore the significance of "self-reported health" as a health indicator. Applying the skills we've developed over the past five weeks to create a GitHub repository project.

## Data Source

### Survey Data:

- Import the survey data from the 1999-2000 National Health and Nutrition Examination Survey (NHANES):
```stata
import sasxport5 "https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/DEMO.XPT", clear
```

### Mortality Follow-up Data:
Obtain follow-up mortality data to analyze over a 20-year period from the [National Center for Health Statistics (NCHS)](https://ftp.cdc.gov/pub/).
- Health Statistics
  - NCHS
    - Datalinkage
      - Linked Mortality
        - ```NHANES_1999_2000_MORT_2019_PUBLIC.dat```
		- ```Stata_ReadInProgramAllSurveys.do```

```stata
 //data
 global mort_1999_2000 https://ftp.cdc.gov/pub/HEALTH_STATISTICS/NCHS/datalinkage/linked_mortality/NHANES_1999_2000_MORT_2019_PUBLIC.dat
 //code
 cat https://ftp.cdc.gov/pub/HEALTH_STATISTICS/NCHS/datalinkage/linked_mortality/Stata_ReadInProgramAllSurveys.do
```
## Code Development
### Edit and Rename Provided Script:
Download, modify, and upload the provided [Stata .do file](https://ftp.cdc.gov/pub/HEALTH_STATISTICS/NCHS/datalinkage/linked_mortality/Stata_ReadInProgramAllSurveys.do) for linking the DEMO.XPT data to mortality follow-up data. 
- Rename this file to ```followup.do```
- modify ```followup.do```
  - modify working directory
  - modify survey data 
- upload modfied ```followup.do```

### Data Merging:
Execute the following Stata code to merge the survey data with the mortality data, ensuring alignment on the unique sequence numbers:
```stata
//use your own username/project repo instead of the class repo below
global repo "https://github.com/jhustata/intermediate/raw/main/"
do ${repo}followup.do
save followup, replace 
import sasxport5 "https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/DEMO.XPT", clear
merge 1:1 seqn using followup
lookfor follow
```

## Following Analysis  

+ Inference
  - Employ 95% CI and p-values basis
  ```stata
  merge 1:1 seqn using demo_mortality, nogen
  sts graph, by(huq010) fail
  stcox i.huq010
  ```
  - Brief conclusion basis
  ```stata
  import sasxport5 "https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/HUQ.XPT", clear 
  huq010 
  desc huq010
  codebook huq010
  ```

## Survival analysis: non-parametric, semi-parametric, parametric
+ click [here](dyndoc1.html) to view nonparametric and semiparametric risk estimates from Stata
