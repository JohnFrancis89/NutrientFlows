---
title: "Assessing the relationship of Nitrogen and Posphorus in freshwater rivers under varied flow conditions"
output:
  html_document:
    df_print: paged
    toc: true
    toc_float: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
```

# Introduction 

In terrestrial and aquatic ecosystems, Nitrogen (N) and Pohsphorus (P) are usually the limiting factors on biotic growth (Marklein et al. 2011, CCME 2004 ). Increased P loading in waters can have negative effects on aquatic ecosystems and biota. The mechanism leading to these impacts is often increased algal growth due to higher P availability, which can cause eutrofication and hypoxia (CCME, 2004). Although this process is more common in polluted lakes, ponds, and wetlands, it can occur in rivers. Increased N loading can also have negative effects on riverine environments. N and P are both highly mobile nutrients when in solution, so it is expected that they would to some extent both increase under high-flow conditions. Because of their detrimental effects on aquatic ecosystems, it is useful to know which flow conditions threaten water quality and aquatic biota.

There is long-established, nonlinear, positive relationship between N and P flux (Rabalais, 2002). In this report, I test the hypothesis that higher flow conditions (low, moderate, high) are significantly associated with higher P (ug/L) levels. This report also considers covariates Nitrogen (mg/l), Chloride (mg/l), sample location (along river or reach), and weather condition (rainy/clear) as potential predictors of P levels. Ultimately, the question this study aims to answer is the following: which of the above variables can best predict Phosphorus levels in the areas sampled?

The null hypothesis states that there is no significant relationship between sample location, flow condition, weather condition, chloride, or nitrogen and phosphorus levels. In other words, the null hypothesis states that random noise outperforms the above variables in predicting Phosphorus levels. 

# Methods

## Data exploration

```{r, results='hide', message=FALSE, warning=FALSE}
library(tidyverse)
library(dplyr)
library(stats)
library(MASS)
library(psych)
library(AICcmodavg)
library(lattice)
library(HH)
library(corrplot)
library(car)
library(pscl)
library(DHARMa)
library(AER)
library(tinytex)
library(ggplot2)
```


```{r, results='hide', message=FALSE, warning=FALSE}
# Read in the data
noot <- read.csv("data/nutrient_flows_ch_n_p.csv",header=TRUE)
noot <- subset (noot, select = -Parameter.Name) # Remove "Parameter Name" column from data frame
noot <- subset (noot, select = -Sampling.Project) 
noot <- subset (noot, select = -Watershed)
noot <- subset (noot, select = -Sample.Date)
  # Remove "." from Flow.Condition column name CODE noot <- rename(noot, FlowCondition = Flow.Condition)

# Reclassify categorical variables as factors
noot$Flow.Condition <- as.factor((noot$Flow.Condition)) 
noot$Site.ID <- as.factor((noot$Site.ID))               
noot$Weather.Status <- as.factor((noot$Weather.Status)) 

# Remove outlier detected using Cooks D after model fitting
# noot <- noot[-13,]   

# Examine the data

summary(noot)
```
#### Initial inspection

This dataset contains 27 observations which were collected over 4 months, in three sites in NH. It has five fields including:

* __Site.ID__, one location along the Wardsboro River (NH), and two locations on the West River (NH)
* __Chl__, total Chloride (mg/l)
* __Nitro__, total Nitrogen (mg/l)
* __Phosph__, total Phosphorus (μg/L)
* __Flow.Condition__, rate of flow at sample location during sampling (low, moderate, high)
* __Weather.Status__, whether it was raining (wet), and not (dry)

All continuous data are bounded above zero, but Phosphorus was measured in micrograms per liter and, Nitrogen and Chloride were measured in milligrams per liter. Nitrogen and Chloride were scaled for modeling, and added to the original data under the column names Chlz for scaled Chloride and Nitroz for scaled Nitrogen. The below plots show data in the original scales unless otherwise noted. 

A histogram of P levels shows a strong right skew, indicating these data are likely not normally-distributed. 
```{r, message=FALSE, warning=FALSE}
ggplot(data = noot, aes(x = Phosph)) +  
  geom_histogram() +
  labs(title = "Histogram of Phosphorus Levels (μg/L)", y = "Number of Phosphorus Samples", x = "Phosphorus Levels (μg/)")
```
*Figure 1. A right-skewed histogram showing total Phosphorus levels (μg/L).*

The Shapiro-Wilk test of normality returns p-value >.05 (p-value 0.00000002108), confirming that the data are not normally distributed. 

```{r, results='hide'}
shapiro.test(noot$Phosph) # Returns low p-value (.00318) so phosph not normal 
```

The plot below shows one data point (77.4 (μg-P/L)) that seems to stand out as much higher than the rest. Initially, data exploration and model building continued with this data point, but it was eventually removed after testing Cooks Distance. This is covered in the section titled "Fit the Model"
```{r, message=FALSE, warning=FALSE}
library(ggplot2)
par(mfrow=c(1,1))

ggplot(data = noot, aes(x = 1:nrow(noot), y = Phosph)) +  
  geom_point(alpha = .4, col = "steelblue", cex = 3) +
  labs(title = "Total Phosphorus Levels (μg/L)", y = "Phosphorus (μg/L)", x = "Index") 
```
*Figure 2. A scatterplot of Phosphorus levels (μg/L) showing a data point (77.4) that was later removed as an outlier*

The conditional boxplots of Phosphorus levels by river flow condition show a much greater variation during low flow conditions. 
```{r}
ggplot(data = noot, aes(x = Flow.Condition, y = Phosph)) +  
  geom_boxplot(alpha = 1, col = "steelblue", cex = 1) +
  labs(title = "Total Phosphorus Levels (μg/L) by Flow Condition", x = "Flow Condition", y = "P (μg/L)") 
```
*Figure 3. A boxplot of Phosphorus (P) levels (μg/L) organized by Flow Condition (low, moderate, high), showing the relative variation among P levels across different flow conditions.*


```{r, results='hide'}
noot$Chlz     <- c(scale(noot$Chl))    #Scale all continuous variables because of differences in units  
noot$Nitroz   <- c(scale(noot$Nitro))
noot$Phosphz  <- c(scale(noot$Phosph))
```

After scaling continuous variables, a scatter plot of Nitrogen and Phosphorus shows what appears to be a positive relationship, with the exception of the probable outlier. Although informal and insufficient, this supports the prediction that Phosphorus and Nitrogen levels would fluctuate together. 
```{r}
ggplot(data = noot, aes(x = Nitroz, y = Phosphz)) +
  geom_point(alpha = .4, col = "steelblue", cex = 3) +
   labs(title = "Total Phosphorus Levels by Nitrogen levels", x = "Phosphorus (scaled)", y = "Nitrogen (scaled)") 
```
*Figure 4. A scatterplot of Phosphorus (scale) by Nitrogen (scaled) showing a potential slight positive association.* 


## Hypotheses

#### Null hypothesis

As stated in the introduction, the null Hypothesis of this study is that there is no significant association between Site.ID, Chloride, Nitrogen, Flow Condition, or Weather Status, and P levels. The statistical formulation of the null hypothesis is as follows 

Null hypothesis $$y_i = 𝛽_{Phosph} + e_i $$
Where:
$$y_i$$ is the y intercept of response variable (P) of the model
$$𝛽_{Phosph}$$ is the mean P
$$e_i$$ is the error component of the model

Importantly, if and only if the above model receives a lower AIC score than any of the alternative hypotheses below, it will be considered the best model, and none of the independent variables above will be considered as significantly related to P levels. 

### Alternative Hypothesis

Many models were tested in the attempt to discover association between P and dependent variables. Not all models will be listed here, but the top three performing alternative hypothesis, and their associated models are presented below. 

The first alternative hypothesis is that N levels and high flow conditions are positively associated with P levels. Considering that there is a wide body of evidence supporting this hypothesis, a model representing this hypothesis is the focus of this study. All other hypotheses are tested as permutations of the variables in the global model, which is appears in the following section. For the sake of brevity, details on only the top three models are included in this report. 

### Model Description

#### Choosing an Error/Residuals Distribution

This dataset includes 3 categorical variables: Site.ID (3 levels), Flow.Condition (3 levels), and Wather.Status (binary). It also includes 3 positive, coninuous variables that are all bounded above 0. A Generalized Linear Model (GLM) is capable of modeling these variales, so I created a GLM and with an inverse gaussian distribution. The model for the first alternative hypothesis that flow conditions and N levels are significantly related to P levels appears below. 

#### Model Structure

Linear Predictor                                $$n_i =𝛽_{0} + 𝛽_{Flow}X_{P_i} + 𝛽_{Nitro_i}X_{P_i}$$

Error/residuals distribution                    $$e_i~InverseGaussian(ŷ,σ)$$            

Link Function                                   $$μ = n = 1/𝛽_{0} + 𝛽_{Flow}X_{P_i} + 𝛽_{Nitro_i}X_{P_i}$$

Therefore                                       $$y_i = 𝛽_{0} + 𝛽_{Flow}X_{P_i} + 𝛽_{Nitro_i}X_{P_i} + e_i  $$

Define the estimates:


   $$𝛽_{0}$$ is mean P when Flow.Condition = High and Nitro = .11 (mg/L)
                                          
   $$𝛽_{Flow}$$ is the slope coefficient showing the effect on P of a 1 unit change Flow.Condition 

   $$𝛽_{Nitro_i}$$ is the slope coefficient showing the effect on P of a 1 unit change in Nitro

#### Model Assumptions
The residuals have an inverse gaussian distribution. The response data is positive, continuous, and bounded above 0. Predictors are measured without error. This data meets all of these assumptions to the best of my knowledge at the time of authorship. 

## Model Fitting

Model fitting was conducted in R Studio Version 1.4.1717. Over the course of model fitting, I used the R packages "MASS" (Venables and Ripley, 2002) and "STATS" (R Core Team, 2021). I fit 18 models in total. Because I expectat the model with predictor variables of flow condition and nitrogen to outperform all others, that model (the first alternative hypothesis), the null model, and my second alternative hypothesis (predicting P with N only) are presented below. The global model is also included. The remaining 15 models share the same format, but are different permutations of the predictor variables. Each model was assigned to a variable with a naming format of "m[model number].[first letter(s) of variable name(s)]".

```{r, echo=TRUE, results='hide'}
m0        <- glm(Phosph ~ 1, 
                 data = noot, 
                 family = gaussian(link = "inverse"))
m9.fn     <- glm(Phosph ~ Flow.Condition + Nitroz, 
                 data = noot, 
                 family = gaussian(link = "inverse"))
m18.n     <- glm(Phosph ~ Nitroz, 
                 data = noot, 
                 family = gaussian(link = "inverse"))
m7.full   <- glm(Phosph ~ Flow.Condition + 
                   Site.ID + Chlz + 
                   Nitroz + Weather.Status, 
                 data = noot, 
                 family = gaussian(link = "inverse"))

```


# Preliminary Results

#### Output evaluation (before outlier removal)

```{r, results='hide'}
# MODEL FITTING
m0        <- glm(Phosph ~ 1, data = noot, family = gaussian(link = "inverse"))

m1.f      <- glm(Phosph ~ Flow.Condition, data = noot, family = gaussian(link = "inverse"))

m2.s      <- glm(Phosph ~ Site.ID, data = noot, family = gaussian(link = "inverse"))

m3.w      <- glm(Phosph ~ Weather.Status, data = noot, family = gaussian(link = "inverse"))

m4.fs     <- glm(Phosph ~ Flow.Condition + Site.ID, data = noot, family = gaussian(link = "inverse"))

m5.fsc    <- glm(Phosph ~ Flow.Condition + Site.ID + Chlz, data = noot, family = gaussian(link = "inverse")) 

m6.fscn   <- glm(Phosph ~ Flow.Condition + Site.ID + Chlz + Nitroz, data = noot, family = gaussian(link = "inverse"))

m7.full   <- glm(Phosph ~ Flow.Condition + Site.ID + Chlz + Nitroz + Weather.Status, data = noot, family = gaussian(link = "inverse"))

m8.fc     <- glm(Phosph ~ Flow.Condition + Chlz, data = noot, family = gaussian(link = "inverse"))

m9.fn     <- glm(Phosph ~ Flow.Condition + Nitroz, data = noot, family = gaussian(link = "inverse"))

m10.fw    <- glm(Phosph ~ Flow.Condition + Weather.Status, data = noot, family = gaussian(link = "inverse"))

m11.sc    <- glm(Phosph ~ Site.ID + Chlz, data = noot, family = gaussian(link = "inverse"))

m12.sn    <- glm(Phosph ~ Site.ID + Nitroz, data = noot, family = gaussian(link = "inverse"))

m13.sw    <- glm(Phosph ~ Site.ID + Weather.Status, data = noot, family = gaussian(link = "inverse"))

m14.swn   <- glm(Phosph ~ Site.ID + Weather.Status + Nitroz, data = noot, family = gaussian(link = "inverse"))

m15.swnc  <- glm(Phosph ~ Site.ID + Weather.Status + Nitroz + Chlz, data = noot, family = gaussian(link = "inverse"))

m16.wc    <- glm(Phosph ~ Weather.Status + Chlz, data = noot, family = gaussian(link = "inverse"))

m17.wn    <- glm(Phosph ~ Weather.Status + Chlz + Nitroz, data = noot, family = gaussian(link = "inverse"))

m18.n     <- glm(Phosph ~ Nitroz, data = noot, family = gaussian(link = "inverse"))


```

```{r, results='hide'}
m.list <- list("m0"       = m0,  
               "m1.f"     = m1.f,
               "m2.s"     = m2.s,
               "m3.w"     = m3.w, 
               "m4.fs"    = m4.fs,
               "m5.fsc"   = m5.fsc,                                                   # BEFORE OUTLIER REMOVAL 
               "m6.fscn"  = m6.fscn, # After outlier removal                         # This model has lowest AIC score 142.10
               "m7.full"  = m7.full, # 5th Lowest Delta AIC = 7.28                   # Lowest Delta AIC                  4.43
               "m8.fc"    = m8.fc,   
               "m9.fn"    = m9.fn,   # LOWEST AIC SCORE = 131.67
               "m10.fw"   = m10.fw, 
               "m11.sc"   = m11.sc,
               "m12.sn"   = m12.sn,  # Third Lowest Delta AIC = 2.88
               "m13.sw"   = m13.sw,
               "m14.swn"  = m14.swn,                                                  # Next lowest Delta AIC             14.09
               "m15.swnc" = m15.swnc,
               "m16.wc"   = m16.wc, 
               "m17.wn"   = m17.wn, 
               "m18.n"    = m18.n)   #Second Lowest DELTA AIC = 2.15
```

The Akaike information criterion (AIC), a  metric commonly used to rank and/or compare the fit of regression models, was used for model selection. Models from a list were processed and organized into an AIC table using the R package "AICcmodavg" (Mazerolle, 2020). Before removing any outliers, model 6, which includes Flow.Condition, Site.ID, Chlz, and Nitroz as predictor variables had the lowest AIC score. The AIC table below shows the scores for each model. 

```{r, warning=FALSE, message=FALSE}
library(AICcmodavg)
aictab(m.list)
```

I used a variance inflation factor, from the R package "car" (Fox and Weisberg, 2019), to check the relative power of each variable in explaining variation of P levels. According to the output values of VIF, N (1.447) had the greatest ability to explain variation in P, followed by Site.ID (2.024), Flow.Condition (2.102), and Chlz (3.458944). 
```{r, warning=FALSE}
library(car)
vif(m6.fscn)
```

I assessed Cooks Distance, using the R package "car" (Fox and Weisberg, 2019) Cooks Distance quantifies the distance of residuals from the estimated mean of the response variable given a paticular model. The residuals vs. leverage plot of m6.fscn shows data points crossing Cooks Distance, indicating poor model fit and that outliers may be present.

```{r, message = FALSE, warning=FALSE}
plot(m6.fscn, which = c(5)) # Add labels
```
*Figure 5. The residuals vs. leverage plot showing data points beyond Cooks Distance, suggesting the prescence of outliers*

The bar graph below shows that observation 13 has a Cooks distance of 353.023. The mean Cooks distance for *m6.fscn* is 14.047. Because 353.023 > 28.148 (twice the mean Cooks D), this observation was removed from the dataset. The value represented by this sample to be removed is 77.4 (μg-P/L). This level of P is within the range of variation for freshwater rivers, and is biologically significant. This sample's status as an outlier (for the sake of this model) is likely the result of flow condition being recorded categorically and subjectively. If continuous, stream gauge data from the sample locations was added to this dataset, better predictive models could be created in the future. 
```{r}
plot(m6.fscn, pch = 16, col = "darkblue", which = c(4))
```
*Figure 6. A bar graph showing Cooks Distance for each observation*

Summary statistics of Cooks Distance for m6.fscn (the model that includes Flow.Condition, Site.ID, Chlz, and Nitro as predictor variables) are presented below. Visual inspection of the Cooks Distance plot and low first quartile, median, and third quartile values indicate overall low dispersion, with the exception of one outlier.   
```{r}
cooksd <- cooks.distance(m6.fscn) #How to interpret ?
cooksd <- round(cooksd, 5)    
summary(cooksd) 
```

```{r, results='hide'}
   # With the below exceptions, all values are below 1 (the below values were for m6, before P=77 was removed)
rev (sort(cooksd))
```



```{r message=FALSE, warning=FALSE, results='hide'}
# Read in the data
noot <- read.csv("data/nutrient_flows_ch_n_p.csv",header=TRUE)
noot <- subset (noot, select = -Parameter.Name) # Remove "Parameter Name" column from data frame
noot <- subset (noot, select = -Sampling.Project) 
noot <- subset (noot, select = -Watershed)
noot <- subset (noot, select = -Sample.Date)

#Scale all continuous variables because of differences in units 
noot$Chlz     <- c(scale(noot$Chl))     
noot$Nitroz   <- c(scale(noot$Nitro))
noot$Phosphz  <- c(scale(noot$Phosph))

# Reclassify categorical variables as factors
noot$Flow.Condition <- as.factor((noot$Flow.Condition)) 
noot$Site.ID <- as.factor((noot$Site.ID))               
noot$Weather.Status <- as.factor((noot$Weather.Status)) 

# Remove outlier detected using Cooks D after model fitting
 noot <- noot[-13,]   

 
 
# Examine the data
summary(noot)
```


```{r, results='hide'}
m0        <- glm(Phosph ~ 1, data = noot, family = gaussian(link = "inverse"))

m1.f      <- glm(Phosph ~ Flow.Condition, data = noot, family = gaussian(link = "inverse"))

m2.s      <- glm(Phosph ~ Site.ID, data = noot, family = gaussian(link = "inverse"))

m3.w      <- glm(Phosph ~ Weather.Status, data = noot, family = gaussian(link = "inverse"))

m4.fs     <- glm(Phosph ~ Flow.Condition + Site.ID, data = noot, family = gaussian(link = "inverse"))

m5.fsc    <- glm(Phosph ~ Flow.Condition + Site.ID + Chlz, data = noot, family = gaussian(link = "inverse")) 

m6.fscn   <- glm(Phosph ~ Flow.Condition + Site.ID + Chlz + Nitroz, data = noot, family = gaussian(link = "inverse"))

m7.full   <- glm(Phosph ~ Flow.Condition + Site.ID + Chlz + Nitroz + Weather.Status, data = noot, family = gaussian(link = "inverse"))

m8.fc     <- glm(Phosph ~ Flow.Condition + Chlz, data = noot, family = gaussian(link = "inverse"))

m9.fn     <- glm(Phosph ~ Flow.Condition + Nitroz, data = noot, family = gaussian(link = "inverse"))

m10.fw    <- glm(Phosph ~ Flow.Condition + Weather.Status, data = noot, family = gaussian(link = "inverse"))

m11.sc    <- glm(Phosph ~ Site.ID + Chlz, data = noot, family = gaussian(link = "inverse"))

m12.sn    <- glm(Phosph ~ Site.ID + Nitroz, data = noot, family = gaussian(link = "inverse"))

m13.sw    <- glm(Phosph ~ Site.ID + Weather.Status, data = noot, family = gaussian(link = "inverse"))

m14.swn   <- glm(Phosph ~ Site.ID + Weather.Status + Nitroz, data = noot, family = gaussian(link = "inverse"))

m15.swnc  <- glm(Phosph ~ Site.ID + Weather.Status + Nitroz + Chlz, data = noot, family = gaussian(link = "inverse"))

m16.wc    <- glm(Phosph ~ Weather.Status + Chlz, data = noot, family = gaussian(link = "inverse"))

m17.wn    <- glm(Phosph ~ Weather.Status + Chlz + Nitroz, data = noot, family = gaussian(link = "inverse"))

m18.n     <- glm(Phosph ~ Nitroz, data = noot, family = gaussian(link = "inverse"))


```

```{r, results='hide'}
m.list <- list("m0"       = m0,  
               "m1.f"     = m1.f,
               "m2.s"     = m2.s,
               "m3.w"     = m3.w, 
               "m4.fs"    = m4.fs,
               "m5.fsc"   = m5.fsc,                                                   # BEFORE OUTLIER REMOVAL 
               "m6.fscn"  = m6.fscn, # After outlier removal                         # This model has lowest AIC score 142.10
               "m7.full"  = m7.full, # 5th Lowest Delta AIC = 7.28                   # Lowest Delta AIC                  4.43
               "m8.fc"    = m8.fc,   
               "m9.fn"    = m9.fn,   # LOWEST AIC SCORE = 131.67
               "m10.fw"   = m10.fw, 
               "m11.sc"   = m11.sc,
               "m12.sn"   = m12.sn,  # Third Lowest Delta AIC = 2.88
               "m13.sw"   = m13.sw,
               "m14.swn"  = m14.swn,                                                  # Next lowest Delta AIC             14.09
               "m15.swnc" = m15.swnc,
               "m16.wc"   = m16.wc, 
               "m17.wn"   = m17.wn, 
               "m18.n"    = m18.n)   #Second Lowest DELTA AIC = 2.15
```


# Results

#### Evaluating Output: Model Selection using AIC 

After removing the outlier, models were re-run, assembled into a list, and organized into a table based on their AIC scores. Model m9.fn recieved the lowest AIC score, followed by models m18.n and m12.sn.The Delta AIC scores for m18.n and m12.sn are less than 10 (as are the next 6 models). For the sake of brevity, only the top model is evaluated and interpreted.

**AIC Table**
```{r, warning=FALSE, message=FALSE}

aictab(m.list)
```

After outlier removal, the summary statistics of Cooks Distance for the model m9.fn show a maximum reduced by 335.754, down from the original 353.023 to the current maximum, 17.269

```{r}
cooksd <- cooks.distance(m9.fn) #How to interpret ?
cooksd <- round(cooksd, 5)    
summary(cooksd)     # With the below exceptions, all values are below 1 (the below values were for m6, before P=77 was removed)
```


##### Model Validation using a Variance Inflation Factor

I used a variance inflation factor again to evaluate the relative power of predictor variables in explaining P variability. According to the VIF output for the top model, m9.fn, GVIF values were 1.11 for Flow.Condition and Nitroz. Because these values are <2, the variables will be kept as predictors. These low VIF values combined with a lower Cooks Distance indicate that this model likely has adequate fit. 

```{r, warning=FALSE, message=FALSE}
library(car)
vif(m9.fn)
```

VIF values for the third best model (Delta AIC < 3 compared to top model) were 2.07 for both Site.ID and Nitro, indicating poorer model fit than the top model. At this point, m9.fn was selected as the best model, re-assigned to a variable called "mbest", used in predictions, and interpreted.

```{r, warning=FALSE, message=FALSE}
vif(m12.sn)
```

#### Interpreting results

```{r, results='hide'}
(mean.chl <- mean(noot$Chl))      #Chloride Summary Stats
(var.chl <- var(noot$Chl))
(sd.chl <- sd(noot$Chl))

(mean.nitro <- mean(noot$Nitro))  #Nitrogen Summary Stats
(var.nitro <- var(noot$Nitro))
(sd.nitro <- sd(noot$Nitro))

(mean.pho <- mean(noot$Phosph))   #Phosphorus Summary Stats
(var.pho <- var(noot$Phosph))
(sd.pho <- sd(noot$Phosph))
```

Summary statistics for the top model, m9.fn, are shown below. 

```{r}
summary(m9.fn)
```

Because nitrogen was scaled before modeling, I back-transformed the Nitrogen estimate of -0.024005 for interpretation. The resulting estimate was 0.2411993. 

```{r, echo=TRUE}
(-0.024005*sd.nitro)+mean.nitro
```


First, we should note that the p-values for Flow.ConditionHigh, Flow.ConditionModerate, and Nitroz, are all below 0.05, indicating that they have a significant relationship with P levels at the sample locations. The p-value for Flow.ConditionLow was 0.05702, which is less than the significance level of .05, so low flow conditions are not significantly related to P levels at the sample locations, according to this model. 

As a reference, the base case is when Flow.Condition is "high" and when Nitrogen levels are at .11 (mg/L). At "high" flow condition, the model coefficient table estimates that mean P levels (the y-intercept) will increase by 0.077322 μg-/L. At moderate levels, the mean  increases by 0.017766. It is notable that the p-value for the "high" flow rate estimate was significantly lower than the p-value for the "moderate" flow rate estimate, indicating that "high" flow rate has a more significant relationship with P levels than moderate. That said, they are both significant. The back transformed estimate for nitrogen indicates thatP levels increase by 0.2373459 μg-/L for a 1 unit increase in Nitrogen. 


#### Predictions

```{r, results='hide'}
#re-assign model to new variable name (not sure why, but this solved problems when first attempting pedictions)
mbest = glm(Phosph ~ Flow.Condition + Nitroz, data = noot, family = gaussian(link = "inverse"))

#Expand grid
noot.newdata <- expand.grid(seq(min(noot$Nitroz),max(noot$Nitroz), .1), levels(noot$Flow.Condition))
colnames(noot.newdata) <- c("Nitroz","Flow.Condition")

# Create prediction
pred.1 <- predict(mbest,noot.newdata,type="response", se.fit=T)


noot.newdata$n.se <- pred.1$se.fit
noot.newdata$n.fit <- pred.1$fit
noot.newdata$resid.scale <- pred.1$residual.scale
```


```{r}
#Return Nitrogen to original scale
noot.newdata$Nitros <- (noot.newdata$Nitroz*sd.nitro)+mean.nitro
```

The regression plot was created using the predict function from the R package "stats". Predictions were created for the range of values of Nitrogen, which was back transformed for plot-labeling. The prediction plot shows a positive relationship between increased N and the response variable, P. Even though Flow.ConditionLow did not have a significant relationship with P levels, this plot indicates that on average, higher flow conditions ehibit greater P levels than moderate flow conditions. 
```{r}
#Plot Prediction with Original Nitrogen Scale (Nitrogen is in mg/l)

ggplot(noot.newdata, aes(x = Nitros, y = n.fit, color = Flow.Condition)) +
  geom_line() +
  labs(title = "Predicted Phosphorus (μg/L)", 
       y = "Predicted Phosphorus Concentration (μg/L)",
       x = "Nitrogen Concentration (mg/L)") +
  theme_classic()
```
*Figure 7. Regression plot for predicted N and P under 3 flow rate scenarios*

Predictive model validation (for m9.fn) was completed using the R package, "DHARMa" (Hartig, 2022). The QQ plot of residuals in the DHARMa residuals plot below (left) shows moderate dispersion of residuals, deviation is not significant, and the Kolmogorov Smirnov Test indicates that the model fits well. 

```{r, warning=FALSE, results='hide', message=FALSE}
library(DHARMa)
prediction_resid <- simulateResiduals(fittedModel = mbest, plot = F)
```

```{r}
#Plot Residuals KS test
plot(prediction_resid)
```
*Figure 8. DHARMa residuals plot for predictive model (m9.fn) showing observed and predicted residuals*


# Discussion

18 models with an inverse gaussian distribution were fit to test the alternative hypothesis that flow rate and N levels are significantly related to P levels. Considering the wide body of evidence supporting the positive relationship between flow rates, N and P levels, it is consistent with landmark ecological and hydrological findings that a model representing flow rate and N as independent variables had the greatest ability to predict P levels. Nitrogen and high flow conditions had the most significant relationship to P levels. The output from this modeling process would likely have been different had flow condition be measured as a continuous variable with precise measurements (to the cubic foot per second). Because nitrogen and phosphorus can both be moved in sediment, and in suspension, and can both be taken up by vegetation, this assessment suggests that terrestrial and riparian restoration and pollutant-mitigation tactics focused on re-vegetation (for increased nitrogen uptake) may also be effective for phosphorus. Additionally, findings from this assessment suggest that terrestrial and riverine/aquatic restoration measures aimed at reducing overland flow and in-stream flowrates have the potential to lower P levels. 










#### References 
       
W. N. Venables and B. D. Ripley, Modern Applied Statistics with S 4th Edition, Springer, New York, 2002, ISBN 0-387-95457-0, https://www.stats.ox.ac.uk/pub/MASS4/

R Core Team, R: A Language and Environment for Statistical Computing, R Foundation for Statistical Computing, Vienna, Austria, 2021, https://www.R-project.org/
  
J Fox, S. Weisberg, An {R} Companion to Applied Regression 3rd edition, 2019, Sage, Thousand Oaks, CA, https://socialsciences.mcmaster.ca/jfox/Books/Companion/
 
M. J. Mazerolle, AICcmodavg: Model selection and multimodel inference based on (Q)AIC(c), 2020, R package version 2.3-1, https://cran.r-project.org/package=AICcmodavg

F. Hartig, DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models, 2022, R package version 0.4.5, https://CRAN.R-project.org/package=DHARMa

Marklein, A. R., & Houlton, B. Z. (2011). Nitrogen inputs accelerate phosphorus cycling rates across a wide variety of terrestrial ecosystems. New Phytologist, 193(3), 696–704. https://doi.org/10.1111/j.1469-8137.2011.03967.x

Rabalais, Nancy N. “Nitrogen in Aquatic Ecosystems.” Ambio, vol. 31, no. 2, 2002, pp. 102–12, http://www.jstor.org/stable/4315222. Accessed 5 May 2022.

Canadian Council of Ministers of the Environment. 2004. Canadian water quality guidelines for the protection of aquatic life: Phosphorus: Canadian Guidance Framework for the Management of Freshwater Systems. In: Canadian environmental quality guidelines, 2004, Canadian Council of Ministers of the Environment, Winnipeg.

Authored by John Slaff, May 5th, 2022

    
    