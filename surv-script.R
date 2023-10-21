library(survival)
install.packages("readxl")
library(readxl)

################
################
## Question 1 ##
################
################

#######
## a ##
#######

# Read the lymphoma data set
lymphoma <- read.csv("data/lymphoma.csv", header = TRUE)
head(lymphoma)
summary(lymphoma)

# Calculate the Kaplan-Meier estimates
lymphoma.km <- survfit(Surv(time, cens) ~ symp, data = lymphoma)
summary(lymphoma.km)

# Plot the Kaplan-Meier estimates
par(pty = "s")
plot(lymphoma.km, mark.time = TRUE,
     xlab ="Time, t", ylab="Estimated survival function, S(t)",
     lty = c(1, 2))
legend("topright", lty = c(1,2), legend = c("Asymptomatic", "Symptomatic"))

#######
## b ##
#######

# View the summary of the KM estimates
summary(lymphoma.km)

# 95% CI for asymptomatic group
# 0.309        0.690

# 95% CI for symptomatic group
# 0.0156        0.642

#######
## c ##
#######

# due to right censoring we can not draw a conclusion as none of our estimates
# fall below 0.3

# t = 276 for symptomatic group

#######
## d ##
#######

# Calculate the Nelson-Aalen estimates
lymphoma.na <- survfit(Surv(time, cens) ~ symp, data = lymphoma, 
                       type = "fleming")
summary(lymphoma.na)


# Plot the Nelson-Aalen estimates
par(pty = "s")
plot(lymphoma.na, mark.time = TRUE,
     xlab="Time, t", ylab="Survival function, S(t)",
     lty = c(1, 2))
legend("topright", lty = c(1,2), legend = c("Asymptomatic", "Symptomatic"))

# Graphically we can see the estimates are very similar (inspecting the 
# numerical results also shows they are very similar)

#######
## e ##
#######

# View the summary of the NA estimates
summary(lymphoma.na)

# 95% CI for asymptomatic group
# 0.319        0.697


# 95% CI for symptomatic group
# 0.0340        0.621

################
################
## Question 2 ##
################
################

#######
## a ##
#######

# Read the duck data set
duck <- read.table("data/duck.txt", header = TRUE)
head(duck)
summary(duck)

# Model with no covariates
duck.wb1 <- survreg(Surv(Time, Observed) ~ 1, data = duck)
summary(duck.wb1)

# Model with just Age
duck.wb2 <- survreg(Surv(Time, Observed) ~ Age + Weight + Length, data = duck)
summary(duck.wb2)
anova(duck.wb1, duck.wb2)
# Do not use Age

# Model with just Weight
duck.wb3 <- survreg(Surv(Time, Observed) ~ log(Weight), data = duck)
summary(duck.wb3)
anova(duck.wb1, duck.wb3)
# Do not use Weight

# Model with just Length
duck.wb4 <- survreg(Surv(Time, Observed) ~ log(Length), data = duck)
summary(duck.wb4)
anova(duck.wb1, duck.wb4)
# Do not use Length

# Fit a main effects model
duck.wb <- survreg(Surv(Time, Observed) ~ Age + Weight + Length, data = duck)
summary(duck.wb)  

#######
## b ##
#######

scaled <- duck$Time/predict(duck.wb1)

# Calculate an estimate of the survivor function
S <- survfit(Surv(scaled) ~ 1)

plot(log(S$time), log(-log(S$surv)), xlab = "log t", ylab= "log (-log S(t))",
     cex.lab = 1.25)

#######
## c ##
#######

duck.cox0 <- coxph(Surv(Time, Observed) ~ 1, data = duck)

# Age
duck.cox1 <- coxph(Surv(Time, Observed) ~ Age + Weight + Length, data = duck)
summary(duck.cox1)

# Weight
duck.cox2 <- coxph(Surv(Time, Observed) ~ log(Weight), data = duck)
summary(duck.cox2)

# Length
duck.cox3 <- coxph(Surv(Time, Observed) ~ log(Length), data = duck)
summary(duck.cox3)

duck.cox4 <- coxph(Surv(Time, Observed) ~ log(Length) + log(Weight), data = duck)
summary(duck.cox4)
anova(duck.cox2, duck.cox4)

duck.cox <- coxph(Surv(Time, Observed) ~ log(Weight), data = duck)
summary(duck.cox)

#######
## d ## 
#######

duck.S1 <- survfit(duck.cox, newdata = data.frame(Age = 1, 
                                                  Weight = 1000, Length = 250))

alpha <- 1 / duck.wb1$scale
theta <- exp(-duck.wb1$coefficients[1])

plot(duck.S1, conf.int = FALSE, ylab= "Survival function", xlab="Time")
curve(expr = exp(-(theta*x)^alpha), from = 0, to = 60, add = TRUE, lty = 2)
legend("topright", lty = c(1,2), legend = c("Cox", "Weibull"))


################
################
## Question 3 ##
################
################

#######
## a ##
#######

# Read the data
mortality <- read.csv("data/mortality.csv", header = TRUE)
summary(mortality)
head(mortality)

# Find m_x for males
mortality$Male.rate <- mortality$Male.deaths / mortality$Male.exposed

# Find m_x for females
mortality$Female.rate <- mortality$Female.deaths / mortality$Female.exposed

# Plot log m_x
plot(mortality$Age, log(mortality$Male.rate), type = 's', 
     xlab = "Age", ylab = "Crude central mortality rate")
lines(mortality$Age, log(mortality$Female.rate), type = 's', lty = 'dashed')
legend("topleft", legend = c("Male", "Female"), lty = c("solid", "dashed"))

#######
## b ##
#######

###################################
### Constant force of mortality ###
###################################

# Find q_x for males under const. force of mortality
mortality$Male.qconst <- 1 - exp(-mortality$Male.rate)

# Find q_x for females under const. force of mortality 
mortality$Female.qconst <- 1 - exp(-mortality$Female.rate)

## Begin constructing the life table for ##
## males under const. force of mortality ##

mortality$Male_l_x1 <- NA

mortality$Male_l_x1[1] <- 100000

for(x in 2:length(mortality$Age)){
  mortality$Male_l_x1[x] <- mortality$Male_l_x1[x-1] * 
    (1 - mortality$Male.qconst[x-1])
}


## Begin constructing the life table for   ##
## females under const. force of mortality ##

mortality$Female_l_x1 <- NA

mortality$Female_l_x1[1] <- 100000

for(x in 2:length(mortality$Age)){
  mortality$Female_l_x1[x] <- mortality$Female_l_x1[x-1] * 
    (1 - mortality$Female.qconst[x-1])
}

for(x in seq(from = 1, to = length(mortality$Age), by = 5)){
  print(c(str(mortality$Age[x]), 
          str(mortality$Male.qconst[x]),
          str(mortality$Female.qconst[x]),
          str(mortality$Male_l_x1[x]),
          str(mortality$Female_l_x1[x])))
}

############################
### Uniform distribution ###
############################

# Find q_x for males under uniform distribution
mortality$Male.qunif <- mortality$Male.rate / (1 + 0.5 * mortality$Male.rate)

# Find q_x for females under uniform distribution
mortality$Female.qunif <- mortality$Female.rate / (1 + 0.5 * mortality$Female.rate)

## Begin constructing the life table for ##
## males under unfiform distribution     ##

mortality$Male_l_x2 <- NA

mortality$Male_l_x2[1] <- 100000

for(x in 2:length(mortality$Age)){
  mortality$Male_l_x2[x] <- mortality$Male_l_x2[x-1] * 
    (1 - mortality$Male.qunif[x-1])
}

## Begin constructing the life table for   ##
## females under const. force of mortality ##

mortality$Female_l_x2 <- NA

mortality$Female_l_x2[1] <- 100000

for(x in 2:length(mortality$Age)){
  mortality$Female_l_x2[x] <- mortality$Female_l_x2[x-1] * 
    (1 - mortality$Female.qunif[x-1])
}

for(x in seq(from = 1, to = length(mortality$Age), by = 5)){
  print(c(str(mortality$Male.qunif[x]),
          str(mortality$Female.qunif[x]),
          str(mortality$Male_l_x2[x]),
          str(mortality$Female_l_x2[x])))
}

#######
## c ##
#######

# Calculate curtate life expectancy for males at age 60

# Initialise e60
Male_e_60 <- 0

# Calculate value with a loop
for(x in 1:(length(mortality$Age) - 1)){
  Male_e_60 <- Male_e_60 + mortality$Male_l_x2[x + 1]
}

Male_e_60 <- 1/mortality$Male_l_x2[1] * Male_e_60
Male_e_60
Male_e_60_complete <- Male_e_60 + 1/2 # Complete life expectancy 
Male_e_60_complete

# Calculate curtate life expectancy for females at age 60

# Initialise e60
Female_e_60 <- 0

# Calculate value with a loop
for(x in 1:(length(mortality$Age) - 1)){
  Female_e_60 <- Female_e_60 + mortality$Female_l_x2[x + 1]
}

Female_e_60 <- 1/mortality$Female_l_x2[1] * Female_e_60
Female_e_60
Female_e_60_complete <- Female_e_60 + 1/2 # Complete life expectancy 
Female_e_60_complete

#######
## d ##
#######

# Download the ELT17
elt17 <- read_excel("data/ELT17.xlsx")
summary(elt17)
head(elt17)

# Calculate the test statistic for males
mX2 <- sum(((mortality$Male.deaths - 
               (mortality$Male.exposed * elt17$male_mx[1:43]))^2) /
             (mortality$Male.exposed * elt17$male_mx[1:43]))
mX2                        # X^2       = 35.218
qchisq(p = 0.95, df = 43)  # Crit. val = 59.304
1 - pchisq(q = mX2, df = 43) # p-value = 0.795

# For males do not reject null hypothesis i.e. standard rates are a reasonable 
# model for study population

# Calculate the test statistic for females
fX2 <- sum(((mortality$Female.deaths - 
               (mortality$Female.exposed * elt17$female_mx[1:43]))^2) /
             (mortality$Female.exposed * elt17$female_mx[1:43]))
fX2                        # X^2       = 57.727
qchisq(p = 0.95, df = 43)  # Crit. val = 59.304
1 - pchisq(q = fX2, df = 43) # p-value = 0.066

# For females do not reject null hypothesis i.e. standard rates are a reasonable 
# model for study population
