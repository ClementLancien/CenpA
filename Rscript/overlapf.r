#Check if their is a significant difference between SRR633614_SRR633612 
#and SRR633615_SRR633613 observed_ration


#################################
#           Load DATA           #
#################################

##Load SRR633614_SRR633612 data after overlapf (cf pipeline.sh)

#data1 -> SRR633614 SRR633612
data1 = read.table("../Result/Align/test.txt", header=TRUE, sep=',')

#data2 -> SRR633615 SRR633613
data2 = read.table("../Result/Align/test2.txt", header=TRUE, sep=',')

#Check data1
summary(data1)
str(data1)

#Check data2
summary(data1)
str(data1)

ratio1 = data1$observed_ratio / data1$theorical_ratio

####################################
#           SHAPIRO TEST           #           
####################################

#Our data are regarded as unpaired. It means there is no possibility of
#the values in one data set being related to or being influenced by the value
#in the other data sets.

#For numerical data, it is important to decide if they follow the 
#parameters of the normal distribution curve (Gaussian curve), 
#in which case parametric tests are applied (else use non-parametric tests).

#The Shapiro–Wilk test tests the null hypothesis that a sample x1, ..., xn came from a normally distributed population.
#(The null-hypothesis of this test is that the population is normally distributed.)
# Thus, if the p-value is less than the chosen alpha level, then the null hypothesis is rejected and 
# there is evidence that the data tested are not from a normally distributed population.
# On the contrary, if the p-value is greater than the chosen alpha level, 
# then the null hypothesis that the data came from a normally distributed population cannot be rejected.
# p <= 0.01 : tres forte présomption contre l'hypothèse nulle
# 0.01 < p <= : forte présomption contre l'hypothèse nulle
# 0.05 < p <= 0,1 : faible présomption contre l'hypothèse nulle
# 1 ^> 0,1 : pas de présemption contre l'hypothèse nulle

#data1
shapiro.test(data1$observed_ratio)
# p-value = 0.2891
# data1 is normally distributed

#data1
shapiro.test(data1$observed_ratio)
# p-value = 0.3322
# data1 is normally distributed

###############################
#           Pearson           #
###############################

cor(data1$observed_ratio, data1$theorical_ratio, method='pearson')
cor(data2$observed_ratio, data2$theorical_ratio, method='pearson')
