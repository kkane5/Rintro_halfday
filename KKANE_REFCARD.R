# This is Kevin Kane's R Reference Card
# Updated 9.3.2013


#########################################
#### FILES AND DIRECTORIES ##############
#########################################

# Set/get working directory
getwd()
setwd("C:\\Users\\Kevin\\Documents\\ASU\\R and GIS")

# This lets you use DBF files
library(foreign)

data = read.dta("data.dta")
write.csv(data, file = "data2.csv")
data = read.csv(file.choose(),header=T)
read.csv(text=readLines('myfile.csv')[-(1:5)])   # omits first 5 lines of file BEFORE reading it 


ls()			# list current objects
rm(object)		# remove object
class(object)	# outputs the type of the object
typeof(object)	# same as above
cbind()		# combines objects by column
rbind()		# combines objects by row
smartbind()		# uses library(gtools) to bind by row - even when names differ

# Remove all but one object from memory... replace QQQ with object name
rm(list=ls()[!grepl('^QQQ$',ls())]) 

###############################################
#### BASIC DATA OPERATIONS ####################
###############################################

#### Making your own data, then adding to it, then transpose, sort, extract ####
data1 = c(1, 2, 3, 4, 5, 6, 7, 8)
data2 = c(43,23,33,64,-1111.1,32, 44, 121)
data3 = c(32,43,-4,323,43,99, 102, 333)
names = c("foo","bar","baz","yak", "ok","len", "bah", "huh")
d = data.frame(data1, data2, data3, row.names=names)

#### OPERATIONS ON A VECTOR ####
first_column = d$data1		# extracts a vector from a data frame
sort(data2)		#this only works for a vector - not a data frame
e = d$data2[1:4]		# extracts first four observations of data2
e = d[1:4,2]		# same as above
d[d$data2==-1111.1,"data2"] = NA	# set variables with a missing value code (-1111.1 in this case) to missing value code in R (NA)
d[d$data3<0,"data3"] = NA		# set variables with a missing value code (less than 0) to R's mising value (NA)

#### OPERATIONS ON A DATA FRAME (MATRIX) ####
d[d<0]=NA		# label all below-zero values in data frame as NA
d2 = na.exclude(d)	# make a new data frame excluding observations which have NA values
d[order(-d$data2), , drop = FALSE]		#sort data frame in descending order by data2 BUT DOES NOT SAVE IT
az3 = az2[order(az2$percentile),]		# another way to order by row
d2 = as.data.frame(aperm(as.matrix(data)))# transpose rows and columns of data frame
d3 = d[1:3,]	# Extract only the first three rows of d
d3 = d[-2,]		# Remove the second observation
d3 = d[-nrow(d),]	# Remove the last observation
d3 = d[(-nrow(d)+2):-nrow(d),]	# Remove the last three observations

# Subset commands (for a data frame) (see TIF)
chi_1997 = subset(mydata, city=="CHICAGO")	# subset of observations where variable city is Chicago
nonChi_1997 = subset(mydata, city!="CHICAGO")	# subset of observations where variable city is NOT Chicago
data2 = split(data,data$gender)			# split data into separate datasets by gender
write.csv(chi_1997, file = "chi_1997.csv")	# write this to a new .csv

# row and column names can be extracted too
> rownames(mydata)
> colnames(mydata)

# Number of rows and columns
> nrow(mydata)
> ncol(mydata)

# Other Useful Data Exploration Tools
str(mydata)			# storage type of each variable
dim(mydata)			# dimension of your data

# Manipulating Variables within a data frame
myvars <- c("pin", "totval", "city", "latitude", "longitud")
newdata <- mydata[myvars]
d2 = subset(d, select=c(data1,data2))	# alternative way to do this using subset
d3 = d[1:2]						# select only the first two columns in data frame d
data2 = data[sapply(data,is.numeric)]    	# makes a new dataset of only the numeric variables

# Coerce variables of unknown type to type-double (see Ancient_Cities):
for(i in 1:ncol(data)){
  data[,i] = as.double(as.character(data[,i]))
}

# Select columns of data frame with ONLY ONE missing value (see Ancient_Cities): 
myvars = NULL
for(i in 1:ncol(data)){
x = data[,i]
if (length(na.exclude(x))>=length(x)-1)
   myvars = append(myvars, colnames(data)[i])
}
d_new = data[myvars]

# create an index, and bind it to an existing data frame (see ECN525 project)
index = seq(1:length(high_gdp2$GDP96))
high_gdp3 = cbind(index,high_gdp2)

# Sum by row entry factor (similar to Stata's collapse or Excel's pivot table) (see TIF)
tif2 = aggregate(tif$Value ~ tif$Entity, data=tif, sum)

#### HOW TO DO A V-LOOKUP IN R #### (see Patents, Rank-Size Paper)
setwd("C://Users//Kevin//Documents//ASU//Patents//RankSize_Paper")
subcats = read.csv("Patent_Subcat_Mapping.csv", header=T)
patents = read.csv("allpatbytotpat2.csv", header=T)
# subcats is the widget... should only have 2 columns - the widget and the new info to put into the main data
# patents is the main data
# column name ('by') must be same - widget and main data
patents2 = merge(subcats, patents, by = "pat_class")

###############################################
#### STATISTICS ###############################
###############################################

summary(data)		# mean, median, etc. of all variables in data frame
head(data)			# gives you the first 6 ish rows of data to check it out
abs(sort(table(n)))	# mode - easy way to get this
mean(data,na.rm=T)	# mean, but after removing NA (missing) values
plot(density(data$var)	# kernel density plot

rnorm()		# random normal distribution
dnorm(1.0)		# pdf at 1.0 standard deviations (pdf = probability density function)
pnorm(1.6)		# cdf up to 1.6 standard deviations (cumulative density function)
quantile(t,0.75) - quantile(t,0.25)		# interquartile range
skewness = sum(((x-mean(x))^3)/(length(x)*sd(x)^3))
kurtosis = sum(((x-mean(x))^4)/(length(x)*sd(x)^4))

# Choose from a binomial distribution
pr = 0.24		# select probability
n = 100		# number of observations
x = 0:n
px = choose(n,x)*pr^x*(1-pr)^(n-x)
plot(x,px)

# draws from uniform distribution [a,b] to illustrate CLT 
a=0  #upper limit of population
b=20 
n=10  #n is the sample size
mu=(b+a)/2  #mu is the mean of the population
sigma=sqrt(((b-a)^2-1)/12) #sigma is the standard deviation of the population
ns=1000      # number of samples
x=1:ns      
for(i in 1:ns){ x[i]=mean(as.integer(runif(n,a,b+1))) } #choose ns samples of size n and calculate the means of the 5000 samples
hist(x,freq=F,breaks=40,border='blue',xlim=c(a,b))    #draw a histogram of the means of the samples
sd(x)
standardError=sigma/sqrt(n)  #calculate the standard error by the central limit theorem
xx=seq(0,20,0.01)  # xx will be used to plot the sampling distribution
fx=dnorm(xx,mu,standardError) #calculate the normal sampling distribution by the central limit theorem
lines(xx,fx,col='red')  #plot the sampling distribution by the central limit theorem

# Correlation
library(Hmisc)
corMat = cor(d, use="pairwise.complete.obs")		# correlation for the entire data frame
cor.test(d$data1,d$data2, method="pearson")		# correlation test - just 2 variables.  Pearson is default.
rcorr(as.matrix(d))			# correlation test - entire data frame
pairs(~var1+var2+var3, data=data)		# gives a matrix of correlation plots

# Basic Statistical Tests
x = c(32.1,34.4,34.9,30.6,38.4,29.4,28.9,32.6,32.9,44.9)
x_male = c(75,5,120,7,10,20,45,10,20,0)
x_female = c(150,75,45,15,30)
n = 10
t.test(x, conf.level=0.99)	# default confidence level is 0.95
t.test(x, mu=30)		# test of whether x's mean is 30
t.test(x_male, x_female, alternative = c("greater"), var.equal=FALSE)	# two-sample t-test
prop.test(38,50)		# test of proportion - is 38/50 greater than 50%?
binom.test(38,50)		# exact version of prop.test

# ANOVA (see CAP/Results/Constr_yr)
contag_half = aov(DevYearLand_csv_CONTAG~factor(Half_Dec), data=phx)
summary(contag_half)
print(model.tables(contag_half,"means"),digits=3)
print(model.tables(contag_half,"effects"),digits=3)
TukeyHSD(contag_half)


#### Cross-tabulation matrix - ####
# 1.) Define Function
crosstab <- function(data, var1, var2){
  ctab = matrix(nrow=2, ncol=2)
  colnames(ctab) = c("Var1 absent", "Var1 present")   
  rownames(ctab) = c("Var2 present", "Var2 absent")   
  ctab[1,1] = nrow(subset(data, var1==0 & var2==1))    # no, yes
  ctab[1,2] = nrow(subset(data, var1==1 & var2==1))    # yes, yes
  ctab[2,1] = nrow(subset(data, var1==0 & var2==0))    # no, no
  ctab[2,2] = nrow(subset(data, var1==1 & var2==0))    # yes, no
  ctab_test = chisq.test(ctab)
  test_result = paste("X-sq=", ctab_test[1], "p-val=", ctab_test[3])
return(ctab)
return(test_result)
}
# 2.) example
crosstab(data, nameofVar1, nameofVar2)
# 3.) Loop through all combinations
com = combn(seq(1:ncol(data)), m=2)
for(i in 1:ncol(com)){
  print(paste("VAR1:", colnames(data)[com[1,i]]))
  print(paste("VAR2:", colnames(data)[com[2,i]]))
  print(crosstab(data, data[,com[1,i]], data[,com[2,i]]))
}






###############################################
#### GEOGRAPHIC STATISTICS ####################
###############################################
x = c(60,45,70,55,65,70,80,45,30,55,70,0)
y = c(80,45,60,60,75,45,60,75,70,50,65,40)
w = c(4,5,6,7,4,3,2,2,2,1,1,1)

plot(x,y,cex=w)		# Plot, scaled by the weights (w)
points(sum(x)/length(x),sum(y)/length(y), col = "red")	# add mean center
points(sum(x*w)/sum(w),sum(y*w)/sum(w),col="blue")		# add weighted mean center
abline(h=median(y),lty=3)	# add Manhattan Median (1/2)
abline(v=median(x),lty=3)	# add Manhattan Median (2/2)
standard_distance = sqrt(sum((x-mean(x))^2)/length(x)+sum((y-mean(y))^2)/length(y))
points(mean(x),mean(y),cex=standard_distance, col="goldenrod")	# add standard distance


###############################################
#### DATA DISPLAY #############################
###############################################

### Simple histogram - do the breaks in a couple of ways, making sure when you define a break sequence that it includes all your data points
h=hist(comm)
hist(comm,breaks=20)
br=seq(0,60,5)
h=hist(comm,breaks=br)
h1 = hist(d1, xlim = c(min(d1),max(d1)), breaks= round((max(d1)-min(d1))/10, digits = 0))		#breaks from min to max, by 10, rounded

### Quantile Plot
sortComm=sort(comm)
rankComm=rank(sortComm)
pComm=rankComm/(length(comm)+1)
plot(sortComm,pComm)

### Boxplot
boxplot(comm)
boxplot(comm ~ commute$gender)	# Boxplot, split by gender
barplot(table(census$gender,census$education),legend=c("male","female"),beside=T)	# Boxplot, split by gender, 2 beside each other

### Fancy Boxplot with counts and trendline (see CAP_GRAD/RESULTS/Constr_yr)
a = boxplot.n(DevYearLand_csv_FRAC_AM ~ Half_Dec, ylab="Fractal Dimension (FRAC_AM)", main="Metrics for Block Groups by Construction Year",data=phx, col=c("gold", "darkgreen"), names, shrink=0.7, top=T)
lines(a$stats[c(3),], lwd=3, col="red")


###############################################
#### GRAPHICS #################################
###############################################

#### GRAPH A CIRCLE - THETA IS IN RADIANS
r = 3.5
x = seq(1,20)
angle = seq(0,2*pi,0.1)
plot(r*cos(angle),r*sin(angle), cex=x)
lines(r*cos(angle),r*sin(angle))


#### A BASIC PLOT
y <- c(1915, 1949,1963)
x <- c(58, 42, 30)
x1 <- c(7, 20, 22)
x2 <- c(21, 10, 20)
x3 <- c(12, 29, 27)
plot(y, x, type = "o", col = "cyan2", pch=0, lwd=3, lty=1, ylim = c(0,70), xlim=c(1915,1963), axes = FALSE, xlab="Year", ylab="Percent of Parcels")
lines(y, x1, type = "o", col = "darkgoldenrod1", pch=1, lwd=3, lty=2)
lines(y, x2, type = "o", col = "darkolivegreen1", pch=2, lwd=3, lty=3)
lines(y, x3, type = "o", col = "maroon", pch=3, lwd=3, lty=4)
axis(1, at=c(1915, 1949, 1963), lab=c(1915, 1949, 1963))
axis(2)
box()
abline(h=c(10, 20, 30, 40, 50, 60, 70), lwd=0.5, lty=3, col="gray50")
title(main="Phoenix Land Use Categories", font.main=4)
legend("topright", inset=0.045, c("R", "C", "V", "N"), pch=0:3, lty=1:4, lwd=3, col=c("cyan2","darkgoldenrod1","darkolivegreen1","maroon"))
title(xlab="Year")
title(ylab="Percent of Parcels")

#### 3-D PLOTTING
library(scatterplot3d)
scatterplot3d(x,x1,y)

#### MULTIPLE PLOTS IN SAME IMAGE EDITOR
par(mfrow = c(1,2))		# sets image editor as a 1 row, 2 column matrix - then populate it!
layout(matrix(c(1,2), 1, 2, byrow = TRUE))  #another way to do the same thing
plot(y, x, type = "o", col = "cyan2", pch=0, lwd=3, lty=1, ylim = c(0,70), xlim=c(1915,1963), axes = FALSE, xlab="Year", ylab="Percent of Parcels")
axis(1, at=c(1915, 1949, 1963), lab=c(1915, 1949, 1963))
axis(2)
plot(y, x1, type = "o", col = "red", ylim = c(0,70), xlim=c(1915,1963), axes = FALSE, xlab="Year", ylab="Percent of Parcels")
box()

#### VARIOUS GRAPHICS PARAMETERS (PAR) (SEE CAPGRANT/CONSTR_YR)
pdf(file="Fig3_Boxplots_Landscape_Metrics_REVISED.pdf",7,14)
par(mfrow=c(4,2))   	 #sets a 2x2 matrix
par(oma=c(0,0,2,0))   	 #set outer margins for entire trellis plot
par(las=3)   		 # sets axis labels perpendicular, vertical, etc.
par(mar=c(3,4,2,2))   	 # sets margins - bottom, left, top, right
par(font.lab=2)		 # font for the labels
mtext("Figure 3: Boxplots of Landscape Metrics", side=3, outer=T, font=2, las=0)    # run this line at the end




#### SAVING GRAPHICS TO PDF
pdf(file="figure1.pdf")
# this is where you execute the code for your graphic - plot below as an example
plot(x,y)
dev.off()



#########################################
#### REGRESSION #########################
#########################################

# TYPICAL REGRESSION
model = lm(y~ x)
model1=lm(y~ (x1 + x2 + x3))
model2 = lm(NET~ LnPAT + CoPAT + ICT + BIO,nl)	#you can put the name of the data frame at the end to avoid extra typing
summary(model)

#extract a variable
beta=summary(model)$coefficients["x","Estimate"]
print(beta)

# Plot a simple regression
plot(y,x, xlab="X Variable", ylab = "Y Variable", main = "Title of Plot")
abline(model, col="blue")
legend("bottomright", inset=0.045, legend=c("R-squared =", round(summary(model)$r.squared,4)), bg="white")

# Residual Plots
layout(matrix(c(1,2,3,4),2,2))
plot(model)

# Diagnostics
library(car)
vif(model1)


#########################################
#### JOIN DATA INTO A SHAPEFILE'S DBF####
#########################################
library(foreign)
library(car)
library(Hmisc)

# join 2006 LU_run4 to shapefile's dbf
# must make sure that both files are in the same order or use below
tracts = read.dbf("tl_2010_36_tract10.dbf") 
tracts2 = tracts[order(tracts$GEOID10),]
work = read.csv("ACS_11_5yr_imputed_work.csv")
work2 = work[order(work$OBJECTID),]
myvars = c("Estimate_Total","Estimate_Not.imputed")
work3 = work2[myvars]
tracts3 = cbind(tracts2,work3)
write.dbf(tracts3,"tl_2010_36_tract10.dbf")

# subset from larger (world) set; put into smaller (country) shapefiles; also exclude some unwanted observations
us = subset(w3,Country=="US")
us1 = subset(us, Region!="US074: Honolulu - Hi")
us2 = subset(us1, Region!="US008: Anchorage - Ak")
write.dbf(us2,file="US.dbf")


#####################################
#### WORD CLOUD #####################
#####################################
library(gdata)
library(wordcloud)
wordcloud(data$tag, data$occurences,scale=c(4,.5), max.words=250, random.order = FALSE, random.color=TRUE)
