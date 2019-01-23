# This file was created by Kevin Kane for an introductory lecture on R for the 2019 School of Data at SCAGFor attribution information contact me via email at kane@scag.ca.gov

# Some R-esources!
# R website (open source): www.r-project.org 
# Rstudio website (freeware): www.rstudio.com 
# A pretty good cheat-sheet: http://www.gardenersown.co.uk/Education/Lectures/R/   #
# Quick-R: www.statmethods.net 

#########################################
######### (1) BASIC SKILLS IN R #########
#########################################

############### (1a) DATA BASICS #################
# Ways to run code:
# (1) Type, or copy/paste into the Console window's prompt (>)
# (2) Place your cursor on a line, or highlight 1 or more lines, and press "Run" at the top-right of the script window
# (3) Place your cursor on a line, or highlight 1 or more lines, and hit Ctrl+R
# Note that the pound/hashtag converts a line of code to a comment

print("hello world")   # sample line of code to run

# csv files are very straightforward in R.  For dta, just replace csv with dta in the code below.  Excel formatted files are trickier - see the linked resources above.
# R is very particular about punctuation and capitalization.  If you have an error, check your syntax for commas, backslashes, quotes, etc.
# And remember, help is just one line of code away:
help(read.csv)
?read.csv

# the easy way. do not use if you want to save your output
data = read.csv(file.choose())

# the harder way. retrieve and set your filepath to wherever your files are.  This is known as your 'working directory.'
getwd()
setwd("C:/Users/Kevin/Documents/a_USC/Rworkshop_Nov2018")
data = read.csv("GDP_Growth.csv")   # Source: OECD Stat Extract 

# to look at your data, just type its name.
data

# but usually all you want are the column names and the number of rows.  Or a summary.
colnames(data)  # this shows column names 
nrow(data)
head(data)
summary(data)
View(data)   # only available in RStudio

# isolate a column using a dollar sign and look at some descriptive statistics another way
area = data$area96
length(area)
mean(area)
median(area)
sd(area)   # standard deviatoin 
min(area)
max(area)
summary(area)

# you can even do a basic linear regression to model which variables determine GDP growth
model = lm(data$GDPGrow ~ data$GDP96 + data$popDens96 + data$unemp99 + data$pctCollege + data$CO2_2005) 
summary(model)

# if we plan to use GDP96 a lot, make it easier by assigning a name
gdp = data$GDP96

# a very basic histogram, then a fancier one.  
hist(gdp)
hist(gdp, breaks = 20, xlim = c(min(gdp),60000), xlab="1996 GDP/capita", ylab = "Number of Regions", col = "blue", border = "red", main = "Regional Output")

# save by right-clicking, or using this cool PDF-making feature
pdf(file="figure1.pdf")
hist(gdp, breaks = 20, xlim = c(min(gdp),60000), xlab="1996 GDP/capita", ylab = "Number of Regions", col = "blue", border = "red", main = "Regional Output")
dev.off()


############## (1b) DATA MANIPULATION ########################

# what if we are interested in looking only at regions with a GDP greater than or equal to $25,000?  This is REALLY handy if your data is too big for Excel (which can only handle 1,048,576 rows) and you need to extract a piece.  
high_gdp = subset(data, GDP96 >= 25000)

# Check the number of rows that remain in your new set.
nrow(high_gdp)

# What if you are only interested in keeping some variables?  Use colnames to see what's contained in the dataset, extract the ones you want, then make a new data item.  Summarize if you want to look at your creation.
colnames(high_gdp)
myvars = c("GDP96", "popDens96", "area96")
high_gdp2 = high_gdp[myvars]
summary(high_gdp2)

# Slicing and dicing a vector of data
gdp = data$GDP96   # first way: use a dollar sign to extract a vector from a data frame
gdp = data[,2]  # another way: Using brackets, extract by [row, column].
gdp[1:10]   # A colon can be used to select a range.
gdpten = gdp[1:10]  # The selection can be assigned a new name
gdpten[gdpten<2000] = NA   # what if you wanted to flag any GDP below a certain threshold as NA or missing. Often this comes up if your data include a missing value code.
gdpten2 = gdpten[!is.na(gdpten)]  # exclude NAs in a vector
length(gdpten2)
sort(gdpten2)
sort(gdpten2, decreasing=TRUE)

# Slicing and dicing a data frame (matrix)
data2 = data[order(data$CO2_2005),] # Sort the entire data frame by the values in one column.  Try adding a negative sign before 'data' here.
head(data2)
data3 = data[-nrow(data),] # remove the last row and save to a new data frame
data4 = data[c(1:4)] # keep only the first four columns
data4$popDens96 = NULL   # remove this variable from the data frame 

# Perhaps you've noticed this data didn't come with an index.  So, we can make a vector from 1 to n
index = seq(1:length(data$GDP96))

# Now we can bind this to our data.  In fact, you can bind any other column like this.
# Note I've named the new data frame high_gdp3. You could name it high_gdp2, and it would overwrite the old, without-an-index version of high_gdp2
data5 = cbind(index, data)

# Save this new dataset for use later... or elsewhere
write.csv(data5, "gdp_subset.csv")



########################################################
######## (2) WORKING WITH MULTIPLE DATA SOURCES ########
########################################################

# Some useful tips 
ls()			# list current objects
rm(object)		# remove object (put the name where "object" is now)
rm(list=ls())   # remove all objects from memory
options(scipen=999)   # avoids printing things in scientific notation
options(stringsAsFactors=FALSE)   # avoids reading in factors when reading in files


######## (2a) MERGE DATASETS FROM DIFFERENT SOURCES AND ANALYZE THEM ####### 
# What predicts domestic migration in California counties? Housing cost, or job growth? 

# (1) Download county migration data. Source: CA Dep't of Finance, http://www.dof.ca.gov/Forecasting/Demographics/Estimates/E-6/ 
dof = read.csv("ca_dof_popchg_e6.csv")

# (2) Dowload median home values from Census' factfinder (factfinder2.census.gov), Advanced Search.
# All California counties, 2017 ACS 5-year sample 
# B25077: Median home value
# Modify table to transpose (so the counties are rows), and download to "Use" the data as a .csv
# Look at it in Excel. Note the 2-part header. Then read into R, omitting the first 2 lines.  
acs = read.csv(text=readLines("./ACS_17_1YR_B25077_with_ann.csv")[(-2)])

# (3) Job info is available from the bureau of labor statistics (https://data.bls.gov/cew/apps/data_views/data_views.htm#tab=Tables)
# For simplicity, I have provided it: 
bls = read.csv("blsemp_ca.csv")

# Join other data to 'dof' data using "match" (will work regardless of whether data are sorted, similar to a v-lookup in Excel)
# creates a new column in 'dof' called 'hoval16,' which is equivalent to the column HD01_VD01 in the acs data
# rows are matched using the unique id in 'dof' (id) and 'acs' (GEO.id2)
dof$hoval16 = acs$HD01_VD01[match(dof$id, acs$GEO.id2)]
dof$emp17 = bls$emp_jun17[match(dof$id, bls$GEOID)]   # same thing for bls data 
dof$emp15 = bls$emp_jun15[match(dof$id, bls$GEOID)]   # one time for each row (bonus: can you turn this into a loop???)

# If we want a variable for percentage population growth or employment growth, we can easily make one:
dof$pctempgr = (dof$emp17 - dof$emp15)/dof$emp15
dof$pctpopgr = (dof$pop17 - dof$pop10)/dof$pop10


######### (2b) ANALYZE THE DATA USING PLOTS AND REGRESSIONS #######

# Compare the size of a county with its median home value.
# A good first step is to always look at how your data are distributed with a histogram
hist(dof$hoval16)
hist(dof$pop17, breaks=20)  # slightly fancier, you can add 'breaks'
boxplot(dof$pctpopgr, dof$pctempgr)   # can also do a boxplot.  This version has two variables side-by-side.
hist(log(dof$pop17))   # these data clearly have right-skew. Transofm this by taking the log.

# Make a scatterplot
x = dof$hoval16
y = log(dof$pop17)
plot(x,y)

# Make a fancier scatterplot with a best-fit line
model = lm(y ~ x)  # stores the results of a simple linear regression in the object 'model'
plot(x, y, ylab="Log of Population, 2017", xlab="Median Home Value, 2016", main="California Counties")   # see help(par) to view all plot options
abline(model, col="blue", lwd=2)   # add a line of 'model' to the existing plot, make it blue with a line width (lwd) of 2
# Run these 5 lines to get the coefficient estimate, significance, and model fit onto the plot
rp = vector('expression', 3)
rp[1] = substitute(expression(italic(B)[1] == MYVALUE3), list(MYVALUE3=format(summary(model)$coefficients[2], dig=4)))[2]
rp[2] = substitute(expression(italic(p) == MYVALUE2), list(MYVALUE2=format(summary(model)$coefficients[2,4], dig=3)))[2]
rp[3] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE=format(summary(model)$r.squared, dig=4)))[2]
legend("bottomright", legend=rp, bty='n')  # can put this at any corner of the plot

# Correlation coefficients can also be done:
cor.test(x, y)

#### Regression analysis - what drives domestic migration? 
# the first command (lm) makes a linear model, and stores the output
# the dependent variable goes first, followed by a tilde and the independent variables. 
# for convenience, you can declare the data frame separately and skip the dollar signs...
m = lm(dommig17 ~ pop10 + hoval16 + pctempgr, data=dof)
summary(m)



#########################################################
####### (3) ADDITIONAL FUNCTIONALITY WITH PACKAGES ######
#########################################################

# In addition to the functions included with the base version of R, there are THOUSANDS of packages with added functionality.
# Making a package is fairly easy, so if someone develops a new method/technique/visualization/etc., it's fairly easy to implement in R.
# Visit a CRAN mirror page to see all the packages available: http://cran.stat.ucla.edu/ 

# Install packages using a single line of code (or, Tools -> Install Packages)
install.packages("sf")   # a basic mapping package
install.packages("doBy")  # allows for some kinds of data summarization
install.packages("scatterplot3d")   # makes 3d scatterplots   
install.packages("foreign")   # packages allows you to read/write .dbf files

# Each time you use R, you'll need to activate the package.
# A good practice is to include this at the beginning of any script which uses the package:
library(sf)
library(scatterplot3d)
library(foreign) 
library(doBy)

# Make a 3D scatterplot
scatterplot3d(log(dof$pop10), log(dof$emp15), dof$hoval16, color="red", pch=10)


######## (3a) MANIPULATING AN ESRI SHAPEFILE ########  
# Often times, you need join additional data to a shapefile to map or analyze, but ESRI's join functionality is slow and not code-based
# By reading and writing the .dbf portion of an ESRI shapefile using the 'foreign' package, this task is a breeze
# BE VERY CAREFUL not to add/delete rows, reorder, or accidentally overwrite the .dbf with something else, as this may ruin the shapefile

# Read in the .dbf portion of the California counties shapefile
shape = read.dbf("cb_2017_ca_county_500k.dbf")   
head(shape)

# When matching, you may have to watch out for strings, leading zeroes, or other issues that may cause ID fields not to match
# Here, we can make a new numeric field in the shapefile to match with our dof file
shape$GEOID2 = as.numeric(as.character(shape$GEOID))

# Match a couple of variables from the DOF data, then summarize to see if it worked
shape$pop17 = dof$pop17[match(shape$GEOID2, dof$id)]
shape$hoval16 = dof$hoval16[match(shape$GEOID2, dof$id)]
summary(shape)

# CAREFULLY write new dbf. Then, open ArcMap. You'll see the new field and can map it!
write.dbf(shape, "cb_2017_ca_county_500k.dbf")


####### (3b) SOME SIMPLE MAPPING USING R & SF PACKAGE ###### 

# Read in the county shapefile
# dsn "." indicates the shapefile resides in the current working directory
# no file extension (e.g., .shp) is required for the layer input
counties = read_sf(dsn = ".", layer = "cb_2017_ca_county_500k")
head(counties)

# Draw a chloropleth map of the "ALAND" field in the counties shapefile
plot(counties["hoval16"])

# A more advanced version. See what you can do at https://cran.r-project.org/web/packages/sf/vignettes/sf5.html#geometry_with_attributes:_sf 
plot(counties["pop17"], breaks="jenks", key.pos=2, pal=sf.colors(10), axes=TRUE)

#################################################
############# (4) PLOTS #########################
#################################################

#### (4a) a fun plot you can do ####
x = seq(1:400)
y = sin(x/10)*exp(x*-0.01)
plot(x,y)

###### (4b) Exercise: Land use types in Phoenix, Arizona #####

# here's some plottable data. note that c() stands for "combine" and renders a list
year = c(1915, 1949, 1963)
residential = c(1371, 995, 282)
commercial = c(178, 471, 530)
vacant = c(547, 232, 483)
nuisance = c(282, 680, 645)
phx = data.frame(year, residential, commercial, vacant, nuisance)

# the very basics.  If you don't specify, R will create default dimensions for a plot.
plot(phx$year, phx$residential)

# that's not very useful. Add x and y limits.  Add more lines to the plot using 'lines.'  Change the type to a line with a point, add color, then width and line type.
plot(phx$year, phx$residential, type='o', xlim=c(1910,1970), ylim=c(100,1400), col='red')
lines(phx$year, phx$commercial, type='o', col='yellow')
lines(phx$year, phx$vacant, type='o', col='blue', lwd=3.5)
lines(phx$year, phx$nuisance, type='o', col='forestgreen', lwd=2, lty=1)

# Here's an even more sophisticated version.  Getting close to publication-worthy!
plot(phx$year, phx$residential, type = "o", col = "cyan2", pch=0, lwd=3, lty=1, ylim = c(100,1400), xlim=c(1915,1963), axes = FALSE, xlab="Year", ylab="Parcel Count")
lines(phx$year, phx$commercial, type = "o", col = "darkgoldenrod1", pch=1, lwd=3, lty=2)
lines(phx$year, phx$vacant, type = "o", col = "darkolivegreen1", pch=2, lwd=3, lty=3)
lines(phx$year, phx$nuisance, type = "o", col = "maroon", pch=3, lwd=3, lty=4)
axis(1, at=c(1915, 1949, 1963), lab=c(1915, 1949, 1963))
axis(2)
box()
abline(h=c(200,400,600,800,1000,1200), lwd=0.5, lty=3, col="gray50")
title(main="Phoenix Land Use Categories", font.main=4)
legend("topright", inset=0.045, c("R", "C", "V", "N"), pch=0:3, lty=1:4, lwd=3, cex=0.7)
title(xlab="Year")
title(ylab="Parcel Count")

######## (4c) Fun with barplots ########

# A handy trick is to be able to summarize data at a different level, or spatial scale.  
# Let's summarize the dof data by the categorical column 'mpo' and put the data you want summarized before it, separated by plus signs
# First, see how many unique values are in the 'mpo' field. This will indicate the number of MPOs in California's 58 counties:
length(unique(unlist(dof$mpo)))
# This relies on the doBy package.
dofsum = summaryBy(pop17 + emp17 ~ mpo, data=dof, FUN=sum)   # takes the sum total of emp17 and pop17 by MPO
dofsum
dofmean = summaryBy(hoval16 + immig17 ~ mpo, data=dof, FUN=mean)   # mean across counties in each MPO. Ensure this makes sense!
dofmean

# Make a single barplot by MPO
barplot(dofsum$pop17.sum, names.arg=dofsum$mpo)

# Make the barplot look better
barplot(dofsum$pop17.sum, names.arg=dofmean$mpo, axes=F, col="dodgerblue", main="2017 Population by MPO (millions)")
axis(2, at=seq(0,20000000,5000000), lab=seq(0,20,5))
abline(h=seq(0,20000000,5000000), lty=3, col="darkgrey")
box()
