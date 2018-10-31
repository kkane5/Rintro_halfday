# This file was created by Kevin Kane for an introductory lecture on R to graduate students at USC The data employed here are not my own. For attribution information contact me via email at kevin@kevinkane.org

# My favorite R cheat-sheets:
http://www.gardenersown.co.uk/Education/Lectures/R/
Quick-R: www.statmethods.net 


############### DATA BASICS #################
# Ways to run code:
# (1) Type, or copy/paste into the Console window's prompt (>)
# (2) Place your cursor on a line, or highlight 1 or more lines, and press "Run" at the top-right of the script window
# (3) Place your cursor on a line, or highlight 1 or more lines, and hit Ctrl+R
# Note that the pound/hashtag converts a line of code to a comment

# csv files are very straightforward in R.  For dta, just replace csv with dta in the code below.  Excel formatted files are trickier.
# R is very particular about punctuation and capitalization.  If you have an error, check your syntax for commas, backslashes, quotes, etc.
# And remember, help is just one line of code away:
help(read.csv)

# the easy way. do not use if you want to save your output
data = read.csv(file.choose())

# the harder way. retrieve and set your filepath to wherever your files are.  This is known as your 'working directory.'
getwd()
setwd("C:/Users/Kevin/Documents/a_USC/Rworkshop_Nov2018")
data = read.csv("GDP_Growth.csv")
setwd("C:/Users/Kevin/Documents/a_USC")   # An alternative way
data = read.csv("./Rworkshop_Nov2018/GDP_Growth.csv")

# to look at your data, just type its name.
data

# but usually all you want are the column names and the number of rows.  Or a summary.
colnames(data)
nrow(data)
head(data)
summary(data)
View(data)   # only available in RStudio

# isolate a column using a dollar sign and look at some descriptive statistics another way
area = data$area96
mean(area)
median(area)
sd(area)
min(area)
max(area)
summary(area)

# you can even do a basic linear regression to model which variables determine GDP growth
model = lm(data$GDPGrow~ (data$GDP96 + data$popDens96 + data$unemp99 + data$pctCollege + data$CO2_2005)) 
summary(model)

# extract the results of your regression model:
beta1=summary(model)$coefficients["data$GDP96","Estimate"]
tvalue2=summary(model)$coefficients["data$popDens96","t value"]
tprob2=summary(model)$coefficients["data$popDens96","Pr(>|t|)"]

# if we plan to use GDP96 a lot, make it easier by assigning a name
gdp = data$GDP96

# a very basic histogram, then a fancier one.  
hist(gdp)
hist(gdp, breaks = 20, xlim = c(min(gdp),60000), xlab="1996 GDP", ylab = "Number of Regions", col = "blue", border = "red", main = "Regional Growth")

# save by right-clicking, or using this cool PDF-making feature
pdf(file="figure1.pdf")
hist(gdp, breaks = 20, xlim = c(min(data$GDP96),60000), xlab="1996 GDP", ylab = "Number of Regions", col = "blue", border = "red", main = "Regional Growth")
dev.off()


############## DATA MANIPULATION ########################

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
gdp = data$GDP96
gdp = data[,2]  # another way to extract a vector from a data frame. Using brackets, extract by [row, column].
gdp[1:10]   # A colon can be used to select a range.
gdpten = gdp[1:10]  # The selection can be assigned a new name
gdpten[gdpten<2000] = NA   # what if you wanted to flag any GDP below a certain threshold as NA or missing. Often this comes up if your data include a missing value code.
gdpten[!is.na(gdpten)]  # exclude NAs in a vector
sort(gdpten)

# Slicing and dicing a data frame (matrix)
data2 = data[order(data$CO2_2005),] # Sort the entire data frame by the values in one column
head(data2)
d3 = d[-nrow(d),] # remove the last row and save to a new data frame
data = data[c(1:4)] # keep only the first four columns
data$variable = NULL   # remove this variable from the data frame 

# Perhaps you've noticed this data didn't come with an index.  So, we can make a vector from 1 to n
index = seq(1:length(high_gdp2$GDP96))

# Now we can bind this to our data.  In fact, you can bind any other column like this.
# Note I've named the new data frame high_gdp3. You could name it high_gdp2, and it would overwrite the old, without-an-index version of high_gdp2
high_gdp3 = cbind(index, high_gdp2)

# Save this new dataset for use later... or elsewhere
write.csv(high_gdp3, "gdp_subset.csv")



####################################
######## MORE FUN WITH DATA (starting w/line 156 in refcard) ########
####################################

# Some useful tips 
ls()			# list current objects
rm(object)		# remove object
rm(list=ls())   # remove all objects from memory

options(scipen=999)   # avoids printing things in scientific notation
options(stringsAsFactors=FALSE)   # avoids reading in factors when reading in files


######## AMERICAN FACT FINDER DATA EXERCISE ####### 

# Dowload THIS from Census. (2017 ACS stats by county ? )

# Look at it in Excel. Note the 2-part header. Then read into R, omitting the first 2 lines.  
d = read.csv(text=readLines("./ACS_14_5YR_B01001_with_ann.csv")[(-2)])

# Collapse this to counties using summaryBy

# Join other county data using "match" 



####### LINEAR REGRESSIONS AND ASSOCIATED PLOTS ##### 




####################################
####### WORKING WITH PACKAGES ######
####################################

install.packages("car")
install.packages("foreign")



######## MANIPULATING AN ESRI SHAPEFILE ########  

# use shapefile data from above 







############# PLOTS #########################

# a fun plot you can do
x = seq(1:400)
y = sin(x/10)*exp(x*-0.01)
plot(x,y)

# here's some new plottable data
phx = read.csv("landuse.csv")

# the very basics.  If you don't specify, R will create default dimensions for a plot.
plot(phx$Year, phx$R)
xlim=c(1910,1970), ylim=c(100,1400)

# that's not very useful. Add limits.  Add lines.  Change the type to a line with a point, add color, then width and line type.
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
legend("topright", inset=0.045, c("R", "C", "V", "N"), pch=0:3, lty=1:4, lwd=3)
title(xlab="Year")
title(ylab="Parcel Count")


######## VERI BRIEF FORAY INTO R SHINY #### 
install.packages("shiny")



