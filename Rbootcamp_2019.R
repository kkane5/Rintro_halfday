###############################
###### Introduction to R ######
###############################
# by Kevin Kane, PhD

# Some R-esources!
# R website (open source): www.r-project.org 
# Rstudio website (freeware): www.rstudio.com 
# Quick-R: www.statmethods.net 

##############################################
######### (1) BASIC DATA SKILLS IN R #########
##############################################

# Ways to run code:
# (1) Type, or copy/paste into the Console window's prompt (>)
# (2) Place your cursor on a line, or highlight 1 or more lines, and press "Run" at the top-right of the script window
# (3) Place your cursor on a line, or highlight 1 or more lines, and hit Ctrl+R or Ctrl+Enter, depending on your system
# Note that the pound/hashtag converts a line of code to a comment

print("hello world")   # sample line of code to run

# csv files are very straightforward in R.  Excel files are tricky and in my opinion not worth the time - just save as .csv from Excel.
# R is very particular about punctuation and capitalization.  If you have an error, check your syntax for commas, backslashes, quotes, etc.
# And remember, help is just one line of code away:
help(read.csv)
?read.csv

# the easy way. do not use if you want to save your output
data = read.csv(file.choose())

# the harder way. retrieve and set your filepath to wherever your files are.  This is known as your 'working directory.'
getwd
setwd("C://Users//Kevin//Dropbox/usc_emup/R_tutorials")
setwd("C://Users//kane//Dropbox//usc_emup//R_tutorials")
data = read.csv("oecd.csv")

# to look at your data, just type its name.
# when you read in data, it becomes a DATA FRAME in R.
data

# but usually all you want are the column names and the number of rows.  Or a summary of the DATA FRAME.
colnames(data)  
nrow(data)
head(data)
summary(data)
View(data)   

# You can isolate a column of a DATA FRAME so that it's just a VECTOR.  
# This is easiest using a dollar sign. Then, you can look at some descriptive statistics of the VECTOR.
pop = data$totpop17
pop
length(pop)
mean(pop)
median(pop)
min(pop)
max(pop)
summary(pop)

# what if we are interested in looking only at regions with unemployment greater than or equal to 10% ? 
# You can create a new DATA FRAME using SUBSET
# Note that == is "equals exactly" and != is "not equals"
hi_unemp = subset(data, unemp17 >= 10)
nrow(hi_unemp)
nrow(hi_unemp)/nrow(data)

# You can make a very basic histogram, then a fancier one by passing more arguments to the command "hist"  
hist(data$unemp17)
hist(data$unemp17, breaks = 20, xlim = c(1,max(data$unemp17, na.rm=T)), xlab="2017 Unemployment Rate", 
     ylab = "Number of Regions", col = "blue", border = "red", main = "Unemployment in World Regions")


#### QUIZ QUESTIONS - CHALLENGE YOURSELF ####
# What is the unemployment rate (rounded to the nearest whole number) in the Mazowia region of Poland, which has a code of PL12 (see http://www.oecd.org/cfe/regional-policy/territorial-grid-2018.pdf for a description of OECD regions)
# Make a histogram of 2016 household income in Italian regions. Italy's country code is "IT" and you can use the subset command with 'exactly equals' . In the plot window, click "export."



############################################
######### (2) DESCRIBING DATA IN R #########
############################################

# Some useful tips 
ls()			# list current objects (data frames, vectors, model outputs, etc.)
rm(object)		# remove an object - the matrix unemp, in this example
rm(list=ls())   # remove all objects from memory
options(scipen=999)   # avoids printing things in scientific notation
options(stringsAsFactors=FALSE)   # avoids reading in factors when reading in files

# Return to the OECD data
setwd("C:/Users/Kevin/Dropbox/usc_emup/R_tutorials")
setwd("C:/Users/kane/Dropbox/usc_emup/R_tutorials")
data = read.csv("oecd.csv")

# Let's make some percentage - and other - variables using mathematical operations
data$pctsenior17 = data$popovr6517/data$totpop17
summary(data$pctsenior17)

# Which has the highest number of unemployed PEOPLE?
data$popunemp17 = data$unemp17/100*data$workagepop17
summary(data$popunemp17)
summary(data$totpop17)
View(data)   # clicking on column headers sorts the data temporarily 

# Compare senior share in 2001 versus 2017 
data$pctsenior01 = data$popovr6501/data$totpop01
data$pctsrincr = data$pctsenior17 - data$pctsenior01
hist(data$pctsrincr)

# Sort in code. Actually, this makes a new data frame which sorts by the desired column.  
data2 = data[order(data$pctsrincr),]   
head(data2)
data2 = data[order(-data$pctsrincr),]  # this repeats the above, but the negative sign makes it in descending order
data2[1:10,]   # data frames can be sliced using brackets. The format is row, column. This isolates rows 1 through 10, and leaves the columns untouched
data2[1:10,1:5]   # this isolates rows 1-10 and columns 1-5 

# Adding a rank field 
srrank = seq(1:length(data2$pctsrincr))   # seq makes a vector from 1 until a specified value
srrank
data2 = cbind(data2, srrank)   # cbind stands for column bind. this pastes srrank to data2 and overwrites data2 with the result.

# Make comparative histograms
par(mfrow=c(1,2))  # par stands for graphical parameter. mfrow specifies the number of rows and columns in the plot window
hist(data$pctsenior01)   # left plot
hist(data$pctsenior17)   # right plot 
par(mfrow=c(1,2))   # doing this again clears the plots 
hist(data$pctsenior01, main="2001", ylab="World Regions", xlab="Percent Senior", col="forestgreen")
hist(data$pctsenior17, main="2017", ylab="World Regions", xlab="Percent Senior", col="dodgerblue")

# Adding a guideline (abline) to each histogram. can be vertical (v) or horizontal (h)
par(mfrow=c(1,2))
hist(data$pctsenior01, main="2001", ylab="World Regions", xlab="Percent Senior", col="forestgreen")
abline(v=mean(data$pctsenior01, na.rm=T), col="black", lty=3, lwd=2)
hist(data$pctsenior17, main="2017", ylab="World Regions", xlab="Percent Senior", col="dodgerblue")
abline(v=mean(data$pctsenior17, na.rm=T), col="black", lty=3, lwd=2)

# Make a side-by-side boxplot to compare  
boxplot(data$pctsenior01, data$pctsenior17)

#### CORRELATION ANALYSIS #### 
# How are vehicle ownership and carbon emissions related?
summary(data)
View(data)
cor.test(data$vehicper100_13, data$co2tot08)   # the entire output of a Pearson's correlation test. 
help(cor.test)

# Perhaps per-capita emissions is a better way to analyze this relationship
cor.test(data$vehicper100_13, data$co2cap08)

# Scatterplot
par(mfrow=c(1,1))   # reset graphics parameters to a single plot
plot(data$vehicper100_13, data$co2cap08)  # plot command will plot 2 vectors against each other 
data3 = subset(data, co2cap08<50)  # make a new data frame to remove high outliers
nrow(data)
nrow(data3)
plot(data3$vehicper100_13, data3$co2cap08)
cor.test(data3$vehicper100_13, data3$co2cap08)
plot(data3$vehicper100_13, data3$co2cap08)

# Scatterplot and best-fit line. Actually requires you to run a simple linear regression (lm command) 
x = data3$vehicper100_13   # assign this vector a shorter, more convenient name. 
y = data3$co2cap08
model = lm(y ~ x)   # runs a 2-variable regression and saves the output to an object I've called 'model'
summary(model)   # regression output 
outcorr = cor.test(data3$vehicper100_13, data3$co2cap08)  # sometimes it's convenient to assign the output of a test to a model so you can manipulate it later 
outcorr
outcorr$estimate

plot(x,y, xlab="Vehicle Ownership per 100 population", ylab = "Carbon emissions per capita", 
     main = "World region vehicle ownership vs. carbon emissions")
abline(model)   # abline can also add model results to a plot 
legend("topright", inset=0.045, legend=c("r =", round(outcorr$estimate,4)))


# A cleaned-up version, starting with the line "plot" above 
plot(x, y, xlab="Vehicle Ownership per 100 population (2013)", ylab = "Carbon emissions per capita", 
     main = "World region vehicle ownership vs. carbon emissions (2008)", pch=20, font.main=4, 
     ylim=c(0,40), xlim=c(0,80))
abline(model, col="blue", lwd=2)
rp = vector('expression', 3)
rp[1] = substitute(expression(italic(B)[1] == MYVALUE3), list(MYVALUE3=format(summary(model)$coefficients[2], dig=4)))[2]
rp[2] = substitute(expression(italic(p) == MYVALUE2), list(MYVALUE2=format(summary(model)$coefficients[2,4], dig=3)))[2]
rp[3] = substitute(expression(italic(r) == MYVALUE), list(MYVALUE=format(sqrt(summary(model)$r.squared), dig=4)))[2]
legend("topright", legend=rp, bty='n')


#### DEALING WITH MISSING VALUES #### 
summary(data)
mean(data$hhinc01)
mean(data$hhinc01, na.rm=T)
data4 = data[!is.na(data$popovr6501),]   # isolate rows wherein popovr6501 is NA, remove them, and save to 'data4'
nrow(data)
nrow(data4)   # notice how many fewer observations are now in data4 

# Other helpful operations on a data frame 
data$pctsrincr = NULL
data5 = data[c(1:4)]


#### QUIZ QUESTIONS - CHALLENGE YOURSELF ####
# (3) How many world regions have a 2017 population density over 200/sqmi? 
# (4) What is the correlation coefficient between income in 2001 and population size in 2001? Round to 2 decimal places





########################################################
######## (3) WORKING WITH MULTIPLE DATA SOURCES ########
########################################################

######## MERGE DATASETS FROM DIFFERENT SOURCES AND ANALYZE THEM ####### 
# What predicts domestic migration in California counties? Housing cost, or job growth? 

# Download county migration data. Source: CA Dep't of Finance, http://www.dof.ca.gov/Forecasting/Demographics/Estimates/E-6/ 
dof = read.csv("ca_dof_popchg_e6.csv")
head(dof)
nrow(dof)

# Download median home values from Census' factfinder (factfinder2.census.gov), Advanced Search.
# All California counties, 2017 ACS 5-year sample 
# B25077: Median home value
# Modify table to transpose (so the counties are rows), and download to "Use" the data as a .csv
# Open first in Excel. Delete line 2 (there are 2 lines' worth of header). Resave and close before reading into R. 
acs = read.csv("ACS_17_5YR_B25077_with_ann.csv")

# Job info is available from the bureau of labor statistics (https://data.bls.gov/cew/apps/data_views/data_views.htm#tab=Tables)
bls = read.csv("blsemp_ca.csv")

# Join other data to 'dof' data using "match" (will work regardless of whether data are sorted, similar to a v-lookup in Excel)
# creates a new column in 'dof' called 'hoval16,' which is equivalent to the column HD01_VD01 in the acs data
# rows are matched using the unique id in 'dof' (id) and 'acs' (GEO.id2)
dof$hoval16 = acs$HD01_VD01[match(dof$id, acs$GEO.id2)]
dof$emp17 = bls$emp_jun17[match(dof$id, bls$GEOID)]   # same thing for bls data 
dof$emp15 = bls$emp_jun15[match(dof$id, bls$GEOID)]   # one time for each row 

# Make a couple of percentage growth variables using mathematical operations:
dof$pctempgr = (dof$emp17 - dof$emp15)/dof$emp15
dof$pctpopgr = (dof$pop17 - dof$pop10)/dof$pop10


######### ANALYZE THE DATA USING PLOTS AND A REGRESSION MODEL #######

# Compare the size of a county with its median home value. 
# Histograms help you understand how the data are distributed!
hist(dof$pop17, breaks=20) 
# Taking the logarithm of non-zero data can get rid of skew issues 
# Right skew is common in social science data - high values pull up the means. Using a log transformation is particularly useful when comparing data, e.g. understanding the correlation between two variables. 
hist(log(dof$pop17))
boxplot(dof$pctpopgr, dof$pctempgr)   # side-by-side boxplot.

# Make a scatterplot and check the correlation
x = dof$hoval16
y = log(dof$pop17)
plot(x,y)
cor.test(x, y)

# Make a fancier scatterplot with a best-fit line
model = lm(y ~ x)  # stores the results of a simple linear regression in the object 'model'
plot(x, y, ylab="Log of Population, 2017", xlab="Median Home Value, 2016", main="California Counties")   # see help(par) to view all plot options
abline(model, col="dodgerblue", lwd=3)   
rp = vector('expression', 3)
rp[1] = substitute(expression(italic(B)[1] == MYVALUE3), list(MYVALUE3=format(summary(model)$coefficients[2], dig=4)))[2]
rp[2] = substitute(expression(italic(p) == MYVALUE2), list(MYVALUE2=format(summary(model)$coefficients[2,4], dig=3)))[2]
rp[3] = substitute(expression(italic(r) == MYVALUE), list(MYVALUE=round(sqrt(summary(model)$r.squared), 4)))[2]
legend("bottomright", legend=rp, bty='n')  # can put this at any corner of the plot




##########################################################
#### (4) SUMMARIZING DATA AT DIFFERENT SPATIAL SCALES ####
##########################################################

####### ADDITIONAL FUNCTIONALITY WITH PACKAGES ######

# In addition to the functions included with the base version of R, there are THOUSANDS of packages with added functionality.
# Making a package is fairly easy, so if someone develops a new method/technique/visualization/etc., it's fairly easy to implement in R.
# Visit a CRAN mirror page to see all the packages available: http://cran.stat.ucla.edu/ 

# Install packages using a single line of code (or, Tools -> Install Packages)
install.packages("sf")   # a basic mapping package
install.packages("doBy")  # allows for some kinds of data summarization
install.packages("foreign")   # packages allows you to read/write .dbf files

# Each time you use R, you'll need to activate the package.
# A good practice is to include this at the beginning of any script which uses the package:
library(sf)
library(foreign) 
library(doBy)


#### SUMMARY BY ####  
# A handy trick is to be able to summarize data at a different level, or spatial scale.  

# First, see how many unique values are in the 'mpo' field.
unique(unlist(dof$mpo))
length(unique(unlist(dof$mpo)))

# Return to the OECD data for a moment to see how many regions are in each country
setwd("C:/Users/Kevin/Dropbox/usc_emup/R_tutorials")
setwd("C:/Users/kane/Dropbox/usc_emup/R_tutorials")
data = read.csv("oecd.csv")

# Make a sorted table and a barplot of the number of regions per country 
table(data$country)
sort(table(data$country), decreasing=TRUE)  # the sort command, with the argument decreasing=TRUE will sort the table in decreasing order.
barplot(sort(table(data$country), decreasing=TRUE), cex.names=0.8,
        main="Number of regions per country")   # a barplot can be made (can also add colors, etc.). cex.names makes the text smaller so all the country names fit.


# Use the summaryBy command (part of the doBy package) to summarize by larger spatial unit 
# a new data frame will be created which summarizes the variables to the left of the tilde at the unit defined to the right of the tilde
View(dof)
dofsum = summaryBy(pop17 + emp17 ~ mpo, data=dof, FUN=sum)   # takes the sum total of emp17 and pop17 by MPO
dofsum
dofmean = summaryBy(hoval16 + immig17 ~ mpo, data=dof, FUN=mean)   # mean across counties in each MPO. Ensure this makes sense!
dofmean

# Write the table output as a .csv to your working directory
write.csv(dofsum, "dofsumtable.csv")

# Make a single barplot by MPO
barplot(dofsum$pop17.sum, names.arg=dofsum$mpo)

# Make the barplot look better
barplot(dofsum$pop17.sum, names.arg=dofmean$mpo, axes=F, col="dodgerblue", main="2017 Population by MPO (millions)")
axis(2, at=seq(0,20000000,5000000), lab=seq(0,20,5))
abline(h=seq(0,20000000,5000000), lty=3, col="darkgrey")
box()

# Summarize some of the OECD data. Note that rate or percentage data doesn't collapse as well. 
out = summaryBy(areasqkm17 + totpop17 + co2tot08 ~ country, data=data, FUN=sum)
out2 = out[order(-out$co2tot08.sum),]  # sort by CO2 
out2
cor.test(out2$co2tot08.sum, out2$totpop17.sum)   # are larger countries higher emitters? 


####  REGRESSION ANALYSIS - THE VERY BASICS ######
# the first command (lm) makes a linear model, and stores the output
# the dependent variable goes first, followed by a tilde and the independent variables. 
# for convenience, you can declare the data frame separately and skip the dollar signs...
m = lm(dommig17 ~ pop10 + hoval16 + pctempgr, data=dof)
summary(m)




###########################################
####### (5) MAPPING MADE SIMPLE(ISH) ###### 
###########################################
setwd("C:/Users/kane/Dropbox/usc_emup/R_tutorials")
library(sf)

#### SIMPLE MAP USING THE SF PACKAGE #### 
# Read in a county shapefile
# dsn "." indicates the shapefile resides in the current working directory
# no file extension (e.g., .shp) is required for the layer input
counties = read_sf(dsn = ".", layer = "cb_2017_ca_county_500k")
head(counties)
cnty = data.frame(counties)   # if you just want to extract the attribute/data table 

# Draw a chloropleth map of the "ALAND" field in the counties shapefile
plot(counties["ALAND"])

# Plot a variable which may actually be meaningful
plot(counties["pop17"])

# A more advanced version
# Review plot options at https://cran.r-project.org/web/packages/sf/vignettes/sf5.html#geometry_with_attributes:_sf 
plot(counties["pop17"], breaks="jenks", key.pos=2, pal=sf.colors(10), 
     main="2017 Population in California counties", axes=TRUE)



#############################################
######## (6) MANIPULATING SHAPEFILES ########  
#############################################
# GIS data often lack efficient or comprehensive data joining capability 
# By reading and writing the .dbf portion of an ESRI shapefile using the 'foreign' package, this task is a breeze
# BE VERY CAREFUL not to add/delete rows, reorder, or accidentally overwrite the .dbf with something else, as this may ruin the shapefile

# Read in the .dbf portion of the California counties shapefile
shape = read.dbf("cb_2017_ca_county_500k.dbf")   
head(shape)

# When matching, you may have to watch out for strings, leading zeroes, or other issues that may cause ID fields not to match
# Here, we can make a new numeric field in the shapefile to match with our dof file
# The nested as.numeric(as.character()) commands pretty reliably coerce things into numeric format
shape$GEOID2 = as.numeric(as.character(shape$GEOID))

# Match a couple of variables from the DOF data, then summarize to see if it worked
shape$pop17 = dof$pop17[match(shape$GEOID2, dof$id)]
shape$hoval16 = dof$hoval16[match(shape$GEOID2, dof$id)]
summary(shape)

# CAREFULLY write new dbf. Then you'll be able to open in a GIS software and map the new field! 
write.dbf(shape, "cb_2017_ca_county_500k.dbf")



######################################
###### (7) Using the Census API ######
######################################

#### Install TidyCensus & tigris, and load other packages we'll need ####
install.packages("tidycensus")
install.packages("tigris")
library(tidycensus)
library(tigris)
library(sf)
library(doBy)

#### Request a key from the Census to access their data #### 
# You'll need to request a Census Key to access their API. 
# This takes only one minute. Then, enter your key below (mine is in there now)
# http://api.census.gov/data/key_signup.html
census_api_key("5de05eeaa3869c3180eed1726a07e0f27e6c9e20", install=TRUE, overwrite=TRUE)

# Extract your first variable from the decennial census using the function get_decennial
medrent90 = get_decennial(geography = "state", variables = "H043A001", year = 1990)

# Make a barplot
barplot(medrent90$value, names.arg=medrent90$NAME, cex.names=0.7, las=2)

# Make a better barplot
medrent90b = medrent90[order(-medrent90$value),]		# sort descending by making a new data frame
barplot(medrent90b$value, names.arg=medrent90b$NAME, cex.names=0.7, las=2,
        main="State median rent, 1990", col="dodgerblue")   #cex.names makes labels smaller, las=2 rotates them.
abline(h=seq(0,700,100), col="darkgrey", lty=3)
box()


# Grab a blank state-level shapefile at https://www2.census.gov/geo/tiger/TIGER2018/STATE/tl_2018_us_state.zip 
# Unzip and place into your working directory
# Restrict to lower 48 states to keep the map looking good (sorry Alaska, Hawaii, and territories!)
states = read_sf(dsn = ".", layer = "tl_2018_us_state")
states = subset(states, GEOID!="78" & GEOID!="69" & GEOID!="66" & GEOID!="02" & GEOID!="60" & GEOID!="72" & GEOID!="15")

# Plot a basic map of the states 
plot(st_geometry(states), col = "grey")

# Match rent data to the shapefile and plot
states$medrent90 = medrent90$value[match(states$GEOID, medrent90$GEOID)]
plot(states["medrent90"])

# pick a good way to display your data
# with 'breaks', make sure you have zero and a high value too. 
hist(states$medrent90)
plot(states["medrent90"], breaks=c(0, 350, 500, 600, 1000))

# You can explicitly make quintiles, or use Jenks optimization 
states$rent90q = cut(states$medrent90, 5)
plot(states["rent90q"])
plot(states["medrent90"], breaks="jenks", main="1990 Median Rent by state")

# You can also get a county-level variable for just counties in California.
# On your own, you might try and bind this to your CA county data from the previous lab, and map it 
medrentcounties90 <- get_decennial(geography = "county", state = "CA", variables = "H043A001", year = 1990)



#### FINDING GOOD CENSUS VARIABLES TO USE #### 
# Focusing on the most recent ACS 5-year estimates (2017 at the time of this writing) is a good way to go. The 2010 Decennial is getting fairly old.
# https://api.census.gov/data/2017/acs/acs5/profile/variables.html 
# All the Census API available datasets are at https://api.census.gov/data.html 
# You can also extract a list of variables using the code here.
acsvars = load_variables(2017, "acs5", cache = TRUE)
head(acsvars)
nrow(acsvars)   # tells us how many variables
length(unique(unlist(acsvars$concept)))   # tells us how many unique "concepts" there are
write.csv(acsvars, "acsvars.csv")   # export to .csv so you explore easily in Excel 


#### ASSEMBLE A GOOD SET OF TRACT-LEVEL VARIABLES FOR A COUNTY #### 
# Get a single variable (Median home value)
vent = get_acs(geography="tract", state="CA", county="Ventura", variables="B25007_001", 
               year=2017, geometry=TRUE)
vent$medhoval = vent$estimate  # create a renamed to avoid future confusion
vent$estimate = NULL   # get rid of the old one 

# Add additional variables with a single extraction by putting them into a list. 
# Total population: B00101_001
# Median age: B01002_001
# Median household income: B19013_001
# Median age of housing: B25035_001

varlist = c("B01001_001", "B01002_001", "B19013_001", "B25035_001")
out =  get_acs(geography="tract", state="CA", county="Ventura", variables=varlist, year=2017)

# Use a match command to bind each new variable to the original data frame (vent). name the variables something intuitive.
# Since the data are stored long, each match operation requires a subset operation first.
venta = subset(out, variable=="B01001_001")
vent$totpop = venta$estimate[match(vent$GEOID, venta$GEOID)]

ventb = subset(out, variable=="B01002_001")
vent$medage = ventb$estimate[match(vent$GEOID, ventb$GEOID)]

ventc = subset(out, variable=="B19013_001")
vent$medhhinc = ventc$estimate[match(vent$GEOID, ventc$GEOID)]

ventd = subset(out, variable=="B25035_001")
vent$medhomeage = ventd$estimate[match(vent$GEOID, ventd$GEOID)]

# Make some more variables using a loop. 
# Total commuters: B08101_001
# Single-occupant vehicle commuters: B08101_009
# Work at home: B08101_001

varlist2 = c("B08101_001", "B08101_009", "B08101_049")
varnames2 = c("totcommute", "SOVcommute", "workathome")
out2 =  get_acs(geography="tract", state="CA", county="Ventura", 
                variables=varlist2, year=2017)

# Try looping through so as to avoid having to repeat for each variable.
# The loop command iterates through the sequence of numbers from 1 to 3 using "i" as the wildcard
# The first line subsets the extraction to just the ith item in varlist2
# The second line adds a new column to vent (one column higher than its current number of columns...). This structure avoids having to declare a specific variable name and allows you to loop through variables.
# The third line renames the last column to the ith element in varnames2, which you made above.  
for(i in 1:3){
  join = subset(out2, variable==varlist2[i])
  vent[,ncol(vent)+1] = join$estimate[match(vent$GEOID, join$GEOID)]
  colnames(vent)[ncol(vent)] = varnames2[i]  
}

# Many census variables have a restricted universe, e.g. those over age 16 are isolated to gauge labor force participation, or the variable just otherwise needs to be a percentage.
vent$sovshare = vent$SOVcommute/vent$totcommute
vent$wahshare = vent$workathome/vent$totcommute


#### AUTOMATICALLY EXTRACT SHAPEFILES FROM CENSUS #### 
# We will use the tigris package to extract Census shapefiles directly. 
# See details at http://rstudio-pubs-static.s3.amazonaws.com/90665_de25062951e540e7b732f21de53001f0.html 

# The function is the name of the geography you want - tracts, in this case. 
# Let's grab all of southern California, then just Ventura county to continue from above 
scagtr = tracts(state = 'CA', county = c('Imperial', 'Los Angeles', 'Orange', 'Riverside', 'San Bernardino', 'Ventura'))
plot(scagtr)
venttr = tracts(state = 'CA', county = 'Ventura')
plot(venttr)


##!?!?!?## RESET FOR A SECOND ##!?!?!?##
# We can also read in the data ... as a shapefile! 
# Previously we'd run the below.  BUT, if you add geometry=TRUE at the end, the upload will include spatial reference data.
ventnew = get_acs(geography="tract", state="CA", county="Ventura", variables="B25007_001", 
                  year=2017, geometry=TRUE)
plot(ventnew['estimate'])

# Now, go back to the section titled 'assembling a good set of tract-level...' and 
# add 'geometry=TRUE' to the end of the get_acs command.
# Re-run everything through the end of that section

# Try some plotting
plot(vent['sovshare'])

# If the island is messing up your view, remove that tract, and plot something else.
vent_noisland = subset(vent, GEOID!="06111980000")
plot(vent_noisland['medhomeage'], main="Median home ages across Ventura County, CA")


# Finally, we can export this to a shapefile so we can use it elsewhere, too.
st_write(vent, "vent_merge.shp")








