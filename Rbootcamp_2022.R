###############################################
###### Introduction to R - Bootcamp 2022 ######
###############################################
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

# the easy way.  the tradeoff is that it's harder to save your output. 
data = read.csv(file.choose())

# the harder way. retrieve and set your filepath to wherever your files are.  This is known as your 'working directory.'
getwd()
setwd("C:/Users/kane/Dropbox/usc_emup/R_bootcamp")
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
pop = data$totpop20
pop
length(pop) 
mean(pop)   
median(pop)
min(pop)
max(pop)
summary(pop)

# what if we are interested in looking only at regions with unemployment greater than or equal to 10% ? 
# You can create a new DATA FRAME using SUBSET
# Note that == is "equals exactly" and != is "not equals" and non-numeric values such as place names need to be placed in quotes when subsetting
hi_unemp = subset(data, unemp20 >= 0.1)
nrow(hi_unemp)
nrow(hi_unemp)/nrow(data)

# You can make a very basic histogram, then a fancier one by passing more arguments to the command "hist"  
hist(data$unemp20)
hist(data$unemp20, breaks = 20, xlim = c(0,max(data$unemp20, na.rm=T)), xlab="2020 Unemployment Rate", 
     ylab = "Number of Regions", col = "blue", border = "red", main = "Unemployment in World Regions")


#### QUIZ QUESTIONS - CHALLENGE YOURSELF ####
# What is the unemployment rate in the Central Bohemia Region (code CZ02) of the Czech Republic (code CZ)? 
# Make a histogram of 2019 income levels in Italian regions. Italy's country code is "IT" and you can use the subset command with 'exactly equals' . In the plot window, click "export."



############################################
######### (2) DESCRIBING DATA IN R #########
############################################

# Some useful tips 
ls()			# list current objects (data frames, vectors, model outputs, etc.)
rm(hi_unemp)		# remove an object - the matrix unemp, in this example
rm(list=ls())   # remove all objects from memory
options(scipen=999)   # avoids printing things in scientific notation
options(stringsAsFactors=FALSE)   # avoids reading in factors when reading in files

# Return to the OECD data
setwd("C:/Users/kane/Dropbox/usc_emup/R_bootcamp")
data = read.csv("oecd.csv")

# Let's make some percentage - and other - variables using mathematical operations
data$pctsenior20 = data$pop65plus20/data$totpop20 
summary(data$pctsenior20)
View(data)   # clicking on column headers sorts the data temporarily 

# Compare senior share in 2005 versus 2020 
data$pctsenior05 = data$pop65plus05/data$totpop05 
data$pctsrincr = data$pctsenior20 - data$pctsenior05
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
hist(data$pctsenior05)   # left plot
hist(data$pctsenior20)   # right plot 
par(mfrow=c(1,2))   # doing this again clears the plots 
hist(data$pctsenior05, breaks=seq(0,0.35,0.05), main="2001", ylab="World Regions", xlab="Percent Senior", col="forestgreen")
hist(data$pctsenior20, breaks=seq(0,0.35,0.05), main="2020", ylab="World Regions", xlab="Percent Senior", col="dodgerblue")

# Adding a guideline (abline) to each histogram. can be vertical (v) or horizontal (h)
par(mfrow=c(1,2))
hist(data$pctsenior05, breaks=seq(0,0.35,0.05), main="2001", ylab="World Regions", xlab="Percent Senior", col="forestgreen")
abline(v=mean(data$pctsenior05, na.rm=T), col="black", lty=3, lwd=2)
hist(data$pctsenior20, breaks=seq(0,0.35,0.05), main="2020", ylab="World Regions", xlab="Percent Senior", col="dodgerblue")
abline(v=mean(data$pctsenior20, na.rm=T), col="black", lty=3, lwd=2)

# Make a side-by-side boxplot to compare  
boxplot(data$pctsenior05, data$pctsenior20)

#### CORRELATION ANALYSIS #### 
# How are income and PM2.5 air pollution levels related in world regions?
summary(data)
View(data)
cor.test(data$inc19, data$pm2pt5_19)   # the entire output of a Pearson's correlation test. 
help(cor.test)

# Scatterplot
par(mfrow=c(1,1))   # reset graphics parameters to a single plot
plot(data$inc19, data$pm2pt5_19)  # plot command will plot 2 vectors against each other 

# Scatterplot and best-fit line. Actually requires you to run a simple linear regression (lm command) 
x = data$inc19   # assign this vector a shorter, more convenient name. 
y = data$pm2pt5_19
model = lm(y ~ x)   # runs a 2-variable regression and saves the output to an object I've called 'model'
summary(model)   # regression output 
outcorr = cor.test(data$inc19, data$pm2pt5_19)  # sometimes it's convenient to assign the output of a test to a model so you can manipulate it later 
outcorr
outcorr$estimate

plot(x,y, xlab="2019 Personal Income ($USD)", ylab = "Average air pollution level (PM2.5)", 
     main = "Income and particulate emissions across world regions")
abline(model)   # abline can also add model results to a plot 
legend("topright", inset=0.045, legend=c("r =", round(outcorr$estimate,4)))


# A cleaned-up version, starting with the line "plot" above 
plot(x, y, xlab="2019 Personal Income ($USD)", ylab = "Average air pollution level (PM2.5)", 
     main = "Income and particulate emissions across world regions", pch=20, font.main=4, 
     ylim=c(0,40), xlim=c(5000,70000))
abline(model, col="blue", lwd=2)
rp = vector('expression', 3)
rp[1] = substitute(expression(italic(B)[1] == MYVALUE3), list(MYVALUE3=format(summary(model)$coefficients[2], dig=4)))[2]
rp[2] = substitute(expression(italic(p) == MYVALUE2), list(MYVALUE2=format(summary(model)$coefficients[2,4], dig=3)))[2]
rp[3] = substitute(expression(italic(r) == MYVALUE), list(MYVALUE=format(sqrt(summary(model)$r.squared), dig=4)))[2]
legend("topright", legend=rp, bty='n')
text(data$inc19, data$pm2pt5_19-1, labels=data$Country, cex=0.5)

#### DEALING WITH MISSING VALUES #### 
summary(data)
mean(data$broadband19)
mean(data$broadband19, na.rm=T)
data3 = data[!is.na(data$broadband19),]   # isolate rows wherein broadband19 is NA, remove them, and save to 'data3'
nrow(data)
nrow(data3)   # notice how many fewer observations are now in data4 

# Other helpful operations on a data frame 
data$pctsrincr = NULL   # remove this column 
data4 = data[c(1:4)]   # extract just the first 4 columns


#### QUIZ QUESTIONS - CHALLENGE YOURSELF ####
# (3) How many world regions have a 2020 population density over 500/sqkm? 
# (4) What is the correlation coefficient between COVID vaccination rates and broadband access in world? Round to 2 decimal places


########################################################
######## (3) WORKING WITH MULTIPLE DATA SOURCES ########
########################################################

######## MERGE DATASETS FROM DIFFERENT SOURCES AND ANALYZE THEM ####### 
# What predicts domestic migration in California counties? Housing cost, or job growth? 

# Download county migration data which I have pre-processed. Data reflect changes ending 7/1 of year listed. Source: CA Dep't of Finance, http://www.dof.ca.gov/Forecasting/Demographics/Estimates/E-6/ 
dof = read.csv("ca_dof_popchg_e6.csv")
head(dof)
nrow(dof)

# Download Census data the old fashioned way (median home values) 
# Go do data.census.gov then choose Advanced search
# Then, Geography -> Counties -> California -> check all counties in California
# Then, enter B25077 in 'Table ID' and hit SEARCH at the bottom-right which is median value of owner-occupied housing units 
# Select the first entry (it's a Table)
# In the bar above the data, select 'Transpose'
# Then, select 'ZIP' to download ('csv' won't work as a method to prevent you from not downloading metadata) 
# Run this fancy version of read.csv to read in WITHOUT the double-line-header the Bureau puts in. MAKE SURE the file name matches yours! 
acs = read.csv(text=readLines("./ACSDT5Y2020.B25077_data_with_overlays_2022-05-24T172640.csv")[(-2)])
head(acs)

# Job info has been processed for you and is available from the bureau of labor statistics (https://data.bls.gov/cew/apps/data_views/data_views.htm#tab=Tables)
bls = read.csv("blsemp_ca.csv")
head(bls)

# Join other data to 'dof' data using "match" (will work regardless of whether data are sorted, similar to a v-lookup in Excel)
# creates a new column in 'dof' called 'hoval20,' which is equivalent to the column B25077_001E in the acs data
# rows are matched using the unique id in 'dof' (id) and 'acs' (GEO.id2)
dof$hoval20 = acs$B25077_001E[match(dof$id, acs$GEO_ID)]
dof$emp21 = bls$emp_jun21[match(dof$id, bls$GEOID)]   # same thing for bls data 
dof$emp15 = bls$emp_jun15[match(dof$id, bls$GEOID)]   # one time for each row 

# You can make a percent growth variable using basic math:
dof$popgr1721 = (dof$pop21 - dof$pop17)/dof$pop17
dof$empgr1521 = (dof$emp21 - dof$emp15)/dof$emp15

######### ANALYZE THE DATA USING PLOTS AND A REGRESSION MODEL #######

# Compare the size of a county with its median home value. 
# Histograms help you understand how the data are distributed!
hist(dof$pop21, breaks=20) 
# Taking the logarithm of non-zero data can get rid of skew issues 
# Right skew is common in social science data - high values pull up the means. Using a log transformation is particularly useful when comparing data, e.g. understanding the correlation between two variables. 
hist(log(dof$pop21))
boxplot(dof$popgr1721, dof$empgr1521)   # side-by-side boxplot.

# Make a scatterplot and check the correlation
x = dof$hoval20
y = log(dof$pop21)
plot(x,y)
cor.test(x, y)

# Make a fancier scatterplot with a best-fit line
model = lm(y ~ x)  # stores the results of a simple linear regression in the object 'model'
plot(x, y, ylab="Log of Population, 2021", xlab="Median Home Value, 2020", main="California Counties")   # see help(par) to view all plot options
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
# Visit a CRAN mirror page to see all the packages available: https://mirror.las.iastate.edu/CRAN/

# Install packages using a single line of code (or, Tools -> Install Packages)
install.packages("sf")   # a basic mapping package
install.packages("doBy")  # allows for some kinds of data summarization
install.packages("foreign")   # packages allows you to read/write .dbf files

# Each time you use R, you'll need to activate the package.
# A good practice is to include this at the beginning of any script which uses the package (this intro script is called a 'preamble')
library(sf)
library(foreign) 
library(doBy)

#### SUMMARY BY ####  
# A handy trick is to be able to summarize data at a different level, or spatial scale.  
# First, see how many unique values are in the 'mpo' field.
unique(unlist(dof$mpo))
length(unique(unlist(dof$mpo)))

# Return to the OECD data for a moment to see how many regions are in each country
setwd("C:/Users/kane/Dropbox/usc_emup/R_tutorials")
data = read.csv("oecd.csv")

# Make a sorted table and a barplot of the number of regions per country 
table(data$country)
sort(table(data$country), decreasing=TRUE)  # the sort command, with the argument decreasing=TRUE will sort the table in decreasing order.
barplot(sort(table(data$country), decreasing=TRUE), cex.names=0.8,
        main="Number of regions per country")   # a barplot can be made (can also add colors, etc.). cex.names makes the text smaller so all the country names fit.

# Use the summaryBy command (part of the doBy package) to summarize by larger spatial unit 
# a new data frame will be created which summarizes the variables to the left of the tilde at the unit defined to the right of the tilde
dofsum = summaryBy(pop21 + emp21 ~ mpo, data=dof, FUN=sum)   # takes the sum total of emp21 and pop21 by MPO
dofsum
dofmean = summaryBy(hoval20 + immig21 ~ mpo, data=dof, FUN=mean)   # mean across counties in each MPO. Ensure this makes sense!
dofmean

# Write the table output as a .csv to your working directory
write.csv(dofsum, "dofsumtable.csv")

# Make a single barplot by MPO
barplot(dofsum$pop21.sum, names.arg=dofsum$mpo)

# Make the barplot look better
barplot(dofsum$pop21.sum, names.arg=dofmean$mpo, axes=F, col="dodgerblue", main="2021 Population by MPO (millions)")
axis(2, at=seq(0,20000000,5000000), lab=seq(0,20,5))
abline(h=seq(0,20000000,5000000), lty=3, col="darkgrey")
box()

# Summarize some of the OECD data. Note that rate or percentage data doesn't collapse as well. 
out = summaryBy(areasqkm + totpop20 + co2cap08 ~ Country, data=data, FUN=sum)
out2 = out[order(-out$co2cap08.sum),]  # sort by CO2 
out2
cor.test(out2$co2cap08.sum, out2$totpop20.sum)   # are larger countries higher emitters? 

####  REGRESSION ANALYSIS - THE VERY BASICS ######
# the first command (lm) makes a linear model, and stores the output
# the dependent variable goes first, followed by a tilde and the independent variables. 
# for convenience, you can declare the data frame separately and skip the dollar signs
# This shows that net domestic migration in 2021 was a function of smaller population, lower home value, and higher job growth.
m = lm(dommig21 ~ pop21 + hoval20 + empgr1521, data=dof)
summary(m)



#############################################
######## (5) MANIPULATING SHAPEFILES ########  
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
shape$pop21 = dof$pop21[match(shape$GEOID2, dof$id_old)]
shape$hoval20 = dof$hoval20[match(shape$GEOID2, dof$id_old)]
shape$emp21 = dof$emp21[match(shape$GEOID2, dof$id_old)]
summary(shape)

# CAREFULLY write new dbf. Then you'll be able to open in a GIS software and map the new field! 
write.dbf(shape, "cb_2017_ca_county_500k.dbf")



###########################################
####### (6) MAPPING MADE SIMPLE(ISH) ###### 
###########################################
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
plot(counties["pop21"])

# A more advanced version
# Review plot options at https://cran.r-project.org/web/packages/sf/vignettes/sf5.html#geometry_with_attributes:_sf 
plot(counties["pop21"], breaks="jenks", key.pos=2, pal=sf.colors(10), 
     main="2021 Population in California counties", axes=TRUE)



######################################
###### (7) Using the Census API ######
######################################

#### Install TidyCensus & tigris, and load other packages we'll need to use the Census API ####
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
medrent00 = get_decennial(geography = "state", variables = "H060001", sumfile="sf3", year = 2000)

# Make a barplot
barplot(medrent00$value, names.arg=medrent00$NAME, cex.names=0.7, las=2)

# Make a better barplot
medrent00b = medrent00[order(-medrent00$value),]		# sort descending by making a new data frame
barplot(medrent00b$value, names.arg=medrent00b$NAME, cex.names=0.7, las=2,
        main="State median rent, 2000", col="dodgerblue")   #cex.names makes labels smaller, las=2 rotates them.
abline(h=seq(0,700,100), col="darkgrey", lty=3)
box()

#### FINDING GOOD CENSUS VARIABLES TO USE #### 
# Focusing on the most recent ACS 5-year estimates (2020 at the time of this writing) is a good way to go. 
# All the Census API available datasets are at https://api.census.gov/data.html 
# See also, "useful_2020_census_vars_KKguide.xlsx"
# You can also extract a list of variables using the code here.
acsvars = load_variables(2020, "acs5", cache = TRUE)
head(acsvars)
nrow(acsvars)   # tells us how many variables
length(unique(unlist(acsvars$concept)))   # tells us how many unique "concepts" there are
write.csv(acsvars, "acsvars.csv")   # export to .csv so you explore easily in Excel 

#### ASSEMBLE A SET OF TRACT-LEVEL VARIABLES FOR A COUNTY #### 
# Get a single variable (Median rent)
tr = get_acs(geography="tract", state="CA", county="Orange", variables="B25031_001", 
               year=2020, geometry=TRUE)
tr$medrent = tr$estimate  # create a renamed to avoid future confusion
tr$estimate = NULL   # get rid of the old one 

# Add additional variables with a single extraction by putting them into a list. 
# Total population: B00101_001
# Median household income: B19013_001
# Median age of housing: B25035_001

varlist = c("B01001_001", "B19013_001", "B25035_001", "B08101_009")
tr_plus =  get_acs(geography="tract", state="CA", county="Orange", variables=varlist, year=2020)

# Use a match command to bind each new variable to the original data frame (vent). name the variables something intuitive.
# Since the data are stored long, each match operation requires a subset operation first.
tr_a = subset(tr_plus, variable=="B01001_001")
tr$totpop = tr_a$estimate[match(tr$GEOID, tr_a$GEOID)]

tr_b = subset(tr_plus, variable=="B19013_001")
tr$medhhinc = tr_b$estimate[match(tr$GEOID, tr_b$GEOID)]

tr_c = subset(tr_plus, variable=="B25035_001")
tr$medhomeage = tr_c$estimate[match(tr$GEOID, tr_c$GEOID)]

# Remember you can plot this as a map too! 
plot(tr['medhhinc'])

# Finally, we can export this to a shapefile so we can use it elsewhere, too.
st_write(tr, "orange_merge.shp")


##############################################
#### 8.) ALL VARIABLES AT THE TRACT LEVEL ####
##############################################
# I've built this loop to extract all of "Kevin's commonly used variables" in one step. 
# Extract first variable (total population) + geometry
# You can declare your parameters below: 
mystate = "CA"
mycounty = "Imperial"
myyear = 2020
mysurvey = "acs5"
cnty = get_acs(geography="tract", state=mystate, county=mycounty, variables="B01001_001", year=myyear, survey=mysurvey, geometry=TRUE)

# Extract all other variables 
colnames(cnty)[4] = "totpop"
varlist = c('B19001_001', 'B01002_001', 'B03002_001', 'B03002_003', 'B03002_004', 'B03002_012', 'B05001_001', 'B05001_006', 
            'B07001_001', 'B07001_017', 'B07001_033', 'B07001_049', 'B07001_065', 'B08014_001', 'B08014_002', 'B08014_003', 
            'B08014_004', 'B08014_005', 'B08014_006', 'B08014_007', 'B08101_001', 'B08101_009', 'B08101_033', 'B08101_017', 
            'B08101_025', 'B08101_041', 'B08101_049', 'B15003_001', 'B15003_002', 'B15003_003', 'B15003_004', 'B15003_005', 
            'B15003_006', 'B15003_007', 'B15003_008', 'B15003_009', 'B15003_010', 'B15003_011', 'B15003_012', 'B15003_013', 
            'B15003_014', 'B15003_015', 'B15003_016', 'B15003_017', 'B15003_018', 'B15003_019', 'B15003_020', 'B15003_021', 
            'B15003_022', 'B15003_023', 'B15003_024', 'B15003_025', 'B19001_002', 'B19001_003', 'B19001_004', 'B19001_005', 
            'B19001_006', 'B19001_007', 'B19001_008', 'B19001_009', 'B19001_010', 'B19001_011', 'B19001_012', 'B19001_013', 
            'B19001_014', 'B19001_015', 'B19001_016', 'B19001_017', 'B19013_001', 'B23025_001', 'B23025_002', 'B25031_001', 
            'B25077_001', 'B25034_001', 'B25034_002', 'B25034_003', 'B25034_004', 'B25034_005', 'B25034_006', 'B25034_007', 
            'B25034_008', 'B25034_009', 'B25034_010', 'B25034_011', 'B25035_001', 'B25024_002', 'B25024_003', 'B25024_004', 
            'B25024_005', 'B25024_006', 'B25024_007', 'B25024_008', 'B25024_009', 'B25024_010', 'B25024_011', 'B25003_001', 
            'B25003_002', 'B25003_003', 'B01001_003', 'B01001_004', 'B01001_005', 'B01001_006', 'B01001_020', 'B01001_021', 
            'B01001_022', 'B01001_023', 'B01001_024', 'B01001_025', 'B01001_027', 'B01001_028', 'B01001_029', 'B01001_030', 
            'B01001_044', 'B01001_045', 'B01001_046', 'B01001_047', 'B01001_048', 'B01001_049')
cnty2 =  get_acs(geography="tract", state=mystate, county=mycounty, variables=varlist, year=myyear, survey=mysurvey, geometry=TRUE)
varnames = c('totHH', 'medage', 'tot4race', 'wnh', 'bnh', 'hisp', 'tot4_cit', 'noncit', 'tot4_mv', 'no_mv', 
             'in_cnty_mv', 'in_st_mv', 'out_st_mv', 'tot4_vehic', 'vehic0', 'vehic1', 'vehic2', 'vehic3', 'vehic4', 'vehic5pl', 
             'tot4_comm', 'comm_sov', 'comm_walk', 'comm_pool', 'comm_trans', 'comm_oth', 'comm_wah', 'tot4_educ', 'educ_none', 
             'educ_nurs', 'educ_kind', 'educ_1st', 'educ_2nd', 'educ_3rd', 'educ_4th', 'educ_5th', 'educ_6th', 'educ_7th', 'educ_8th', 
             'educ_9th', 'educ_10th', 'educ_11th', 'educ_12th', 'educ_hs', 'educ_ged', 'educ_cnud1', 'educ_some', 'educ_asso', 
             'educ_bach', 'educ_mast', 'educ_prof', 'educ_doct', 'hincund10', 'hinc1014', 'hinc1519', 'hinc2024', 
             'hinc2529', 'hinc3034', 'hinc3539', 'hinc4044', 'hinc4549', 'hinc5059', 'hinc6074', 'hinc7599', 'hinc100124', 'hinc125149',
             'hinc150199', 'hinc200pl', 'medhhinc', 'tot4_work', 'emp', 'medrent', 'medhoval', 'totHU', 'HU2014', 'HU1013', 'HU9', 
             'HU9099', 'HU8089', 'HU7079', 'HU6069', 'HU5059', 'HU4049', 'HU1939', 'medyrblt', 'HU_sfd', 'HU_sfa', 'HU_mf2', 'HU_mf3_4', 
             'HU_mf5_9', 'HU_mf1019', 'HU_mf2049', 'HU_mf50pl', 'HU_mob', 'HU_boatRV', 'tot4_tenu', 'ownerHH', 'renterHH',
             'male0005', 'male0509', 'male1014', 'male1517', 'male6566', 'male6769', 'male7074', 'male7579', 'male8084',
             'male8500', 'female0005','female0509', 'female1014', 'female1517', 'female6566', 'female6769', 'female7074',
             'female7579', 'female8084', 'female8500')

for(i in 1:length(unique(unlist(cnty2$variable)))){
  join = subset(cnty2, variable==varlist[i])
  cnty[,ncol(cnty)+1] = join$estimate[match(cnty$GEOID, join$GEOID)]
  colnames(cnty)[ncol(cnty)] = varnames[i]  
}

# Prep some useful summary variables
cnty$hhinc_sub20 = cnty$hincund10 + cnty$hinc1014 + cnty$hinc1519
cnty$hhinc_over100 = cnty$hinc100124 + cnty$hinc125149 + cnty$hinc150199 + cnty$hinc200pl
cnty$BAplus = cnty$educ_doct + cnty$educ_prof + cnty$educ_mast + cnty$educ_bach
cnty$noHS = cnty$tot4_educ - cnty$BAplus - cnty$educ_asso - cnty$educ_some - cnty$educ_cnud1 - cnty$educ_ged - cnty$educ_hs
cnty$age0017 = cnty$male0005 + cnty$male0509 + cnty$male1014 + cnty$male1517 + cnty$female0005 + cnty$female0509 + cnty$female1014 + cnty$female1517
cnty$age65plus = cnty$male6566 + cnty$male6769 + cnty$male7074 + cnty$male7579 + cnty$male8084 + cnty$male8500 + cnty$female6566 + cnty$female6769 + cnty$female7074 + cnty$female7579 + cnty$female8084 + cnty$female8500

# Carefully write to a .shp
st_write(cnty, "IM_tracts_ACS20.shp")

# Carefully write to a .csv
write.csv(cnty, "IM_tracts_ACS20.csv")

# View a variable
plot(cnty['medhhinc'])


