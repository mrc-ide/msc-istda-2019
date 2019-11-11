## Exercises for week six of introduction to R

#1) Make a plot of the mean BMI and Weight against age decade data for the NHANES data set.
#a) Select only the relevant data from the data set
#b) Calculate the mean BMI and Weight for each age decade (you may need to exclude missing data)
#c) Use gather to make the data format long
#d) Make the plot

#2) For each race in the NHANES data set make a facet wrapped plot of poverty against BMI  Plot a linear 
#model to see if there is a positive or negative relationship.
#a) Select the relevant columns
#b) Plot Poverty against BMI
#c) Use geom_smooth to add a linear model
#d) Facet wrap

#3) Make a linear models of the NHANES data set between nested data frames.  Choose three or 
#four variables that you think will predict the weight of individuals the best.  Are your chosen
#variables significant?

#4) Join and tidy up these data_sets (kids.csv, household.csv)
#a) Load in the two csv files
#b) Investigate which is the best join to use on this data set - why?
#c) Clean up the data