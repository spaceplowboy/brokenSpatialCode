rm(list=ls(all=TRUE))
setwd('~/')
list.files()
library(ggfortify)
library(maptools)
library(rgdal)
library(rgeos)
library(spdep)
library(tidyverse)


#############################################
# Stage 1: Retrieve and load data using API #
#############################################

## API HTTP request components
nass      <- 'https://quickstats.nass.usda.gov/api/api_GET/'
nasskey   <- '?key=XXXXXXXXXXXXXXXXXXXXXXXXXX'
fmt_type  <- '&format=CSVE'
src_desc  <- '&source_desc=SURVEY'
year      <- '&year__GE=1950'
ref_desc  <- '&reference_period_desc=YEAR'
ag_desc  <- '&agg_level_desc=STATE'
comm_desc <- '&commodity_desc=AG LAND'
descript  <- 
    '&short_desc=AGLAND, INCL BUILDINGS - ASSET VALUE, MEASURED IN $ / ACRE'
dom_desc  <- '&domain_desc=TOTAL '
freq_desc <- '&freq_desc=ANNUAL'

## Construct the API HTTP request URL
http.request <- paste(nass, nasskey, fmt_type, src_desc, year, ref_desc,
                      agg_desc, comm_desc, descript, dom_desc, freq_desc,
                      sep = ' ') 

## Issue API HTTP request, download data to CSV file, and load into R workspace
download.file(http.request, destfile="farmlandvalues.csv")
mydata.tbl     <- read_csv(file='farmlandvalues.csv')
mydata.tbl     <- mydata.tbl %>%
    select(state_fips_code, state_alpha, state_name, year, Value) %>%
    rename(val_acre = Value)
mydata2014.tbl <- mydata.tbl %>% 
  filter(year==2014) %>%
  arrange(year, state_fips_code)
mydata2014.tbl



#########################################
# Stage 2: Construct Weighting Matrices #
#########################################

# Convert our shapefile into a useable data frame
states <- readOGR(dsn='.', layer='Lower48tests')
states <- subset(states, STATE_FIPS != 11)
states <- states[order(states@data$STATE_FIPS),] # this creates polygon list, states
    
states.q1 <-poly2nb(state)                 #  contruct neightbours list from polygon list (states)

states.q2 <- nblag(states.q1, maxlag=2) #creates second order lag connectivity 
states.q2 <- nblag_cumul(states.q2)

states.q3 <- nblag(states.q1, maxlag=3) # creates third order lag connectivity 
states.q3 <- nblag_cumul(states.q3)



########################################################
# Stage 3: Generate Chloropleth Map of Farmland Values #
########################################################

# Make the shapefile a useable data frame for mapping
states.points <- fortify(states, region='STATE_FIPS')
states.points <- as_tibble(states.points)

# Join the farmland value data with the map data
outmap.tbl <- inner_join(states.points, mydata2014.tbl,
                        by=c('id' = 'state_fips_code'))
outmap.tbl$Value <- as.numeric(as.character(outmap.tbl$val_acre))
outmap.tbl <- outmap.tbl %>% filter(!is.na(val_acre))

whichstrata <- function(x) {
    if (x <= 1850)                    '<= 1850'
    else if (x > 1850 && x <= 3600)   '1851 - 3600'
    else if (x > 3601 && x <= 5600)   '3601 - 5600'
    else if (x > 5601 && x <= 8500)   '5601 - 8500'
    else if (x >= 8500 && x <= 13700) '8500 - 13700'
}
outmap.tbl$ValStrat <- sapply(outmap.tbl$val_acre, whichstrata)

# Plot the choropleth map!
ggplot(outmap.tbl) +
    aes(x=long, y=lat, group=group, fill=ValStrat) + 
    geom_polygon() +
    geom_path(color="black") +
    coord_map() +
    scale_fill_manual(values=c('#e6ccff', '#b98ae6', '#9052cc',
                              '#6d24b2', '#4e0099')) +
    labs(fill='Values($/acre)') +
    labs(x='Longitude (°W)', y='Latitude (°N)') +
    theme(legend.position='bottom')

ggsave('farmlandmap.png', width=6.4, height=4.2, units='in', dpi=600)



#####################################################
# Stage 4: test Moran's I on farmland values        #
#####################################################

# Order our data table by state_fips_code
mydata.tbl <- mydata.tbl %>%
  arrange(year, state_fips_code)

# Test for 2014
mydata2014.tbl$Value <- 
    as.numeric(as.character(mydata2014.tbl$val_acre))
mm <- moran.test(filter(mydata2014.tbl, year==2014)$Value, states.q1)  
mm

