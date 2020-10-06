# Calculate solar terrain corrections for WRF data
require(raster)
require(insol)
require(lutz) # for timezones
#require(tidyverse)

# read in parameter file
params <- '/home/abby/DATA/WRF/xxxxxxx.csv'

# Inputs:
#demres = params[3,] # resolution of fine res dem (m)
#demfn = params[4,] # fine resolution dem
#lonlat = params[5:6,] # lon and lat to calculate corrections at
#deltat = params[7,] # temporal resolution (hrs)
#outfn = params[8,] # output filename
# to try out on laptop:
demres = 210
demfn = '/Volumes/WDPassport/DATA/DEM/NED/210m/WUS_NED_210m.tif'
lonlat = rbind(c(-120, 41), c(-121, 42))
pts_to_model = lonlat; # this will be different in reality
deltat = 4;


# Output:
# a terrain correction value for each prediction point, each mid-month, each hour
nsites = dim(pts_to_model)[1]
solar_tc <- array(rep(NaN, nsites*12*(24/deltat)), c(nsites, 12, (24/deltat)));   # site, month, hour



# 1. prepare mountainous DEM:
#-------------------------
print('preparing DEMs')
demorig <- raster(demfn)
buf <- .5
demorig <- crop(demorig,c(min(lonlat[,1])-buf, max(lonlat[,1])+buf, min(lonlat[,2])-buf, max(lonlat[,2])+buf))
demfine <- projectRaster(demorig,res=demres,crs="+proj=utm +zone=11 +ellps=WGS84 +units=m +no_defs") # zone 11 is middle of WUS
#writeRaster(dem,paste0(demdir,demres,'m/WUS_',demres,'m_utm.tif'),overwrite=T)
#rm(dem);gc()
#demfine <- raster(paste0(demdir,demres,'m/WUS_',demres,'m_utm.tif')) # takes up slightly less space if not stored in memory

height =  mean(values(demfine),na.rm=T)

# 2. prepare 'flat' DEM: 
#-------------------------
# It is difficult (impossible?) to import WRF DEM into R because it is on an irregular grid.
# Instead, aggregate fine DEM to 4km to use as the flat DEM

demflat <- aggregate(demfine, 
                     round(4000/demres), 
                     fun=mean, 
                     expand=TRUE, na.rm=TRUE)#, 
                     #filename=paste0(demdir,demres,'m/WUS_',demres,'m_utm_agg_to_4km.tif'))
#rm(demflat);gc()
#demflat <- raster(paste0(demdir,demres,'m/WUS_',demres,'m_utm_agg_to_4km.tif'))


# 3. Set Parameters:
#-------------------------
visibility=30 # km
RH=50
tempK=273.15
year=2007 # mean year of WRF time period
day=15
#timeh=12
#buf_m = 10000 # buffer in meters around each site to consider in terrain correction
mlon = mean(lonlat[,1])
mlat = mean(lonlat[,2])
ts = seq(0,23,deltat) # time steps each day

# 4. Run 'insol' solar monthly on mtn and flat dems:
#-------------------------
#print('running solar radiation algorithm over DEMs')
demm <- raster:::as.matrix(demfine)
demf <- raster:::as.matrix(demflat)
dl <- res(demfine)[1]
dl2 <- res(demflat)[1]
cgr <- cgrad(demfine) # compute unit vector normal to every grid cell
cgrf <- cgrad(demflat)

for (mm in 1:12){
  tz <- tz_lookup_coords(mlat, mlon, method = "fast", warn = F)
  tz <- tz_offset(paste0("2007-",mm,"-15"), tz = tz)
  tmz <- tz$utc_offset_h
  
  jd=JDymd(year,mm,day,hour=12) # compute Julian Day from a date
  day1=insol::daylength(mlat,mlon,jd,tmz) # length of daylight
  
  for (tt in 1:length(ts)){ # for each time interval b/n sunrise and sunset:
    srs = ts[tt]
    
    if (srs<day1[1] | srs>day1[2]){ # if hr is during dark, no terrain correction
      solar_tc[, mm, tt] <- 1
    
    } else if (srs>=day1[1] & srs<=day1[2]) { # if hr is during daylight then run solar routine
        
      jd=JDymd(year,mm,day,hour=srs,minute=0)
      sv=sunvector(jd,mlat,mlon,tmz)
      zenith=sunpos(sv)[2]
      
      # for mtn dem:
      Idirdif=insolation(zenith,jd,demm,visibility,RH,tempK,0.002,0.15)
      Idirdif = array(Idirdif,c(dim(Idirdif)[1],dim(Idirdif)[2]/2,2))
      # plot(raster(Idirdif[,,1]))
      hsh=hillshading(cgr,sv) # seems to account for surface slope
      sh=doshade(demm,sv,dl) # terrain shading, 0=shade, 1=sun
      ## direct radiation modified by terrain + diffuse irradiation 
      ## values in J/m2
      Irr=(Idirdif[,,1] * hsh *sh + Idirdif[,,2]) * 3600 * deltat
      # plot(raster(Irr))
      
      # for flat dem:
      Idirdif=insolation(zenith,jd,demf,visibility,RH,tempK,0.002,0.15)
      Idirdif = array(Idirdif,c(dim(Idirdif)[1],dim(Idirdif)[2]/2,2))
      hsh=hillshading(cgrf,sv) # seems to account for surface slope
      sh=doshade(demf,sv,dl2) # terrain shading, 0=shade, 1=sun
      ## direct radiation modified by terrain + diffuse irradiation 
      ## values in J/m2
      Irrflat=(Idirdif[,,1] * hsh *sh + Idirdif[,,2]) * 3600 * deltat
      # plot(raster(Irrflat))
      
      # rasterize outputs and make comparable
      Irr <- raster(Irr)
      crs(Irr) <- crs(demfine)
      extent(Irr) <- extent(demfine)
      Irrflat <- raster(Irrflat)
      crs(Irrflat) <- crs(demflat)
      extent(Irrflat) <- extent(demflat)
      
      Ifnew <- disaggregate(Irrflat,4000/demres)
  
      # calculate ratio and assign it to output
      rat <- Irr/Ifnew
      
      # transform back to lat lon
      rat <- projectRaster(rat, demorig)
      
      # extract points to model at
      rat <- raster::extract(rat, pts_to_model)
      
      # output
      solar_tc[, mm, tt] <- rat
    } # end if during daylight
  } # end daylight hour
} # end month

# reshape solar_tc to be space x time
solar_tc = array(solar_tc, dim=c(nsites, 12*length(ts)))

# save terrain corrections
writeMat(outfn, solar_tc = solar_tc)

