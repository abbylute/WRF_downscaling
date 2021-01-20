# Calculate solar terrain corrections for WRF data
#paramfile <-  '/Users/abbylute/Downloads/chunk_140coarse.mat'

# to run this script from command line:
# Rscript --vanilla pathtothisRscript.R pathtoparamfile.mat

# to run this script from matlab:
#pathtoR = '/Library/Frameworks/R.framework/Resources/bin/Rscript';
#pathtoRscript = '/Volumes/WDPassport/DATA/WRF/Downscaled/Code/get_solar_terrain_corrections.R';
#pathtoparamfile = '/Users/abbylute/Downloads/chunk_140coarse.mat';
#[status] = system([pathtoR,' --vanilla ',pathtoRscript,' ',pathtoparamfile]);


  
  args = commandArgs(trailingOnly=TRUE)
  
  # test if there is at least one argument: if not, return an error
  if (length(args)==0) {
    stop("At least one argument must be supplied (paramfile.mat)", call.=FALSE)
  } else if (length(args)==1) {
    # default output file
    paramfile = args[1]
  }

# so that rgdal will load:
#  dyn.load("/opt/modules/climatology/gdal/3.0.2/lib/libgdal.so") # on thunder
  
require(raster)
require(insol)
require(lutz) # for timezones
require(R.matlab)
require(Rcpp)
require(rgdal)

  # to install additional package to use in R from Matlab on thunder, log onto 
  # matlabuser on thunder, module load R, R, install.packages("package")

# read in parameter file
params <- readMat(paramfile)

# Inputs:
ch = params$ch[1,1]
demres = params$outSR[1,1] # resolution of fine res dem (m)
demfn = params$outDEM[1,1] # fine resolution dem
outlon = params$outlon # lons to save
outlat = params$outlat # lats to save
deltat = params$outTR[1,1] # temporal resolution (hrs)
outfn = params$outfile[1,1] # output filename
outgmttz = params$finalGMTtz[1,1] # time zone GMT of output (of WRF data)
reggmttz = params$regionGMTtz[1,1] # time zone GMT of full modeling domain

tmz = reggmttz

pts_to_model = matrix(c(outlon, outlat), ncol = 2); # lon and lat to save corrections at

# Output:
# a terrain correction value for each prediction point, each mid-month, each hour
nsites = dim(pts_to_model)[1]
solar_tc <- array(rep(NaN, nsites*12*(24/deltat)), c(nsites, 12, (24/deltat)));   # site, month, hour



# 1. prepare mountainous DEM:
#-------------------------
demorig <- raster(demfn)
buf <- .5
demorig <- crop(demorig,c(min(pts_to_model[,1])-buf, max(pts_to_model[,1])+buf, min(pts_to_model[,2])-buf, max(pts_to_model[,2])+buf))
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

# then interpolate the flat dem to fine resolution to avoid spatial chunkiness in the terrain corrections
#demflat <- projectRaster(demflat, demfine, res = res(demfine))
demflat <- disaggregate(demflat, fact=round(4000/demres), method='bilinear')

# 3. Set Parameters:
#-------------------------
visibility=30 # km
RH=50
tempK=273.15
year=2007 # mean year of WRF time period
day=15
#timeh=12
#buf_m = 10000 # buffer in meters around each site to consider in terrain correction
mlon = mean(pts_to_model[,1])
mlat = mean(pts_to_model[,2])
#ts = seq(0,23,deltat) # time steps each day in WRF time
ts = seq(deltat/2, 23, deltat) # use midpoint of each timestep instead of start time

# translate wrf times to local times:
tslocal = ts + reggmttz
tslocal[tslocal<0] <- tslocal[tslocal<0] + 24
ts = tslocal 
# so output will be for WRF hours of interest starting at 0:00am GMT0


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
  #tz <- tz_lookup_coords(mlat, mlon, method = "fast", warn = F)
  #tz <- tz_offset(paste0("2007-",mm,"-15"), tz = tz)
  #tmz <- tz$utc_offset_h
  
  jd=JDymd(year,mm,day,hour=12) # compute Julian Day from a date
  day1=insol::daylength(mlat,mlon,jd,tmz) # length of daylight
  
  for (tt in 1:length(ts)){ # for each time interval b/n sunrise and sunset:
    srs = ts[tt]
    
    if (srs<day1[1] | srs>day1[2]){ # if hr is during dark, no terrain correction
      solar_tc[, mm, tt] <- 1
    
    } else if (srs>=day1[1] & srs<=day1[2]) { # if hr is during daylight then run solar routine
        hr_min = strsplit(as.character(srs),"\\.")
        hr = as.numeric(hr_min[[1]][1])
        minu = as.numeric(hr_min[[1]][2]); minu[is.na(minu)] <- 0;
        
      jd=JDymd(year,mm,day,hour=hr,minute=minu)
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
      
      #Irrflat <- crop(Irrflat, extent(Irr))
      
      #Ifnew <- disaggregate(Irrflat,4000/demres)
      Ifnew <- crop(Irrflat, extent(Irr))
  
      # calculate ratio and assign it to output
      rat <- Irr/Ifnew
      
      # transform back to lat lon
      rat <- projectRaster(rat, demorig)
      
      # extract points to model at
      rat <- raster::extract(rat, pts_to_model)
      
      # output
      solar_tc[, mm, tt] <- rat
    } # end if during daylight
  } # end timesteps
} # end month

# save terrain corrections
writeMat(outfn, solar_tc = solar_tc)

