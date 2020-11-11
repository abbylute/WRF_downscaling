pts = matfile('/Volumes/WDPassport/DATA/WRF/Downscaled/coarse_fine_cell_locations.mat');
lonc=pts.lonc;
latc=pts.latc;
lonf=pts.lonf;
latf=pts.latf;


latlim = [31 49.5];
lonlim = [-125 -104];
states = shaperead('usastatehi', 'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
lat = [states.Lat];
lon = [states.Lon];


% check that points are inside US
infine = inpoly2([lonf,latf],[lon',lat']);
inc = inpoly2([lonc,latc],[lon',lat']);

% texas
tlon = states(12).Lon;
tlat = states(12).Lat;
intexf = inpoly2([lonf,latf],[tlon',tlat']);
intexc = inpoly2([lonc,latc],[tlon',tlat']);

% exclude texas too
inc(intexc==1) = 0;
infine(intexf==1) = 0;


figure(1);clf; fig = gcf; fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 12];
ax1 = usamap(latlim, lonlim);
plotm(latf(infine),lonf(infine),'.','MarkerSize',.001);
plotm(lat,lon,'k','linewidth',1.6);
gridm off
print('/Volumes/WDPassport/DATA/WRF/Downscaled/fine_cell_coverage.png','-dpng','-r400');      
        
        
figure(2);clf; fig = gcf; fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 12];
ax1 = usamap(latlim, lonlim);
plotm(latc(inc),lonc(inc),'.','MarkerSize',.6);
plotm(lat,lon,'k','linewidth',1.6);
gridm off
print('/Volumes/WDPassport/DATA/WRF/Downscaled/coarse_cell_coverage.png','-dpng','-r400');      
        


% What percent of domain is covered by fine vs. coarse?
finearea = (sum(infine).*210*210)/1000/1000; % km2
coarsearea = (sum(inc).*1050*1050)/1000/1000; % km2
totalarea = finearea+coarsearea;
finepercent = (finearea/totalarea)*100
coarsepercent = (coarsearea/totalarea)*100







