pts = matfile('/Volumes/WDPassport/DATA/WRF/Downscaled/coarse_fine_cell_locations.mat');



latlim = [30 49.5];
lonlim = [-125 -104];
states = shaperead('usastatehi', 'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
lat = [states.Lat];
lon = [states.Lon];


figure(1);clf; fig = gcf; fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 12 12];
ax1 = usamap(latlim, lonlim);
plotm(pts.latsf,pts.lonsf,'.','MarkerSize',.001);
plotm(lat,lon,'k','linewidth',.7);

print('/Volumes/WDPassport/DATA/WRF/Downscaled/fine_cell_coverage.png','-dpng','-r300');      
        
        
