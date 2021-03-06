% Create deltas for each WRF downscaling chunk to use to perturb the
% historical downscaled WRF data to preindustrial conditions

% load 23-GCM median delta values
cmip = matfile('/Volumes/WDPassport/DATA/CMIP5/NAMERICA_median_delta.mat');
% variables are 30x70x12
clon = cmip.lon -360;
clat = cmip.lat;
cyrs = cmip.yrs;
g = cmip.g;

% calculate the change in global temperature between pre-industrial and 2000-2013:
% g is the change in mean temperature from a 1850-1879 baseline (pre-industrial) starting at year 1861
dt = mean(g(ismember(cyrs,2000:2013))); % about 1C

% calculate changes in each variable
dhuss = cmip.deltaHUSS * dt; % specific humidity, fractional change so multiply
dlrad = cmip.deltaLRAD * dt; % longwave
dppt = cmip.deltaPR * dt; % precip
dsrad = cmip.deltaSRAD * dt; % solar (W/m2), subtract
dtmax = cmip.deltaTMAX * dt;
dtmin = cmip.deltaTMIN * dt;
dwind = cmip.deltaWIND * dt; % wind (m/s), subtract
% don't bother perturbing psfc

% use delta average temperature, not tmax and tmin
dtmean = (dtmax+dtmin)./2;
clear dtmax dtmin



%%%%% interpolate deltas to chunk grid %%%%%
%===========================================
% load chunk grid
ch = matfile('/Volumes/WDPassport/DATA/WRF/Downscaled/chunk_coordinates_210m.mat');
ch = ch.chunk_coords;
nchunk = size(ch.in_us,1);
% this file contains the row and col coordinates within the dem matfile for
% each chunk
dem = matfile('/Volumes/WDPassport/DATA/DEM/NED/new/WUS_NED_210m.mat');

% find the center lat/lon of each chunk:
latcen = ones(nchunk,1).*NaN;
loncen = ones(nchunk,1).*NaN;
for chu = 1:nchunk
    latcen(chu) = mean(dem.lat(ch.st_row(chu):ch.en_row(chu), ch.st_col(chu):ch.en_col(chu)),'all');
    loncen(chu) = mean(dem.lon(ch.st_row(chu):ch.en_row(chu), ch.st_col(chu):ch.en_col(chu)),'all');
end


dhusswrf  = ones(nchunk,12).*NaN;
dlradwrf  = ones(nchunk,12).*NaN;
dpptwrf   = ones(nchunk,12).*NaN;
dsradwrf  = ones(nchunk,12).*NaN;
dwindwrf  = ones(nchunk,12).*NaN;
dtmeanwrf = ones(nchunk,12).*NaN;
for chu = 1:nchunk
    for mm = 1:12
        dhusswrf(chu,mm)  = interp2(clon, clat, dhuss(:,:,mm), loncen(chu), latcen(chu));
        dlradwrf(chu,mm)  = interp2(clon, clat, dlrad(:,:,mm), loncen(chu), latcen(chu));
        dpptwrf(chu,mm)   = interp2(clon, clat, dppt(:,:,mm), loncen(chu), latcen(chu));
        dsradwrf(chu,mm)  = interp2(clon, clat, dsrad(:,:,mm), loncen(chu), latcen(chu));
        dwindwrf(chu,mm)  = interp2(clon, clat, dwind(:,:,mm), loncen(chu), latcen(chu));
        dtmeanwrf(chu,mm) = interp2(clon, clat, dtmean(:,:,mm), loncen(chu), latcen(chu));
    end
end


% For solar,
% convert the absolute delta to a percent data in order to apply to hourly
% WRF data. Do this by first creating monthly historical solar
% climatologies at chunk locations. Interpolation has to be done
% differently since WRF is not on a regular meshgrid.
msw = matfile('/Volumes/WDPassport/DATA/WRF/Monthly_means/SWDOWN_monthly_CTRL.mat');
msw = msw.sw_monthly; % units = W/m2
wlonlat = matfile('/Volumes/WDPassport/DATA/WRF/Monthly_means/lon_lat_hgt.mat');
wlon = double(wlonlat.lon);
wlat = double(wlonlat.lat);
msw = permute(msw,[2,1,3]);
wlonl = reshape(wlon,size(wlon,1)*size(wlon,2),1);
wlatl = reshape(wlat,size(wlat,1)*size(wlat,2),1);

msw1 = ones(nchunk,12).*NaN;
for mm = 1:12
    F = scatteredInterpolant(wlonl, wlatl, reshape(msw(:,:,mm),size(msw,1)*size(msw,2),1));
    msw1(:,mm) = F(loncen, latcen); 
end

% calculate percent deltas:
dsradwrfp = dsradwrf./msw1;
figure(1);clf;scatter(loncen,latcen,100,dsradwrfp(:,1),'filled');colorbar();
%xlabel('lon');ylabel('lat');title('january fractional change in solar (ratio of delta to historical climatology)');


save('/Volumes/WDPassport/DATA/WRF/Downscaled/cmip5_preindustrial_median_deltas_for_WRF_chunks.mat',...
    'dhusswrf','dlradwrf','dpptwrf','dsradwrfp','dwindwrf','dtmeanwrf','-v7.3')

