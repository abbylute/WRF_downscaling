function[] = downscale_WRF_lapse_rates(ch, outDEM, outSR, outTR, window, outdir, wrfhdir, wrfmdir)

% ch = chunk number to run
% outDEM =  filename for fine resolution DEM
% outSR = output spatial resolution (m)
% outTR = output temporal resolution (hrs)
% window = number of WRF grid cells to use for lapse rate downscaling
% outdir = output directory for downscaled data
% wrfhdir = hourly WRF data
% wrfmdir = monthly WRF data

%NEED TO BE ABLE TO DO THIS FOR PREIND, HIST, FUTURE
%------x8aadhfayfoadfladfjoiayfafndsfyoayf

    chunks = matfile([outdir,'chunk_coordinates.mat']); chunks = chunks.chunk_coords;
    varnms = {'ACLWDNB','PREC_ACC_NC','PSFC','Q2','WIND'};
    
    side = (window-1)/2;
    nmonths = 156;
    
    % load high res dem
    hidem = matfile(outDEM);
    hilon = hidem.lon(chunks.st_row(ch):chunks.en_row(ch), chunks.st_col(ch):chunks.en_col(ch));
    hilat = hidem.lat(chunks.st_row(ch):chunks.en_row(ch), chunks.st_col(ch):chunks.en_col(ch));
    hielev = hidem.elev(chunks.st_row(ch):chunks.en_row(ch), chunks.st_col(ch):chunks.en_col(ch));
    inner_coords = [min(hilon,[],'all'), min(hilat,[],'all');
                    min(hilon,[],'all'), max(hilat,[],'all');
                    max(hilon,[],'all'), max(hilat,[],'all');
                    max(hilon,[],'all'), min(hilat,[],'all');
                    min(hilon,[],'all'), min(hilat,[],'all')];
    hilon = hidem.lon(chunks.st_row_buf(ch):chunks.en_row_buf(ch), chunks.st_col_buf(ch):chunks.en_col_buf(ch));
    hilat = hidem.lat(chunks.st_row_buf(ch):chunks.en_row_buf(ch), chunks.st_col_buf(ch):chunks.en_col_buf(ch));
    hielev = hidem.elev(chunks.st_row_buf(ch):chunks.en_row_buf(ch), chunks.st_col_buf(ch):chunks.en_col_buf(ch));
    outer_coords = [min(hilon,[],'all'), min(hilat,[],'all');
                    min(hilon,[],'all'), max(hilat,[],'all');
                    max(hilon,[],'all'), max(hilat,[],'all');
                    max(hilon,[],'all'), min(hilat,[],'all');
                    min(hilon,[],'all'), min(hilat,[],'all')];
    
    % load wrf dem
    wrfdem = matfile([wrfmdir,'lon_lat_hgt.mat']);
    wrfelev = wrfdem.hgt;
    wrfland = wrfdem.land;
    wrflon = wrfdem.lon; 
    wrflonl = reshape(wrflon, size(wrflon,1)*size(wrflon,2), 1);
    wrflat = wrfdem.lat; 
    wrflatl = reshape(wrflat, size(wrflat,1)*size(wrflat,2), 1);
    
    % find wrf points within chunk domain
    in = inpolygon(wrflonl, wrflatl, outer_coords(:,1), outer_coords(:,2));
    nwrf = sum(in);
    fin = find(in);
    
    % for reading in hourly nc data later, identify a rectangle to import
    [wrfrow, wrfcol] = ind2sub(wrfelev, in);
    wrfminrow = min(wrfrow); wrfmaxrow = max(wrfrow);
    wrfmincol = min(wrfcol); wrfmaxcol = max(wrfcol);
    
    
    

    % check that we are picking the right points:
    %figure(1);clf;
    %plot(outer_coords(:,1),outer_coords(:,2)); hold on;
    %plot(inner_coords(:,1),inner_coords(:,2));
    %plot(wrflonl(in), wrflatl(in), '.');
    
    % set up nc file dates
    yrmonth = [2000*ones(3,1) (10:12)';
        repelem((2001:2012)',12) repmat((1:12)',12,1);
        2013*ones(9,1) (1:9)'];
    filestym = yrmonth(1:3,size(yrmonth,1),:);
    fileenym = yrmonth(3:3:size(yrmonth,1),:);
    cal = datevec(datetime(2000,10,1,0,0,0):hours(1):datetime(2013,9,30,23,0,0));
    calout = cal(1:outTR:size(cal,1),:);
    
% --- Downscale each variable ---%

    for vv = 1:length(varnms)
        if vv==1
            mon = matfile([wrfmdir,char(varnms(vv)),'_monthly_CTRL.mat']);
            mon = mon.lw_monthly;
        elseif vv ==2
            mon = matfile([wrfmdir,char(varnms(vv)),'_monthly_CTRL.mat']);
            mon = mon.pr_monthly;
        elseif vv == 3
            mon = matfile([wrfmdir,char(varnms(vv)),'_CTRL.mat']);
            mon = mon.psfc_monthly;
        elseif vv == 4
            mon = matfile([wrfmdir,char(varnms(vv)),'_CTRL.mat']);
            mon = squeeze(mon.QDATA);
        elseif vv == 5
            mon = matfile([wrfmdir,'vs_monthly_CTRL.mat']);
            mon = mon.vs_monthly;
        else
            error('vv should be 1 through 5.');
        end
        
        
        %--- Calculate Lapse Rates for each WRF point in chunk domain ---%
        
        % preallocate
        lr = ones(nwrf, nmonths) .* NaN;

        for pp = 1:nwrf
            [rowpicks, colpicks] = find(wrflon == wrflonl(fin(pp)) & wrflat == wrflatl(fin(pp)));
            rowpicks = (rowpicks - side):(rowpicks + side);
            rowpicks = rowpicks(rowpicks > 0 & rowpicks <= size(wrfelev,1));
            colpicks = (colpicks - side):(colpicks + side);
            colpicks = colpicks(colpicks > 0 & colpicks <= size(wrfelev,2));
            ncell = length(rowpicks)*length(colpicks);
            
           for mm = 1:nmonths
                dat = mon(rowpicks,colpicks,mm);
                elev = wrfelev(rowpicks,colpicks);
                land = wrfland(rowpicks,colpicks);

                dat = reshape(dat, ncell, 1);
                elev = reshape(elev, ncell, 1);
                land = reshape(land, ncell, 1);
                dat = dat(land==1);
                elev = elev(land==1);

                % calculate lapse rate and y intercept
                X = [ones(length(elev),1) elev];
                b = X\dat;
                lr(pp,mm) = b(2); % per m
            end % end months
        end % end wrf points
        
        %figure(1);clf;scatter(wrflonl(fin),wrflatl(fin),25,lr(:,1),'filled');colorbar();


        %--- Interpolate lapse rates to high resolution ---%   
        
        xout = reshape(hilon, size(hilon,1)*size(hilon,2),1);
        yout = reshape(hilat, size(hilat,1)*size(hilat,2),1);
        lr_fine = ones(size(hilon,1),size(hilon,2),nmonths)*NaN;
        
        for mm = 1:nmonths
            F = scatteredInterpolant(double([wrflonl(fin), wrflatl(fin)]), lr(:,mm));
            lr_f = F(xout,yout);
            %figure(2);clf; scatter(xout, yout, 24, lr_fine,'filled');colorbar();     
            lr_fine(:,:,mm) = reshape(lr_f, size(hilon));
        end
        clear F xout yout lr_f
        
        
        
        %--- Compute hourly downscaled values using lapse rates ---%

        
        
        if ismember(vv, 1:4) % not wind
            for tt = 1:52 % for each 3 month file
                stym = [num2str(filestym(tt,1)), char(sprintfc('%02i',filestym(tt,2)))];
                enym = [num2str(fileenym(tt,1)), char(sprintfc('%02i',fileenym(tt,2)))];
                filenm = [wrfhdir,char(varnms(vv)),'/wrf2d_d01_CTRL_',char(varnms(vv)),'_',stym,'-',enym,'.nc'];
                
                ym = find(yrmonth(:,1) == filestym(tt,1) & ismember(yrmonth(:,2), filestym(tt,2):fileenym(tt,2)));
                ymdh = find(cal(:,1) == filestym(tt,1) & ismember(cal(:,2), filestym(tt,2):fileenym(tt,2)));
                ymdh_grp = findgroups(cal(ymdh,2));
                
                ncid = netcdf.open(filenm);
                datall = squeeze(netcdf.getVar(ncid, 0, [wrfminrow, wrfmincol, 0], [wrfmaxrow, wrfmaxcol, size(ymdh,1)]));

                % thin to desired temporal resolution:
                datall = datall(:,:,1:outTR:size(datall,3));
                
                % extract desired spatial points:                
                datlong = datall((wrfrow-wrfminrow+1), (wrfcol-wrfmincol+1), :); % NEED TO CHECK, RESHAPE THIs to space x time
                
                % for each time step
                % - spatially interpolate hourly data
                % - apply lapse rate correction
                
                for tt = 1:size(datlong,2)
                    % identify month to use for lapse rates
                    mm = ;% THIS NEEDS TO BE BETWEEN 1 and 156
                    
                    F = scatteredInterpolant(double([wrflonl(fin), wrflatl(fin)]), datlong(:,tt));
                    dat_f = F(xout,yout);
                    %figure(2);clf; scatter(xout, yout, 24, lr_fine,'filled');colorbar();     
                    dat_fine(:,:,tt) = reshape(dat_f, size(hilon));
                    
                    dat(:,:,tt) = (lr_fine(:,:,mm) .* hielev) + dat_fine(:,:,tt);
                    

                end % end time steps
                
                % Save downscaled data
                save([outdir,char(varnms(vv)),'_CTRL_',ch,'.mat'],'dat','-v7.3');
                
            end
        elseif vv==5 % for wind
        end
        
    end % end variables
    
end % end function


