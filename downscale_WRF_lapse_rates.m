function[] = downscale_WRF_lapse_rates(ch, outSR, inDEM, outDEM, outDEMtif,...
    outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era, resname, solarparamdir, pathtoR, solartcRscript, reggmttz)

% ch = chunk number to run
% outSR = spatial resolution to downscale to (m)
% inDEM = original WRF DEM
% outDEM =  filename for DEM to downscale to
% outTR = output temporal resolution (hrs)
% window = number of WRF grid cells to use for lapse rate downscaling
% outdir = output directory for downscaled data
% wrfhdir = hourly WRF data
% wrfmdir = monthly WRF data
% prismppt = filenames for PRISM ppt to use for bias correction
% era = CTRL or PGW
% resname = 'fine' or 'coarse'

    chunks = matfile([outdir,'chunks/chunk_coordinates_',num2str(outSR),'m.mat']); chunks = chunks.chunk_coords;
    pts_to_model = matfile([outdir,'chunks/points_to_model_chunk_',num2str(ch),'.mat']);
    if isequal(resname,'fine')
        outlon = pts_to_model.outlonf;
        outlat = pts_to_model.outlatf;
    elseif isequal(resname,'coarse')
        outlon = pts_to_model.outlonc;
        outlat = pts_to_model.outlatc;
    else
        error('argument resname is invalid. should be fine or coarse');
    end
    
    varnms = {'ACLWDNB','PREC_ACC_NC','Q2','WIND'};
    
    side = (window-1)/2;
    nmonths = 156;
    yrs = 2000:2013;
    
    % load high res dem
    hidem = matfile(outDEM);
    hilon = hidem.lon(chunks.st_row(ch):chunks.en_row(ch), chunks.st_col(ch):chunks.en_col(ch));
    hilat = hidem.lat(chunks.st_row(ch):chunks.en_row(ch), chunks.st_col(ch):chunks.en_col(ch));
    hielev = hidem.elev(chunks.st_row(ch):chunks.en_row(ch), chunks.st_col(ch):chunks.en_col(ch));
    inner_coords = [min(min(hilon)), min(min(hilat));
                    min(min(hilon)), max(max(hilat));
                    max(max(hilon)), max(max(hilat));
                    max(max(hilon)), min(min(hilat));
                    min(min(hilon)), min(min(hilat))];
                        
    xout = reshape(hilon, size(hilon,1)*size(hilon,2),1);
    yout = reshape(hilat, size(hilat,1)*size(hilat,2),1);
    
    hilonl = reshape(hilon, size(hilon,1)*size(hilon,2),1);
    hilatl = reshape(hilat, size(hilat,1)*size(hilat,2),1);

    hilonb = hidem.lon(chunks.st_row_buf(ch):chunks.en_row_buf(ch), chunks.st_col_buf(ch):chunks.en_col_buf(ch));
    hilatb = hidem.lat(chunks.st_row_buf(ch):chunks.en_row_buf(ch), chunks.st_col_buf(ch):chunks.en_col_buf(ch));
    outer_coords = [min(min(hilonb)), min(min(hilatb));
                    min(min(hilonb)), max(max(hilatb));
                    max(max(hilonb)), max(max(hilatb));
                    max(max(hilonb)), min(min(hilatb));
                    min(min(hilonb)), min(min(hilatb))];
    clear hilonb hilatb 
    
    % create coarse (~4km) res grid to interpolate wrf hourly too first,
    % before fine res
    xoutc = outer_coords(1,1):0.036:outer_coords(3,1);
    youtc = outer_coords(1,2):0.036:outer_coords(2,2);
    xoutc = repmat(xoutc,size(youtc,2),1);
    youtc = repmat(youtc',1,size(xoutc,2));

    % load wrf dem
    wrfdem = matfile(inDEM);
    wrfelev = wrfdem.hgt;
    wrfland = wrfdem.land;
    wrflon = wrfdem.lon; 
    wrflonl = reshape(wrflon, size(wrflon,1)*size(wrflon,2), 1);
    wrflat = wrfdem.lat; 
    wrflatl = reshape(wrflat, size(wrflat,1)*size(wrflat,2), 1);
    
    % create interpolated WRF dem
    F = scatteredInterpolant(double(wrflonl), double(wrflatl), double(reshape(wrfelev,size(wrfelev,1)*size(wrfelev,2),1)));
    wrfelevfine = F(xoutc,youtc);
    wrfelevfine = interp2(xoutc, youtc, wrfelevfine, xout, yout);
        

    % find wrf points within chunk domain
    in = inpolygon(wrflonl, wrflatl, outer_coords(:,1), outer_coords(:,2));
    nwrf = sum(in);
    fin = find(in);
    
    % for reading in hourly nc data later, identify a rectangle to import
    [wrfrow, wrfcol] = ind2sub(size(wrfelev), fin);
    wrfminrow = min(wrfrow); wrfmaxrow = max(wrfrow);
    wrfmincol = min(wrfcol); wrfmaxcol = max(wrfcol);
    
    % check that we are picking the right points:
    %figure(1);clf;
    %plot(outer_coords(:,1),outer_coords(:,2)); hold on;
    %plot(inner_coords(:,1),inner_coords(:,2));
    %plot(wrflonl(in), wrflatl(in), '.');
    
    
    %plot(reshape(xoutc,size(xoutc,1)*size(xoutc,2),1), reshape(youtc,size(youtc,1)*size(youtc,2),1),'.');
    
     cal = datevec(datetime(2000,10,1,0,0,0):hours(outTR):datetime(2013,9,30,23,0,0));
     ym = cal(:,1:2);
     ymgrp = findgroups(cellstr([num2str(ym(:,1)), num2str(ym(:,2))]));
   
% --- Downscale each variable ---%

%% Downscale Solar
vartime = tic;
[paramfilename,tcfilename] = write_solar_paramfile(ch, inDEM, outSR, outDEMtif, outlon, outlat, outTR, era, solarparamdir, outdir, reggmttz);
[status] = system([pathtoR,' --vanilla ',solartcRscript,' ',paramfilename]);

downscale_WRF_solar(ch, tcfilename, wrfhdir, outTR, outSR, era, ...
    wrflon, wrflat, wrfminrow, wrfmaxrow, wrfmincol, wrfmaxcol, ...
    outlon, outlat, outdir, xoutc, youtc, xout, yout)

vtoc = toc(vartime);
disp(['running ACSWDNB took ',num2str(vtoc/60),' minutes.']);

    
%% Downscale temperature
varnm = 'T2';
vartime = tic;

% preallocate output
datdown = ones(size(hilon,1), size(hilon,2), size(cal,1)) * NaN;
        
    for yy = 1:length(yrs) % for each yearly file
      %tic
        filenm = [wrfhdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_trimmed_',num2str(outTR),'hr_',num2str(yrs(yy)),'.mat'];
        datall = matfile(filenm);
        datall = datall.outdata((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side),:);

        ymdh = find(cal(:,1) == yrs(yy));

        lr = ones(nwrf,length(ymdh))*NaN;
    %toc % 4s

% create subsets of lon, lat, elev, land to match the subset used for
% datall
    wrflont = wrflon((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side));
    wrflatt = wrflat((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side));
    wrfelevt = wrfelev((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side));
    wrflandt = wrfland((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side));
    
        % Calculate hourly lapse rates
        %tic
        for pp = 1:nwrf
            [rowpicks, colpicks] = find(wrflont == wrflonl(fin(pp)) & wrflatt == wrflatl(fin(pp)));
            rowpicks = (rowpicks - side):(rowpicks + side);
            rowpicks = rowpicks(rowpicks > 0 & rowpicks <= size(wrfelevt,1));
            colpicks = (colpicks - side):(colpicks + side);
            colpicks = colpicks(colpicks > 0 & colpicks <= size(wrfelevt,2));
            ncell = length(rowpicks)*length(colpicks);
            
            elev = wrfelevt(rowpicks,colpicks);
            land = wrflandt(rowpicks,colpicks);
            elev = reshape(elev, ncell, 1);
            land = reshape(land, ncell, 1);
            elev = elev(land==1);
            X = [ones(length(elev),1) elev];

           for tt = 1:length(ymdh)
                dat = datall(rowpicks,colpicks,tt);
                dat = reshape(dat, ncell, 1);
                dat = dat(land==1);
                % calculate lapse rate
                b = X\dat;
                lr(pp,tt) = b(2); % per m
            end % end months
        end % end wrf points
        clear dat elev land X b pp rowpicks colpicks ncell
        %figure(1);clf;scatter(wrflonl(fin),wrflatl(fin),45,lr(:,1),'filled');colorbar();
%toc % 6.3s
        
        %--- Interpolate lapse rates to high resolution ---%   
        %tic
        lr_fine = ones(size(hilon,1),size(hilon,2),length(ymdh))*NaN;
        
        for tt = 1:length(ymdh) % ~ 10 minutes for yy=1
            F = scatteredInterpolant(double([wrflonl(fin), wrflatl(fin)]), lr(:,tt));
            lr_f = F(xoutc,youtc);
            lr_f = interp2(xoutc, youtc, lr_f, xout, yout);
            %figure(2);clf; scatter(xout, yout, 24, lr_fine,'filled');colorbar();     
            lr_fine(:,:,tt) = reshape(lr_f, size(hilon));
        end
        clear F lr_f lr tt
         %toc % 5s
         
         %tic
        % trim data to chunk
        %datall = datall(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,:);
        datall = datall((side+1):(size(datall,1)-side), (side+1):(size(datall,2)-side),:);
        
        % extract desired spatial points: 
        datlong = ones(nwrf, size(ymdh,1)) * NaN;
        for ii= 1:nwrf
            datlong(ii,:) = datall((wrfrow(ii)-wrfminrow+1), (wrfcol(ii)-wrfmincol+1), :); 
        end
        clear datall ii
       % toc % 0.04s
        
        % for each time step
        % - spatially interpolate hourly data
        % - apply lapse rate correction
        %tic
        for tt = 1:size(datlong,2) % for each time step
            % interpolate raw outTR hourly data to finer spatial resolution
            F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), datlong(:,tt));
            dat_f = F(xoutc,youtc); 
            dat_f = interp2(xoutc, youtc, dat_f, xout,yout);
            
            %figure(2);clf; scatter(xout, yout, 24, dat_f,'filled');colorbar();
            dat_fine = reshape(dat_f, size(hilon));

            % apply lapse rate correction
            datdown(:,:,ymdh(tt)) = (lr_fine(:,:,tt) .* (hielev-reshape(wrfelevfine,size(hielev)))) + dat_fine;
            clear F dat_f dat_fine
        end % end time steps
        clear datlong lr_fine
       % toc % 5.2s
    end % end years
    %tic
    datdown = single(datdown);

    
    % extract points to model at
    datdown = reshape(datdown, size(datdown,1)*size(datdown,2), size(datdown,3));
    [~,i] = ismember([outlon,outlat], [hilonl,hilatl], 'rows');
    datdown = datdown(i,:); 

    % round to desired precision
    datdown = round(datdown, 1); % to the 10th
%toc % 1.5s

    % Save downscaled data
    savetime = tic;
    save([outdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'],'datdown','-v7.3');
    stoc = toc(savetime);
    vtoc = toc(vartime);
    disp(['running ',char(varnm),' took ',num2str(vtoc/60),' minutes.']);
    disp(['saving ',char(varnm),' took ',num2str(stoc/60),' minutes.']);


    
    
%% Downscale LW, PPT, Q2, WIND
r1 = 126;
r2 = 604;
c1 = 286;
c2 = 834;

for vv = 1:length(varnms)
        vartime = tic;
        %--- Calculate monthly Lapse Rates for each WRF point in chunk domain ---%
        % read in monthly data and trim spatially to match trimmed hourly data
       %tic
        if vv==1
            mon = matfile([wrfmdir,char(varnms(vv)),'_',era,'.mat']);
            mon = mon.lw_monthly(r1:r2,c1:c2,:);
        elseif vv ==2
            wrfppt = matfile([wrfmdir,char(varnms(vv)),'_monthly_',era,'.mat']);
            wrfppt = wrfppt.pr_monthly(r1:r2,c1:c2,:);
            % Prep PRISM ppt
            prism = matfile(prismppt);
            plonlat = prism.lonlat;
            prism = prism.prism_ppt_monthly;
            % interpolate PRISM to WRF grid
            % we could save time if we did this once ahead of time, then
            % just read in the PRISM ppt at the WRF grid points
            mon = ones(size(wrflonl,1), nmonths)*NaN;
            good = ~isnan(prism(:,1));
            for mm = 1:nmonths
                F = scatteredInterpolant(plonlat(good,:), prism(good,mm));
                mon(:,mm) = F(double(wrflonl), double(wrflatl));
            end
            
            mon = reshape(mon, size(wrflon,1), size(wrflon,2), nmonths);
                     
            % remove edge effect along US-CAN border
            for mm = 1:nmonths
                f = wrflat>49.05;
                monm = mon(:,:,mm);
                wrfm = wrfppt(:,:,mm);
                monm(f) = wrfm(f);
                mon(:,:,mm) = monm;
            end
            clear f monm wrfm
                
            %mon(wrflat>49.05) = wrfppt(wrflat>49.05);
            
            % infill portions of domain that prism doesn't cover (e.g. near and above canadian border)
            %mon(isnan(mon)) = wrfppt(isnan(mon)); 
            
            % calculate correction factors to apply to hourly data
            ratppt = mon./wrfppt;
            ratppt(ratppt > 5) = 5;
            ratppt(isnan(ratppt)) = 1;
            clear wrfppt plonlat prism F mm good
        elseif vv == 3
            mon = matfile([wrfmdir,char(varnms(vv)),'_',era,'.mat']);
            mon = squeeze(mon.QDATA(r1:r2,c1:c2,1,:));
        elseif vv == 4
            mon = matfile([wrfmdir,'vs_monthly_',era,'.mat']);
            mon = mon.vs_monthly(r1:r2,c1:c2,:);
        end
        %toc % 2s
        %tic
        % preallocate
        lr = ones(nwrf, nmonths) .* NaN;

        for pp = 1:nwrf
            [rowpicks, colpicks] = find(wrflon == wrflonl(fin(pp)) & wrflat == wrflatl(fin(pp)));
            rowpicks = (rowpicks - side):(rowpicks + side);
            rowpicks = rowpicks(rowpicks > 0 & rowpicks <= size(wrfelev,1));
            colpicks = (colpicks - side):(colpicks + side);
            colpicks = colpicks(colpicks > 0 & colpicks <= size(wrfelev,2));
            ncell = length(rowpicks)*length(colpicks);
            
            elev = wrfelev(rowpicks,colpicks);
            land = wrfland(rowpicks,colpicks);
            elev = reshape(elev, ncell, 1);
            land = reshape(land, ncell, 1);
            elev = elev(land==1);
            X = [ones(length(elev),1) elev];

           for mm = 1:nmonths
                dat = mon(rowpicks,colpicks,mm);
                dat = reshape(dat, ncell, 1);
                dat = dat(land==1);
                % calculate lapse rate and y intercept
                b = X\dat;
                lr(pp,mm) = b(2); % per m
            end % end months
        end % end wrf points
        clear rowpicks colpicks ncell mon dat elev land X b
        %figure(1);clf;scatter(wrflonl(fin),wrflatl(fin),25,lr(:,1),'filled');colorbar();
%toc % 2s

%tic
        %--- Interpolate lapse rates to high resolution ---%   
        lr_fine = ones(size(hilon,1),size(hilon,2),nmonths)*NaN;
        
        for mm = 1:nmonths
            F = scatteredInterpolant(double([wrflonl(fin), wrflatl(fin)]), lr(:,mm));
            lr_f = F(xoutc,youtc);
            lr_f = interp2(xoutc, youtc, lr_f, xout, yout);           
            %figure(2);clf; scatter(xout, yout, 24, lr_fine,'filled');colorbar();     
            lr_fine(:,:,mm) = reshape(lr_f, size(hilon));
        end
        clear F lr_f mm lr
       %toc % 1.5s
       
        % if ppt, extract spatial points from correction data and
        % interpolate to output resolution
        if vv==2 
            ratpptl = ones(nwrf, nmonths) * NaN;
            ratppt = ratppt(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,:);

            for ii = 1:nwrf
                ratpptl(ii,:) = ratppt((wrfrow(ii)-wrfminrow+1), (wrfcol(ii)-wrfmincol+1), :);
            end
            ratppt_fine = ones(size(hilon,1),size(hilon,2),nmonths)*NaN;
            for mm = 1:nmonths
                F = scatteredInterpolant(double([wrflonl(fin), wrflatl(fin)]), ratpptl(:,mm));
                rp_f = F(xoutc,youtc);
                rp_f = interp2(xoutc, youtc, rp_f, xout, yout);               
                ratppt_fine(:,:,mm) = reshape(rp_f, size(hilon));
            end
            clear ratppt ratpptl ii mm F rp_f
        end
                
        
        %--- Compute outTR hourly downscaled values using lapse rates ---%

        % preallocate output
        datdown = ones(size(hilon,1), size(hilon,2), size(cal,1)) * NaN;
     %tic   
    for yy = 1:length(yrs) % for each yearly file
        %tic
        filenm = [wrfhdir,era,'/',char(varnms(vv)),'/',char(varnms(vv)),'_',era,'_trimmed_',num2str(outTR),'hr_',num2str(yrs(yy)),'.mat'];

        yrcalind = find(cal(:,1)==yrs(yy));
        yrcal = cal(yrcalind,:);
        ymdh = find(cal(:,1) == yrs(yy));

        datall = matfile(filenm);
        datall = datall.outdata(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,:);
        
        % extract desired spatial points: 
        datlong = ones(nwrf, size(ymdh,1)) * NaN;
        for ii= 1:nwrf
            datlong(ii,:) = datall((wrfrow(ii)-wrfminrow+1), (wrfcol(ii)-wrfmincol+1), :); 
        end
        clear datall ii
        %toc %.65s
        
        % for each time step
        % - spatially interpolate hourly data
        % - apply lapse rate correction
        
        % if ppt, compute number of timesteps/month with ppt for each
        % location
        if vv==2
            % compute number of timesteps/month with ppt>0
            mgrp = findgroups(yrcal(:,2));
            mgrpu = unique(mgrp);
            nd = zeros(size(hilon,1),size(hilon,2),length(mgrpu));

            %for m = 1:length(mgrpu)
            %   nd(:,m) = sum(datlong(:,mgrp==mgrpu(m))>0,2);
            %end  
            for tt = 1:size(datlong,2)
                % interpolate raw outTR hourly data to finer spatial resolution
                F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), datlong(:,tt));
                dat_f = F(xoutc,youtc); 
                dat_f = interp2(xoutc, youtc, dat_f, xout,yout);
                dat_fine = reshape(dat_f, size(hilon));
                %figure(1);clf;imagesc(dat_fine);colorbar();
                dat_fine = dat_fine > 0;
                nd(:,:,mgrp(tt)) = nd(:,:,mgrp(tt)) + dat_fine;
            end
            clear F dat_f dat_fine mgrpu
        end
        
       % tic
        for tt = 1:size(datlong,2) % for each time step
            % identify month to use for lapse rates
            mm = ymgrp(yrcalind(tt));

            % interpolate raw outTR hourly data to finer spatial resolution
            F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), datlong(:,tt));
            dat_f = F(xoutc,youtc); 
            dat_f = interp2(xoutc, youtc, dat_f, xout,yout);
            %figure(2);clf; scatter(xout, yout, 24, dat_f,'filled');colorbar();

            dat_fine = reshape(dat_f, size(hilon));

            % apply bias and lapse rate corrections
            %tic
            if vv==2 % if ppt
                dat_fine = dat_fine .* ratppt_fine(:,:,mm);
                %nds = nd(:,mgrp(tt));
                % interpolate # of ppt time steps to fine res
                %F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), nds);
                %nds = F(xoutc,youtc);
                %nds = interp2(xoutc, youtc, nds, xout, yout);
                %nds = reshape(nds, size(lr_fine,1),size(lr_fine,2));
                nds = nd(:,:,mgrp(tt));
                lapse = lr_fine(:,:,mm)./nds;
                datdown1 = (lapse .* (hielev-reshape(wrfelevfine,size(hielev)))) + dat_fine; %0.008 
                datdown1(dat_fine ==0) = 0;
                datdown1(datdown1<0) = 0;
                datdown(:,:,ymdh(tt)) = datdown1;
                clear nds F lapse datdown1
            else
             %   tic
                datdown(:,:,ymdh(tt)) = (lr_fine(:,:,mm) .* (hielev-reshape(wrfelevfine,size(hielev)))) + dat_fine; %0.008 
            %toc % 0.004s
            end
            clear dat_f dat_fine
            %toc      %    0.003s
        end % end time steps
        %toc     %10s
        clear datlong mgrp mgrpu nd mm
    end % end years
    %toc
    
    datdown = single(datdown);
    
    if vv==3 % Q2
        datdown(datdown < 0.00005) = 0.00005;
    elseif ismember(vv, [1,2,4]) % LW, PPT or WIND
        datdown(datdown<0)=0;
    end
    
    % extract points to model at
    datdown = reshape(datdown, size(datdown,1)*size(datdown,2), size(datdown,3));
    [~,i] = ismember([outlon,outlat], [hilonl,hilatl], 'rows');
    datdown = datdown(i,:); 

    % round to desired precision
    if ismember(vv, [2,4]) % for ppt, wind
        datdown = round(datdown, 1); % to the 10th
    elseif ismember(vv,3) % for q2
        datdown = round(datdown,4); % round to 10th of a gram/kg
    elseif ismember(vv,1) % for longwave
        datdown = round(datdown, 5, 'significant'); % 5 sig figs
    end
    
    % Save downscaled data
    savetime = tic;
    save([outdir,era,'/',char(varnms(vv)),'/',char(varnms(vv)),'_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'],'datdown','-v7.3');
    stoc = toc(savetime); % 81 minutes

    vtoc = toc(vartime);
    disp(['running ',char(varnms(vv)),' took ',num2str(vtoc/60),' minutes.']);
    disp(['saving ',char(varnms(vv)),' took ',num2str(stoc/60),' minutes.']);

    clear datdown lr_fine
end % end variables   




    
%% Downscale PSFC    
% compute temporal average from monthly WRF for each location
% interpolate to output resolution
% extract points to model at
% output should be points x 1 (no temporal variability)

    varnm = 'PSFC';
    vartime = tic;

    % preallocate output
    %datdown = ones(size(hilon,1), size(hilon,2), 1) * NaN;
   
    mon = matfile([wrfmdir,char(varnm),'_',era,'.mat']);
    mon = mon.psfc_monthly(r1:r2,c1:c2,:);
    
    % temporal average:
    mon = mean(mon,3);
    
    % trim to chunk
    for pp = 1:nwrf
        [rowpicks, colpicks] = find(wrflon == wrflonl(fin(pp)) & wrflat == wrflatl(fin(pp)));
        montrim(pp,:) = mon(rowpicks,colpicks);
    end

    % interpolate to output resolution:
    F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), montrim);
    datfine = F(xoutc,youtc);
    datfine = interp2(xoutc, youtc, datfine, xout, yout);

    datdown = single(datfine);
    
    % extract points to model at
    [~,i] = ismember([outlon,outlat], [hilonl,hilatl], 'rows');
    datdown = datdown(i,:); 
    
    % round to desired precision
    datdown = round(datdown/100,1)*100; % round to 10th of hPa or mb
    
    % Save downscaled data
    savetime = tic;
    save([outdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'],'datdown','-v7.3');
    stoc = toc(savetime);
    vtoc = toc(vartime);
    disp(['running ',char(varnm),' took ',num2str(vtoc/60),' minutes.']);
    disp(['saving ',char(varnm),' took ',num2str(stoc/60),' minutes.']);


end % end function


