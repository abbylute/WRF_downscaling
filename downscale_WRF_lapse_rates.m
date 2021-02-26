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
    clear pts_to_model
    
    varnms = {'ACLWDNB','PREC_ACC_NC','Q2','WIND'};
    
    side = (window-1)/2;
    nmonths = 156;
    yrs = 2000:2013;
    
    % load high res dem
    hidem = matfile(outDEM);
    hilon = single(hidem.lon(chunks.st_row(ch):chunks.en_row(ch), chunks.st_col(ch):chunks.en_col(ch)));
    hilat = single(hidem.lat(chunks.st_row(ch):chunks.en_row(ch), chunks.st_col(ch):chunks.en_col(ch)));
    hielev = single(hidem.elev(chunks.st_row(ch):chunks.en_row(ch), chunks.st_col(ch):chunks.en_col(ch)));
    %inner_coords = [min(min(hilon)), min(min(hilat));
    %                min(min(hilon)), max(max(hilat));
    %                max(max(hilon)), max(max(hilat));
    %                max(max(hilon)), min(min(hilat));
    %                min(min(hilon)), min(min(hilat))];
                        
    xout = reshape(hilon, size(hilon,1)*size(hilon,2),1);
    yout = reshape(hilat, size(hilat,1)*size(hilat,2),1);
    
    hilonl = reshape(hilon, size(hilon,1)*size(hilon,2),1);
    hilatl = reshape(hilat, size(hilat,1)*size(hilat,2),1);


    hilonb = single(hidem.lon(chunks.st_row_buf(ch):chunks.en_row_buf(ch), chunks.st_col_buf(ch):chunks.en_col_buf(ch)));
    hilatb = single(hidem.lat(chunks.st_row_buf(ch):chunks.en_row_buf(ch), chunks.st_col_buf(ch):chunks.en_col_buf(ch)));
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
    wrfelevfine = F(double(xoutc),double(youtc));
    wrfelevfine = interp2(xoutc, youtc, wrfelevfine, xout, yout);
    clear F

    % find wrf points within chunk domain
    in = inpolygon(wrflonl, wrflatl, outer_coords(:,1), outer_coords(:,2));
    nwrf = sum(in);
    fin = single(find(in));
    clear in
    
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
    
     cal = int16(datevec(datetime(2000,10,1,0,0,0):hours(outTR):datetime(2013,9,30,23,0,0)));
     ym = cal(:,1:2);
     ymgrp = int16(findgroups(cellstr([num2str(ym(:,1)), num2str(ym(:,2))])));
     clear ym
     
% --- Downscale each variable ---%

%% Downscale Solar
vartime = tic;
[paramfilename,tcfilename] = write_solar_paramfile(ch, inDEM, outSR, outDEMtif, outlon, outlat, outTR, era, solarparamdir, outdir, reggmttz);
system([pathtoR,' --vanilla ',solartcRscript,' ',paramfilename]);

 downscale_WRF_solar(ch, tcfilename, wrfhdir, outTR, outSR, era, ...
     wrflon, wrflat, wrfminrow, wrfmaxrow, wrfmincol, wrfmaxcol, ...
     outlon, outlat, outdir, xoutc, youtc, xout, yout)
clear paramfilename tcfilename

vtoc = toc(vartime);
disp(['running ACSWDNB took ',num2str(vtoc/60),' minutes.']);

    
% %% Downscale temperature
% varnm = 'T2';
% vartime = tic;
% 
% % create file to save to:
% m = matfile([outdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'],'Writable',true');
%         
%     for yy = 1:length(yrs) % for each yearly file
% 
%         filenm = [wrfhdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_trimmed_',num2str(outTR),'hr_',num2str(yrs(yy)),'.mat'];
%         datall = matfile(filenm);
%         datall = single(datall.outdata((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side),:));
% 
%         ymdh = find(cal(:,1) == yrs(yy));
% 
%         lr = single(ones(nwrf,length(ymdh)))*NaN;
% 
% 
%         % create subsets of lon, lat, elev, land to match the subset used for datall
%         wrflont = wrflon((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side));
%         wrflatt = wrflat((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side));
%         wrfelevt = wrfelev((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side));
%         wrflandt = wrfland((wrfminrow-side):(wrfmaxrow+side), (wrfmincol-side):(wrfmaxcol+side));
% 
%         % Calculate hourly lapse rates
% 
%         for pp = 1:nwrf
%             [rowpicks, colpicks] = find(wrflont == wrflonl(fin(pp)) & wrflatt == wrflatl(fin(pp)));
%             rowpicks = (rowpicks - side):(rowpicks + side);
%             rowpicks = rowpicks(rowpicks > 0 & rowpicks <= size(wrfelevt,1));
%             colpicks = (colpicks - side):(colpicks + side);
%             colpicks = colpicks(colpicks > 0 & colpicks <= size(wrfelevt,2));
%             ncell = length(rowpicks)*length(colpicks);
%             
%             elev = wrfelevt(rowpicks,colpicks);
%             land = wrflandt(rowpicks,colpicks);
%             elev = reshape(elev, ncell, 1);
%             land = reshape(land, ncell, 1);
%             elev = elev(land==1);
%             %X = [ones(length(elev),1) elev];
% 
%            for tt = 1:length(ymdh)
%                 dat = datall(rowpicks,colpicks,tt);
%                 dat = reshape(dat, ncell, 1);
%                 dat = dat(land==1);
%                 % calculate lapse rate
%                 [p,~,mu]=polyfit(elev,dat,1);
%                 lr(pp,tt) = p(1)/mu(2); % per m
% 
%                 %b = X\dat;
%                 %lr(pp,tt) = b(2); % per m
%             end % end months
%         end % end wrf points
%         clear dat elev land X b pp rowpicks colpicks ncell
%         %figure(1);clf;scatter(wrflonl(fin),wrflatl(fin),45,lr(:,1),'filled');colorbar();
% 
%         
%         %--- Interpolate lapse rates to high resolution ---%   
%         lr_fine = single(ones(size(hilon,1),size(hilon,2),length(ymdh)))*NaN;
%         
%         for tt = 1:length(ymdh) % ~ 10 minutes for yy=1
%             F = scatteredInterpolant(double([wrflonl(fin), wrflatl(fin)]), double(lr(:,tt)));
%             lr_f = F(double(xoutc),double(youtc));
%             lr_f = interp2(xoutc, youtc, lr_f, xout, yout);
%             %figure(2);clf; scatter(xout, yout, 24, lr_fine,'filled');colorbar();     
%             lr_fine(:,:,tt) = reshape(lr_f, size(hilon));
%         end
%         clear F lr_f lr tt
% 
% 
%         % trim data to chunk
%         %datall = datall(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,:);
%         datall = datall((side+1):(size(datall,1)-side), (side+1):(size(datall,2)-side),:);
%         
%         % extract desired spatial points: 
%         datlong = single(ones(nwrf, size(ymdh,1))) * NaN;
%         for ii= 1:nwrf
%             datlong(ii,:) = datall((wrfrow(ii)-wrfminrow+1), (wrfcol(ii)-wrfmincol+1), :); 
%         end
%         clear datall ii
% 
%         % preallocate output
%         datdown = single(ones(size(hilon,1), size(hilon,2), size(ymdh,1))) * NaN;
%         
%         % for each time step
%         % - spatially interpolate hourly data
%         % - apply lapse rate correction
% 
%         for tt = 1:size(datlong,2) % for each time step
%             % interpolate raw outTR hourly data to finer spatial resolution
%             F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), double(datlong(:,tt)));
%             dat_f = F(double(xoutc),double(youtc)); 
%             dat_f = interp2(xoutc, youtc, dat_f, xout,yout);
%             
%             %figure(2);clf; scatter(xout, yout, 24, dat_f,'filled');colorbar();
%             dat_fine = reshape(dat_f, size(hilon));
% 
%             % apply lapse rate correction
%             datdown(:,:,tt) = (lr_fine(:,:,tt) .* (hielev-reshape(wrfelevfine,size(hielev)))) + dat_fine;
%             clear F dat_f dat_fine
%         end % end time steps
%         clear datlong lr_fine
%         
%             
%         % extract points to model at
%         datdown = reshape(datdown, size(datdown,1)*size(datdown,2), size(datdown,3));
%         [~,i] = ismember([outlon,outlat], [hilonl,hilatl], 'rows');
%         datdown = datdown(i,:); 
% 
%         % round to desired precision
%         datdown = round(datdown, 1); % to the 10th
% 
%         % save outputs
%         m.datdown(1:size(outlon,1), ymdh) = datdown;
%         
%         clear datdown
%     end % end years
% 
%     %datdown = single(datdown);
% 
%     
% %     % extract points to model at
% %     datdown = reshape(datdown, size(datdown,1)*size(datdown,2), size(datdown,3));
% %     [~,i] = ismember([outlon,outlat], [hilonl,hilatl], 'rows');
% %     datdown = datdown(i,:); 
% % 
% %     % round to desired precision
% %     datdown = round(datdown, 1); % to the 10th
% % 
% % 
% %     % Save downscaled data
% %     savetime = tic;
% %     save([outdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'],'datdown','-v7.3');
% %     stoc = toc(savetime);
%      vtoc = toc(vartime);
%      disp(['running ',char(varnm),' took ',num2str(vtoc/60),' minutes.']);
% %     disp(['saving ',char(varnm),' took ',num2str(stoc/60),' minutes.']);
% 
% % 
% %     
% %     
%% Downscale LW, PPT, Q2, WIND
r1 = 110;
r2 = 609;
c1 = 255;
c2 = 860;

% for vv = 1:length(varnms)
%         vartime = tic;
%         %--- Calculate monthly Lapse Rates for each WRF point in chunk domain ---%
%         % read in monthly data and trim spatially to match trimmed hourly data
% 
%         if vv==1
%             mon = matfile([wrfmdir,char(varnms(vv)),'_',era,'.mat']);
%             mon = single(mon.lw_monthly(r1:r2,c1:c2,:));
%         elseif vv ==2
%             wrfppt = matfile([wrfmdir,char(varnms(vv)),'_monthly_',era,'.mat']);
%             wrfppt = single(wrfppt.pr_monthly(r1:r2,c1:c2,:));
%             % Prep PRISM ppt
%             prism = matfile(prismppt);
%             plonlat = prism.lonlat;
%             prism = prism.prism_ppt_monthly;
%             % interpolate PRISM to WRF grid
%             % we could save time if we did this once ahead of time, then
%             % just read in the PRISM ppt at the WRF grid points
%             mon = ones(size(wrflonl,1), nmonths)*NaN;
%             good = ~isnan(prism(:,1));
%             for mm = 1:nmonths
%                 F = scatteredInterpolant(plonlat(good,:), prism(good,mm));
%                 mon(:,mm) = F(double(wrflonl), double(wrflatl));
%             end
%             
%             mon = reshape(mon, size(wrflon,1), size(wrflon,2), nmonths);
%                      
%             % remove edge effect along US-CAN border
%             for mm = 1:nmonths
%                 f = wrflat>49.05;
%                 monm = mon(:,:,mm);
%                 wrfm = wrfppt(:,:,mm);
%                 monm(f) = wrfm(f);
%                 mon(:,:,mm) = monm;
%             end
%             clear f monm wrfm
%                 
%             %mon(wrflat>49.05) = wrfppt(wrflat>49.05);
%             
%             % infill portions of domain that prism doesn't cover (e.g. near and above canadian border)
%             %mon(isnan(mon)) = wrfppt(isnan(mon)); 
%             
%             % calculate correction factors to apply to hourly data
%             ratppt = single(mon./wrfppt);
%             ratppt(ratppt > 5) = 5;
%             ratppt(isnan(ratppt)) = 1;
%             clear wrfppt plonlat prism F mm good
%         elseif vv == 3
%             mon = matfile([wrfmdir,char(varnms(vv)),'_',era,'.mat']);
%             mon = single(squeeze(mon.QDATA(r1:r2,c1:c2,1,:)));
%         elseif vv == 4
%             mon = matfile([wrfmdir,'vs_monthly_',era,'.mat']);
%             mon = single(mon.vs_monthly(r1:r2,c1:c2,:));
%         end
% 
%         % preallocate
%         lr = single(ones(nwrf, nmonths)) .* NaN;
% 
%         for pp = 1:nwrf
%             [rowpicks, colpicks] = find(wrflon == wrflonl(fin(pp)) & wrflat == wrflatl(fin(pp)));
%             rowpicks = (rowpicks - side):(rowpicks + side);
%             rowpicks = rowpicks(rowpicks > 0 & rowpicks <= size(wrfelev,1));
%             colpicks = (colpicks - side):(colpicks + side);
%             colpicks = colpicks(colpicks > 0 & colpicks <= size(wrfelev,2));
%             ncell = length(rowpicks)*length(colpicks);
%             
%             elev = wrfelev(rowpicks,colpicks);
%             land = wrfland(rowpicks,colpicks);
%             elev = reshape(elev, ncell, 1);
%             land = reshape(land, ncell, 1);
%             elev = elev(land==1);
%             %X = [ones(length(elev),1) elev];
% 
%            for mm = 1:nmonths
%                 dat = mon(rowpicks,colpicks,mm);
%                 dat = reshape(dat, ncell, 1);
%                 dat = dat(land==1);
%                 % calculate lapse rate and y intercept
%                    [p,~,mu]=polyfit(elev,dat,1);
%                    lr(pp,mm) = p(1)/mu(2); % per m
% 
%                 %b = X\dat;
%                 %lr(pp,mm) = b(2); % per m
%             end % end months
%         end % end wrf points
%         clear rowpicks colpicks ncell mon dat elev land X b
%         %figure(1);clf;scatter(wrflonl(fin),wrflatl(fin),25,lr(:,1),'filled');colorbar();
% 
% 
%         %--- Interpolate lapse rates to high resolution ---%   
%         lr_fine = single(ones(size(hilon,1),size(hilon,2),nmonths))*NaN;
%         
%         for mm = 1:nmonths
%             F = scatteredInterpolant(double([wrflonl(fin), wrflatl(fin)]), double(lr(:,mm)));
%             lr_f = F(double(xoutc),double(youtc));
%             lr_f = interp2(xoutc, youtc, lr_f, xout, yout);           
%             %figure(2);clf; scatter(xout, yout, 24, lr_fine,'filled');colorbar();     
%             lr_fine(:,:,mm) = reshape(lr_f, size(hilon));
%         end
%         clear F lr_f mm lr
% 
%        
%         % if ppt, extract spatial points from correction data and
%         % interpolate to output resolution
%         if vv==2 
%             ratpptl = single(ones(nwrf, nmonths)) * NaN;
%             ratppt = ratppt(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,:);
% 
%             for ii = 1:nwrf
%                 ratpptl(ii,:) = ratppt((wrfrow(ii)-wrfminrow+1), (wrfcol(ii)-wrfmincol+1), :);
%             end
%             ratppt_fine = single(ones(size(hilon,1),size(hilon,2),nmonths))*NaN;
%             for mm = 1:nmonths
%                 F = scatteredInterpolant(double([wrflonl(fin), wrflatl(fin)]), double(ratpptl(:,mm)));
%                 rp_f = F(double(xoutc),double(youtc));
%                 rp_f = interp2(xoutc, youtc, rp_f, xout, yout);               
%                 ratppt_fine(:,:,mm) = reshape(rp_f, size(hilon));
%             end
%             clear ratppt ratpptl ii mm F rp_f
%         end
%                 
%         
%         %--- Compute outTR hourly downscaled values using lapse rates ---%
% 
%         % create file to save to:
%         m = matfile([outdir,era,'/',char(varnms(vv)),'/',char(varnms(vv)),'_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'],'Writable',true);
% 
%         % preallocate output
%         %datdown = ones(size(hilon,1), size(hilon,2), size(cal,1)) * NaN;
%   
%     for yy = 1:length(yrs) % for each yearly file
% 
%         filenm = [wrfhdir,era,'/',char(varnms(vv)),'/',char(varnms(vv)),'_',era,'_trimmed_',num2str(outTR),'hr_',num2str(yrs(yy)),'.mat'];
% 
%         yrcalind = find(cal(:,1)==yrs(yy));
%         yrcal = cal(yrcalind,:);
%         ymdh = find(cal(:,1) == yrs(yy));
% 
%         datall = matfile(filenm);
%         datall = single(datall.outdata(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,:));
%         
%         % extract desired spatial points: 
%         datlong = single(ones(nwrf, size(ymdh,1))) * NaN;
%         for ii= 1:nwrf
%             datlong(ii,:) = datall((wrfrow(ii)-wrfminrow+1), (wrfcol(ii)-wrfmincol+1), :); 
%         end
%         clear datall ii
%         
%         % for each time step
%         % - spatially interpolate hourly data
%         % - apply lapse rate correction
%         
%         % if ppt, compute number of timesteps/month with ppt for each
%         % location
%         if vv==2
%             % compute number of timesteps/month with ppt>0
%             mgrp = findgroups(yrcal(:,2));
%             mgrpu = unique(mgrp);
%             nd = single(zeros(size(hilon,1),size(hilon,2),length(mgrpu)));
% 
%             %for m = 1:length(mgrpu)
%             %   nd(:,m) = sum(datlong(:,mgrp==mgrpu(m))>0,2);
%             %end  
%             for tt = 1:size(datlong,2)
%                 % interpolate raw outTR hourly data to finer spatial resolution
%                 F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), double(datlong(:,tt)));
%                 dat_f = F(double(xoutc),double(youtc)); 
%                 dat_f = interp2(xoutc, youtc, dat_f, xout,yout);
%                 dat_fine = reshape(dat_f, size(hilon));
%                 %figure(1);clf;imagesc(dat_fine);colorbar();
%                 dat_fine = dat_fine > 0;
%                 nd(:,:,mgrp(tt)) = nd(:,:,mgrp(tt)) + dat_fine;
%             end
%             clear F dat_f dat_fine mgrpu
%         end
%         
%         % preallocate output
%         datdown = single(ones(size(hilon,1), size(hilon,2), size(ymdh,1))) * NaN;
% 
%         
%         for tt = 1:size(datlong,2) % for each time step
%             % identify month to use for lapse rates
%             mm = ymgrp(yrcalind(tt));
% 
%             % interpolate raw outTR hourly data to finer spatial resolution
%             F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), double(datlong(:,tt)));
%             dat_f = F(double(xoutc),double(youtc)); 
%             dat_f = interp2(xoutc, youtc, dat_f, xout,yout);
%             %figure(2);clf; scatter(xout, yout, 24, dat_f,'filled');colorbar();
% 
%             dat_fine = single(reshape(dat_f, size(hilon)));
%             clear F dat_f
%             
%             % apply bias and lapse rate corrections
%             if vv==2 % if ppt
%                 dat_fine = dat_fine .* ratppt_fine(:,:,mm);
%                 %nds = nd(:,mgrp(tt));
%                 % interpolate # of ppt time steps to fine res
%                 %F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), nds);
%                 %nds = F(xoutc,youtc);
%                 %nds = interp2(xoutc, youtc, nds, xout, yout);
%                 %nds = reshape(nds, size(lr_fine,1),size(lr_fine,2));
%                 nds = nd(:,:,mgrp(tt));
%                 lapse = lr_fine(:,:,mm)./nds;
%                 datdown1 = (lapse .* (hielev-reshape(wrfelevfine,size(hielev)))) + dat_fine; %0.008 
%                 datdown1(dat_fine ==0) = 0;
%                 datdown1(datdown1<0) = 0;
%                 datdown(:,:,tt) = datdown1;
%                 clear nds F lapse datdown1
%             else
%                 datdown(:,:,tt) = (lr_fine(:,:,mm) .* (hielev-reshape(wrfelevfine,size(hielev)))) + dat_fine; %0.008 
%             end
%             
%             clear dat_fine dat_long
%         end % end time steps
%         
%         
%             
%         if vv==3 % Q2
%             datdown(datdown < 0.00001) = 0.00001;
%         elseif ismember(vv, [1,2,4]) % LW, PPT or WIND
%             datdown(datdown<0)=0;
%         end
%     
%         % extract points to model at
%         datdown = reshape(datdown, size(datdown,1)*size(datdown,2), size(datdown,3));
%         [~,i] = ismember([outlon,outlat], [hilonl,hilatl], 'rows');
%         datdown = datdown(i,:); 
% 
%         % round to desired precision
%         if ismember(vv, [2,4]) % for ppt, wind
%             datdown = round(datdown, 1); % to the 10th
%         elseif ismember(vv,3) % for q2
%             datdown = round(datdown,5); % round to 100th of a gram/kg
%         elseif ismember(vv,1) % for longwave
%             datdown = round(datdown, 5, 'significant'); % 5 sig figs
%         end
%     
%         % Save downscaled data
%         m.datdown(1:size(outlon,1),ymdh) = datdown;
% 
%         clear datdown mm
%         %clear mgrp mgrpu 
%         
%     end % end years
%     
%     clear datdown lr_fine
%     vtoc = toc(vartime);
%     disp(['running ',char(varnms(vv)),' took ',num2str(vtoc/60),' minutes.']);
% 
% end % end variables   




    
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
    mon = single(mon.psfc_monthly(r1:r2,c1:c2,:));
    
    % temporal average:
    mon = mean(mon,3);
      
    % calculate lapse rates
    lr = single(ones(nwrf, 1)) .* NaN;
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

        dat = mon(rowpicks,colpicks,1);
        dat = reshape(dat, ncell, 1);
        dat = dat(land==1);
        % calculate lapse rate and y intercept
        [p,~,mu]=polyfit(elev,dat,1);
        lr(pp,1) = p(1)/mu(2); % per m
    end % end wrf points
    clear rowpicks colpicks ncell dat elev land p mu
    %figure(1);clf;scatter(wrflonl(fin),wrflatl(fin),25,lr(:,1),'filled');colorbar();

    
    %--- Interpolate lapse rates to high resolution ---%   
    lr_fine = single(ones(size(hilon,1),size(hilon,2)))*NaN;
    F = scatteredInterpolant(double([wrflonl(fin), wrflatl(fin)]), double(lr(:,1)));
    lr_f = F(double(xoutc),double(youtc));
    lr_f = interp2(xoutc, youtc, lr_f, xout, yout);           
    %figure(2);clf; scatter(xout, yout, 24, lr_fine,'filled');colorbar();     
    lr_fine(:,:) = reshape(lr_f, size(hilon));
    clear F lr_f lr

    % trim to chunk
    montrim = ones(nwrf,1).*NaN;
    for pp = 1:nwrf  
        [rowpicks, colpicks] = find(wrflon == wrflonl(fin(pp)) & wrflat == wrflatl(fin(pp)));
        montrim(pp,:) = mon(rowpicks,colpicks); % datlong
    end
    
    % interpolate to output resolution:
    F = scatteredInterpolant(double(wrflonl(fin)), double(wrflatl(fin)), double(montrim));
    dat_f = F(double(xoutc),double(youtc));
    dat_f = interp2(xoutc, youtc, dat_f, xout, yout);

    dat_fine = single(reshape(dat_f, size(hilon)));
    clear F dat_f
    %datdown = single(datfine);
    
    % preallocate output
    datdown = single(ones(size(hilon,1), size(hilon,2))) * NaN;

    % apply bias and lapse rate corrections
    datdown(:,:) = (lr_fine(:,:) .* (hielev-reshape(wrfelevfine,size(hielev)))) + dat_fine; %0.008 
    clear dat_fine montrim


    % extract points to model at
    datdown = reshape(datdown, size(datdown,1)*size(datdown,2), 1);
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
    disp(['saving ',char(varnm),' took ',num2str(stoc/60),' minutescd ../W.']);


end % end function


