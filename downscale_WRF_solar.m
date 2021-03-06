function[] = downscale_WRF_solar(ch, tcfilename, wrfhdir, solar_outTR, finaloutTR, outSR, era, ...
    wrflon, wrflat, wrfminrow, wrfmaxrow, wrfmincol, wrfmaxcol, ...
    outlon, outlat, outdir, xoutc, youtc, xout, yout)

    outTR = solar_outTR;
    
    varnm = 'ACSWDNB';
    
    % create file to save to
    m = matfile([outdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'],'Writable',true);
    
    % set up cal
    cal = int16(datevec(datetime(2000,10,1,0,0,0):hours(outTR):datetime(2013,9,30,23,0,0)));
    yrs = unique(cal(:,1));
    calout = int16(datevec(datetime(2000,10,1,0,0,0):hours(finaloutTR):datetime(2013,9,30,23,0,0)));

    % load terrain corrections:
    tc = matfile(tcfilename);
    tc = single(tc.solar_tc);
    nsites = size(tc,1);
    
    % tc has terrain corrections for the 15th day of every month (1 to 12)
    % at every outTR hour
    nt = 24/outTR; % # of time steps per day (month)
    
    % set NaN's to 1, since this will create a terrain correction of 1 = no
    % change
    tc(isnan(tc)) = 1;
    % quality check- don't allow super outrageous correction factors
    tc(tc > 5) = 5;
    tc(tc < .1) = .1;
    
    
    % roughly trim wrf lon/lat to chunk:
    wrflonsub = wrflon(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol, :);
    wrflatsub = wrflat(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol, :);
    wrflonsub = reshape(wrflonsub, size(wrflonsub,1)*size(wrflonsub,2),1);
    wrflatsub = reshape(wrflatsub, size(wrflatsub,1)*size(wrflatsub,2),1);

    % define points to model at
    [~,i] = ismember([outlon,outlat], [xout,yout], 'rows');

    % interpolate each hour from mid-month days to all days of the year
    % to create an interpolation that works for all years, use a year with
    % 366 days
    ymdh = (find(cal(:,1) == 2004));
    tccal = cal(ymdh,:);
    tcint = single(ones(nsites, size(ymdh,1)/nt, nt)) * NaN; % site, day, hour
    dailycal =unique(cal(ymdh,1:3),'rows');
    mos_in_yy = unique(cal(ymdh,2));
    mos_ext_end = mos_in_yy(end)+1; mos_ext_end(mos_ext_end==13) =1;
    mos_ext_st = mos_in_yy(1)-1; mos_ext_st(mos_ext_st==0) = 12;
    tic
    for hr = 1:nt
           dat = single(ones(nsites,(size(ymdh,1)/nt+30)))*NaN;
           f = find(dailycal(:,3) == 15);% &  cal(:,4) == (hr-1));

           dat(:,f+15) = tc(:,mos_in_yy,hr);%repmat(tc(:,:,hr),1,13); % 13 = # years
           dat(:,1) = tc(:,mos_ext_st,hr);
           dat(:,end) = tc(:,mos_ext_end,hr);
           dat = fillmissing(dat,'spline',2);

       tcint(:,:,hr) = dat(:, 16:(end-15));      
    end
    clear dailycal hr f dat mos_in_yy mos_ext_end mos_ext_st
    toc
    % reshape terrain corrections to match cal
    tcint = reshape(permute(tcint, [1,3,2]), nsites, size(ymdh,1)); % sites, time


    
    for yy = 1:14
        ymdh = (find(cal(:,1) == yrs(yy)));
        %ymdhout = (find(calout(:,1) == yrs(yy)));
        
%         % interpolate each hour from mid-month days to all days of the year
%         tcint = single(ones(nsites, size(ymdh,1)/nt, nt)) * NaN; % site, day, hour
%         dailycal =unique(cal(ymdh,1:3),'rows');
%         mos_in_yy = unique(cal(ymdh,2));
%         mos_ext_end = mos_in_yy(end)+1; mos_ext_end(mos_ext_end==13) =1;
%         mos_ext_st = mos_in_yy(1)-1; mos_ext_st(mos_ext_st==0) = 12;
%         tic
%         for hr = 1:nt
%                dat = single(ones(nsites,(size(ymdh,1)/nt+30)))*NaN;
%                f = find(dailycal(:,3) == 15);% &  cal(:,4) == (hr-1));
% 
%                dat(:,f+15) = tc(:,mos_in_yy,hr);%repmat(tc(:,:,hr),1,13); % 13 = # years
%                dat(:,1) = tc(:,mos_ext_st,hr);
%                dat(:,end) = tc(:,mos_ext_end,hr);
%                dat = fillmissing(dat,'spline',2);
% 
%            tcint(:,:,hr) = dat(:, 16:(end-15));      
%         end
%         clear dailycal hr f dat mos_in_yy mos_ext_end mos_ext_st
%         toc
%         % reshape terrain corrections to match cal
%         tcint = reshape(permute(tcint, [1,3,2]), nsites, size(ymdh,1)); % sites, time

        % load raw WRF solar data
        filenm = [wrfhdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_trimmed_',num2str(outTR),'hr_',num2str(yrs(yy)),'.mat'];
        datall = matfile(filenm);
        
        mos_in_yr = unique(cal(ymdh,2));
        
        for mo_ct = 1:length(mos_in_yr)
           picks_full = find(calout(:,1)==yrs(yy) & calout(:,2) == mos_in_yr(mo_ct));
           picks_in_yr = find(cal(ymdh,2) == mos_in_yr(mo_ct));
           month_cal = cal(ymdh(picks_in_yr),:);
           picks_tc = ismember(tccal(:,2:4),month_cal(:,2:4),'rows');
           %find(tccal(:,2)==month_cal(:,2) & tccal(:,3)==month_cal(:,3) & tccal(:,4)==month_cal(:,4));
           st = picks_in_yr(1);
           en = picks_in_yr(end);
           
           outdata = single(datall.outdata(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,st:en));
           outdata = reshape(outdata, size(outdata,1)*size(outdata,2),size(outdata,3));

           dat_fine = ones(size(xout,1),size(outdata,2)).*NaN;
           for ti = 1:size(outdata,2)
                F = scatteredInterpolant(double(wrflonsub), double(wrflatsub), double(outdata(:,ti)));
                dat_f = F(double(xoutc),double(youtc)); 
                dat_fine(:,ti) = interp2(xoutc, youtc, dat_f, xout,yout);
           end
           clear F dat_f outdata

           % extract points to model at
           dat_fine = dat_fine(i,:); 

           % apply terrain correction
           datdown = dat_fine .* tcint(:,picks_tc);%tcint(:,st:en);
           clear dat_fine
           %figure(1);clf;scatter(wrflonsub(fsm),wrflatsub(fsm),45,datdown(:,1),'filled');colorbar();
         
           % aggregate to desired temporal resolution
           datdown = reshape(datdown,size(datdown,1),finaloutTR,size(datdown,2)/finaloutTR);
           %datdown = squeeze(sum(datdown,2)); % this fails when nsites =1
           datdown = reshape(sum(datdown,2),size(datdown,1),size(datdown,3));
           
           % round to desired precision
           datdown = round(datdown, 5, 'significant'); % 5 sig figs
           datdown(datdown < 0) = 0;

           % save downscaled data
           m.datdown(1:nsites,picks_full) = datdown;

        end
        
%         % run this over each temporal chunk to aggregate to (e.g. 4 hours)
%         for tt = 1:(length(ymdh)/finaloutTR)
%             st = finaloutTR*tt-3;
%             en = finaloutTR*tt;
%             
%             outdata = single(datall.outdata(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,st:en));
%             outdata = reshape(outdata, size(outdata,1)*size(outdata,2),size(outdata,3));
%             
%             dat_fine = ones(size(xout,1),size(outdata,2)).*NaN;
%             for ti = 1:finaloutTR
%                 F = scatteredInterpolant(double(wrflonsub), double(wrflatsub), double(outdata(:,ti)));
%                 dat_f = F(double(xoutc),double(youtc)); 
%                 dat_fine(:,ti) = interp2(xoutc, youtc, dat_f, xout,yout);
%             end
%             clear F dat_f outdata
%             
%             % extract points to model at
%             dat_fine = dat_fine(i,:); 
% 
%             % apply terrain correction
%             datdown = dat_fine .* tcint(:,st:en);
%             clear dat_fine
%             %figure(1);clf;scatter(wrflonsub(fsm),wrflatsub(fsm),45,datdown(:,1),'filled');colorbar();
% 
%             % aggregate to desired temporal resolution
%             datdown = sum(datdown,2);
% 
%             % round to desired precision
%             datdown = round(datdown, 5, 'significant'); % 5 sig figs
%             datdown(datdown < 0) = 0;
% 
%             % save downscaled data
%             m.datdown(1:nsites,ymdhout(tt)) = datdown;
% 
%         end
%         
        
%         datall = single(datall.outdata(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,:));
%         datall = reshape(datall, size(datall,1)*size(datall,2),size(datall,3));
%         
%         dat_fine = single(ones(length(xout),length(ymdh)))*NaN;
%         for tt = 1:length(ymdh)
%             F = scatteredInterpolant(double(wrflonsub), double(wrflatsub), double(datall(:,tt)));
%             dat_f = F(double(xoutc),double(youtc)); 
%             dat_fine(:,tt) = interp2(xoutc, youtc, dat_f, xout,yout);
%         end
%         clear F dat_f
%         
%         % extract points to model at
%         [~,i] = ismember([outlon,outlat], [xout,yout], 'rows');
%         dat_fine = dat_fine(i,:); 
%         clear i
% 
%         % apply terrain correction
%         datdown = dat_fine .* tcint;%tcint(:,ymdh);
%         clear dat_fine
%         %figure(1);clf;scatter(wrflonsub(fsm),wrflatsub(fsm),45,datdown(:,1),'filled');colorbar();
% 
%         % aggregate to desired temporal resolution
%         datdown = reshape(datdown,size(datdown,1),finaloutTR, size(datdown,2)/finaloutTR);
%         datdown = squeeze(sum(datdown,2));
%         
%         % round to desired precision
%         datdown = round(datdown, 5, 'significant'); % 5 sig figs
%         datdown(datdown < 0) = 0;
% 
%         % save downscaled data
%         m.datdown(1:nsites,ymdhout) = datdown;
%         
    end % end years

end % ACSWDNB

