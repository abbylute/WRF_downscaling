function[] = downscale_WRF_solar(ch, tcfilename, wrfhdir, outTR, outSR, era, ...
    wrflon, wrflat, wrfminrow, wrfmaxrow, wrfmincol, wrfmaxcol, ...
    outlon, outlat, outdir, xoutc, youtc, xout, yout)


    varnm = 'ACSWDNB';
    
    % create file to save to
    m = matfile([outdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'],'Writable',true);
    
    % set up cal
    cal = int16(datevec(datetime(2000,10,1,0,0,0):hours(4):datetime(2013,9,30,23,0,0)));
    yrs = unique(cal(:,1));
    
    % load terrain corrections:
    tc = matfile(tcfilename);
    tc = single(tc.solar_tc);
    nsites = size(tc,1);
    
    % tc has terrain corrections for the 15th day of every month (1 to 12)
    % at every outTR hour
    % shift tc months (jan-dec) to match cal (oct-sep)
    nt = 24/outTR; % # of time steps per day (month)
    %tc = [tc(:,(size(tc,2)-(nt*3)+1):size(tc,2)) tc(:,1:(nt*9))];
    %%%%tc = [tc(:,10:12,:) tc(:,1:9,:)];
    
    % set NaN's to 1, since this will create a terrain correction of 1 = no
    % change
    tc(isnan(tc)) = 1;
    % quality check- don't allow super outrageous correction factors
    tc(tc > 5) = 5;
    tc(tc < .1) = .1;
    
    % reshape tc
    %tc = reshape(tc, nsites, size(tc,2)/nt, nt); % sites, month/day, hour
    
     
    % interpolate each hour from mid-month days to all days of the year
    %tcint = single(ones(nsites, size(cal,1)/nt, nt)) * NaN; % site, day, hour
    %dailycal =unique(cal(:,1:3),'rows');
    %for hr = 1:nt
        % interpolate over extended calendar to avoid end effects
       %dat = outernestedfunction(nsites,cal,nt,dailycal,tc,hr); % try using nested functions to save memory. Not 
       % sure if avoiding copies of dat helps though because we have to
       % duplicate tc (by passing it to the function).
       
    %       dat = single(ones(nsites,(size(cal,1)/nt+30)))*NaN;
    %       f = find(dailycal(:,3) == 15);% &  cal(:,4) == (hr-1));

    %       dat(:,f+15) = repmat(tc(:,:,hr),1,13); % 13 = # years
    %       dat(:,1) = tc(:,end,hr);
    %       dat(:,end) = tc(:,1,hr);
    %       dat = fillmissing(dat,'spline',2);

    %   tcint(:,:,hr) = dat(:, 16:(end-15));      
    %end
    %clear tc dailycal hr f dat
        
    % reshape terrain corrections to match cal
    %tcint = reshape(permute(tcint, [1,3,2]), nsites, 4748*nt); % sites, time
    
    % preallocate downscaled data
    %datdown = single(ones(size(tcint)))*NaN; % sites x time
    
    % roughly trim wrf lon/lat to chunk:
    wrflonsub = wrflon(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol, :);
    wrflatsub = wrflat(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol, :);
    wrflonsub = reshape(wrflonsub, size(wrflonsub,1)*size(wrflonsub,2),1);
    wrflatsub = reshape(wrflatsub, size(wrflatsub,1)*size(wrflatsub,2),1);

    % find wrf cells that are closest to each point of interest
    %dists = pdist2(single([outlon outlat]), single([wrflonsub wrflatsub]));
    %for dd=1:size(dists,1) 
    %    fsm(dd) = find(dists(dd,:) == min(dists(dd,:)));
    %end
    
    for yy = 1:14
        ymdh = int16(find(cal(:,1) == yrs(yy)));
        
        % interpolate each hour from mid-month days to all days of the year
        tcint = single(ones(nsites, size(ymdh,1)/nt, nt)) * NaN; % site, day, hour
        dailycal =unique(cal(ymdh,1:3),'rows');
        mos_in_yy = unique(cal(ymdh,2));
        mos_ext_end = mos_in_yy(end)+1; mos_ext_end(mos_ext_end==13) =1;
        mos_ext_st = mos_in_yy(1)-1; mos_ext_st(mos_ext_st==0) = 12;
        for hr = 1:nt
            % interpolate over extended calendar to avoid end effects
           %dat = outernestedfunction(nsites,cal,nt,dailycal,tc,hr); % try using nested functions to save memory. Not 
           % sure if avoiding copies of dat helps though because we have to
           % duplicate tc (by passing it to the function).

               dat = single(ones(nsites,(size(ymdh,1)/nt+30)))*NaN;
               f = find(dailycal(:,3) == 15);% &  cal(:,4) == (hr-1));

               dat(:,f+15) = tc(:,mos_in_yy,hr);%repmat(tc(:,:,hr),1,13); % 13 = # years
               dat(:,1) = tc(:,mos_ext_st,hr);
               dat(:,end) = tc(:,mos_ext_end,hr);
               dat = fillmissing(dat,'spline',2);

           tcint(:,:,hr) = dat(:, 16:(end-15));      
        end
        clear dailycal hr f dat
        
        % reshape terrain corrections to match cal
        tcint = reshape(permute(tcint, [1,3,2]), nsites, size(ymdh,1)); % sites, time

        % preallocate downscaled data
        %datdown = single(ones(nsites,size(ymdh,1)))*NaN; % sites x time

        filenm = [wrfhdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_trimmed_',num2str(outTR),'hr_',num2str(yrs(yy)),'.mat'];
        datall = matfile(filenm);
        datall = single(datall.outdata(wrfminrow:wrfmaxrow, wrfmincol:wrfmaxcol,:));
        datall = reshape(datall, size(datall,1)*size(datall,2),size(datall,3));
        
        dat_fine = single(ones(length(xout),length(ymdh)))*NaN;
        for tt = 1:length(ymdh)
            F = scatteredInterpolant(double(wrflonsub), double(wrflatsub), double(datall(:,tt)));
            dat_f = F(double(xoutc),double(youtc)); 
            dat_fine(:,tt) = interp2(xoutc, youtc, dat_f, xout,yout);
        end
        clear F dat_f
        
        % extract points to model at
        [~,i] = ismember([outlon,outlat], [xout,yout], 'rows');
        dat_fine = dat_fine(i,:); 
        clear i

        % apply terrain correction
        datdown = dat_fine .* tcint;%tcint(:,ymdh);
        clear dat_fine
        %figure(1);clf;scatter(wrflonsub(fsm),wrflatsub(fsm),45,datdown(:,1),'filled');colorbar();

        % round to desired precision
        datdown = round(datdown, 5, 'significant'); % 5 sig figs
        datdown(datdown < 0) = 0;

        % save downscaled data
        m.datdown(1:nsites,ymdh) = datdown;
        
    end % end years
    
    %datdown = single(datdown);

    % round to desired precision
    %datdown = round(datdown, 5, 'significant'); % 5 sig figs
    %datdown(datdown < 0) = 0;

    % save downscaled data
    %savetime = tic;
    %save([outdir,era,'/',char(varnm),'/',char(varnm),'_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'],'datdown','-v7.3');
    %stoc = toc(savetime);
    %disp(['saving ',char(varnm),' took ',num2str(stoc/60),' minutes.']);

end % ACSWDNB


%figure(1);clf;plot(dat); hold on; plot(varout(ymdh,ss));
%figure(1);clf;plot(dat, varout(ymdh,ss),'.');

%function[dat]= outernestedfunction(nsites,cal,nt,dailycal,tc,hr)
%           dat = single(ones(nsites,(size(cal,1)/nt+30)))*NaN;%

%           f = find(dailycal(:,3) == 15);% &  cal(:,4) == (hr-1));
           
%           innernestedfunction
%           function innernestedfunction
%               dat(:,f+15) = repmat(tc(:,:,hr),1,13); % 13 = # years
%               dat(:,1) = tc(:,end,hr);
%               dat(:,end) = tc(:,1,hr);
%               dat = fillmissing(dat,'spline',2);
%           end
%end




