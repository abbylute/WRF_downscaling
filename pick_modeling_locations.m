function[outlonf, outlatf, outlonc, outlatc] = pick_modeling_locations(ch, outDEMf, outDEMc, outSRf, outSRc, elev_dif_thres, solar_dif_thres, outdir)

% ch = chunk to evaluate
% outDEMf = fine resolution output DEM
% outDEMc = coarse resolution output DEM
% outSRf = fine resolution
% outSRc = coarse resolution
% outdir = output directory


% preallocate output
outlonc = ones(0,1,'single');
outlatc = outlonc;
outelevc = outlonc;

outlonf = outlonc;
outlatf = outlonc;
outelevf = outlonc;



    chunks = matfile([outdir,'chunks/chunk_coordinates_',num2str(outSRf),'m.mat']); chunks = chunks.chunk_coords;
    st_row = chunks.st_row(ch);
    en_row = chunks.en_row(ch);
    st_col = chunks.st_col(ch);
    en_col = chunks.en_col(ch);
    
    
    coarse = matfile(outDEMc);
    elevco = coarse.elev;
    lonco = coarse.lon;
    latco = coarse.lat;
    clear coarse
    
    fine = matfile(outDEMf);
    elevfine = fine.elev(st_row:en_row, st_col:en_col);
    lonfine = fine.lon(st_row:en_row, st_col:en_col);
    latfine = fine.lat(st_row:en_row, st_col:en_col);
    clear fine

    %figure(1);clf;
    %imagesc(elevfine);colorbar();
        
    
    % chunk extent
    chlon = [min(min(lonfine)), max(max(lonfine))];
    chlat = [min(min(latfine)), max(max(latfine))];
    
    
    % crop coarse dem to chunk extent
    f = (lonco >= chlon(1) & lonco <= chlon(2) & latco >= chlat(1) & latco <= chlat(2));
    lonco = lonco(f);
    latco = latco(f);
    elevco = elevco(f);
    
   % figure(2);clf;
   % scatter(lonco,latco,25,elevco,'filled');colorbar();
  
    %figure(3);clf;
    %plot(lonfine(1:10,1:10),latfine(1:10,1:10),'.r'); hold on;
    %plot(lonco([1,2,81,82],1), latco([1,2,81,82],1),'.k');
    
    
    % For each coarse grid cell,
    % - fine the associated fine grid cells
    % - if elevs are all nan (ocean), don't model that cell
    % - quantify topographic complexity
    % - decide to run at coarse or fine resolution
    tol = 10^-8;
    side = ((outSRc/outSRf)-1)/2;
    for cc = 1:size(elevco,1)
        
        % find the fine grid cell that matches the coarse one.
        cent = find(abs(lonfine - lonco(cc)) < tol & abs(latfine - latco(cc)) < tol);
        [r, c] = ind2sub(size(lonfine), cent);
        rpicks = (r-side):(r+side);
        rpicks = rpicks(rpicks>=1 & rpicks <=size(lonfine,1));
        cpicks = (c-side):(c+side);
        cpicks = cpicks(cpicks>=1 & cpicks <=size(lonfine,2));
        fine_elevs = elevfine(rpicks, cpicks);
        %figure(3);clf;imagesc(fine_elevs);colorbar();
        
        % check if in ocean
        if (sum(sum(~isnan(fine_elevs)))==0) % all nans
            % skip to next grid cell, don't model this one
        else

            elevdif = max(max(fine_elevs)) - min(min(fine_elevs));
            
            solar = solarradiation_days(fine_elevs, latfine(rpicks,1), outSRf, .25, 1);
            solarpdif = (max(max(solar))-min(min(solar)))/mean(mean(solar))*100;

            if (elevdif < elev_dif_thres && solarpdif < solar_dif_thres) 
                % append lon/lat of coarse grid cell to list of locations to model
                outlonc = [outlonc; lonco(cc)];
                outlatc = [outlatc; latco(cc)];
                outelevc = [outelevc; elevco(cc)];
            else
                %if (solarpdif < solar_dif_thres) 
                %    % append lon/lat of coarse grid cell to list of locations to model
                %    outlonc = [outlonc; lonco(cc)];
                %    outlatc = [outlatc; latco(cc)];
                %    outelevc = [outelevc; elevco(cc)];
                %else
                   % append lon/lat of all fine grid cells in this coarse grid cell 
                   % to the list of locations to model 
                   lon1 = reshape(lonfine(rpicks,cpicks), length(rpicks)*length(cpicks), 1);
                   lat1 = reshape(latfine(rpicks,cpicks), length(rpicks)*length(cpicks), 1);
                   elev1 = reshape(elevfine(rpicks,cpicks), length(rpicks)*length(cpicks), 1);
                   outlonf = [outlonf; lon1];
                   outlatf = [outlatf; lat1];
                   outelevf = [outelevf; elev1];
                %end
            end
        end % end if all nans
    end % end cells
    
    %figure(2);clf;
    %plot(outlonf,outlatf,'.k'); hold on;
    %plot(outlonc,outlatc,'.r');
    


    % Save lat lon and resolution of each point to model at
    save([outdir,'chunks/points_to_model_chunk_',num2str(ch),'.mat'],'outSRf','outlonf','outlatf','outelevf',...
        'outSRc','outlonc','outlatc','outelevc');
    
end