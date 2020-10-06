%% WRF Downscaling Master Script

% parameters
outTR = 4; % output temporal resolution (hrs)
outSRf = 210; % fine output spatial resolution (m)
outSRc = 1050; % coarse output spatial resolution (m)
chunksize = outSRf*400; % chunk size without buffer (m). should be a multiple of outSR.
buffer = 16000; % size of buffer around each chunk (m)
    % buffer should allow > than (window-1)/2 WRF grid cells to be
    % available on all sides for lapse rate calculations
window = 7; % number of WRF grid cells to use for lapse rate downscaling
era = 'CTRL'; % time period
elev_dif_thres = 50; % maximum elevation difference across outSRf resolution cells within each outSRc cell (m)
solar_dif_thres = 10; % as above, but for solar radiation from matlab (%)


% directories
mdir = '/home/abby/'; % master directory
wrfhdir = [mdir,'DATA/WRF/WUS_',num2str(outTR),'hr/']; % hourly WRF data
wrfmdir = [mdir,'DATA/WRF/monthly/']; % monthly WRF data
prismppt = [mdir,'DATA/PRISM/Monthly_ppt_4km_Oct2000_Sep2013.mat'];
outdir = [mdir,'DATA/WRF/downscaled/WUS/']; % output directory for downscaled data
inDEM = [wrfhdir,'lon_lat_hgt_trimmed.mat']; % filename for 4km WRF DEM
outDEMf = [mdir,'DATA/Mapping/WUS_NED_',num2str(outSRf),'m.mat']; % filename for fine resolution output DEM
outDEMc = [mdir,'DATA/Mapping/WUS_NED_',num2str(outSRc),'m.mat']; % filename for coarse resolution output DEM
us_latlon = [mdir,'DATA/Mapping/US_latlon.mat']; % location of US latlon
addpath([mdir,'DATA/WRF/downscaled/Code/'])
tmpparamfile = []; % file to store parameters temporarily for solar downscaling


%% Define spatial chunks

% at finer resolution:
define_spatial_chunks(outDEMf, outSRf, chunksize, buffer, outdir, us_latlon)

% at medium resolution:
define_spatial_chunks(outDEMc, outSRc, chunksize, buffer, outdir, us_latlon)


%% For each spatial chunk,

finechunks = matfile([outdir,'chunks/chunk_coordinates_',num2str(outSRf),'m.mat']);
finechunks = finechunks.chunk_coords;
%coarsechunks = matfile([outdir,'chunks/chunk_coordinates_',num2str(outSRc),'m.mat']);
nchunk = size(finechunks.st_col,2);

for ch = 1:nchunk
    
    % define where to model at what resolution
    [outlonf, outlatf, outlonc, outlatc] = pick_modeling_locations(ch, outDEMf, outDEMc, outSRf, outSRc, elev_dif_thres, solar_dif_thres, outdir);
    % AT THIS CHUNK SIZE, HOW MANY CHUNKS ARE COMPLETELY COARSE RES?
    
    
    % downscale outTR hourly datasets dependent on spatial scale
    % only downscale at fine res
    if (isempty(outlonc) && ~isempty(outlonf)) 
        downscale_WRF_lapse_rates(ch, outSRf, inDEM, outDEMf, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era, 'fine')
    
        
    % only downscale at coarse res
    elseif (~isempty(outlonc) && isempty(outlonf)) 
        downscale_WRF_lapse_rates(ch, outSRc, inDEM, outDEMc, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era, 'coarse')
        
        
    % downscale at both res    
    elseif (~isempty(outlonc) && ~isempty(outlonf)) 
        downscale_WRF_lapse_rates(ch, outSRf, inDEM, outDEMf, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era, 'fine')
        downscale_WRF_lapse_rates(ch, outSRc, inDEM, outDEMc, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era, 'coarse')
    else
        warning(['chunk ',num2str(ch),' did not have any points to downscale'])
    end

    
    
            paramtext = {ch; inDEM; outSRf; outDEMf; outlonf; outlatf; outTR;
                [outdir,era,'/ADSWDNB/solartc_',era,'_chunk',num2str(ch),'_fine.mat']};
        csvwrite(outputfilename, paramtext);

    % - calculate solar terrain corrections in R, save
    % write a text file to say what chunk we're on, where to find the datasets,
    % where to output them. Have teh R script read in this text file
    paramtext = {ch;
                inDEM;
                outDEMf;
                outDEMc;
                [outdir,'chunks/points_to_model_chunk_',num2str(ch),'.mat'];
                [outdir,era,'/ADSWDNB/solartc_',era,'_chunk',num2str(ch),'.mat']};
    csvwrite(outputfilename, paramtext);
    % call the R script
    !R CMD BATCH path2Rscript.r

% - downscale solar using terrain corrections, save


end

% OUTPUT
% outTR x outSR datasets for each chunk for each variable








