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
elev_dif_thres = 10; % maximum elevation difference across outSRf resolution cells within each outSRc cell (m)
aspect_dif_thes = 10; % as above, but for aspect (degrees)


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

define_spatial_chunks(outDEM, outSRf, chunksize, buffer, outdir, us_latlon)


%% For each spatial chunk,

load([outdir,'chunks/chunk_coordinates.mat']);
nchunk = size(chunk_coords.st_col,2);

for ch = 1:nchunk
    
    % define where to model at what resolution
    [outlonf, outlatf, outlonc, outlatc] = pick_modeling_locations(ch, outDEMf, outDEMc, outSRf, outSRc, elev_dif_thres, aspect_dif_thres, outdir);
    
    

    
    % downscale outTR hourly datasets for all variables except solar
    % dependent on spatial scale
    if (length(outlonc)==0 & length(outlonf)>0) % only downscale at fine res
        downscale_WRF_lapse_rates(ch, inDEM, outDEMf, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era)
    elseif (length(outlonc)>0 & length(outlonf)==0) % only downscale at coarse res
        downscale_WRF_lapse_rates(ch, inDEM, outDEMc, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era)
    elseif (length(outlonc)>0 & length(outlonf)>0) % downscale at both
        downscale_WRF_lapse_rates(ch, inDEM, outDEMf, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era)
        downscale_WRF_lapse_rates(ch, inDEM, outDEMc, outTR, window, outdir, wrfhdir, wrfmdir, prismppt, era)
    else
        warning(['chunk ',num2str(ch),' did not have any points to downscale'])
    end
    % need to modify downscale_WRF_lapse_rates to work at multiple output
    % spatial resolutions and to only save downscaled data at points, so
    % output will be 2d not 3d
    
    
    
% - calculate solar terrain corrections in R, save
% write a text file to say what chunk we're on, where to find the datasets,
% where to output them. Have teh R script read in this text file
paramtext = {ch;
            inDEM;
            outDEMf;
            [outdir,era,'/ADSWDNB/solartc_',era,'_chunk',num2str(ch),'.mat']};
csvwrite(outputfilename, paramtext);
% call the R script
!R CMD BATCH path2Rscript.r

% - downscale solar using terrain corrections, save


end

% OUTPUT
% outTR x outSR datasets for each chunk for each variable








