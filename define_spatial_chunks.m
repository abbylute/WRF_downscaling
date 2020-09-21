function[] = define_spatial_chunks(outDEM, outSR, chunksize, buffer, outdir)

% inDEM = matfile of lon, lat, elev for WRF
% outDEM = matfile of lon, lat, elev for the finer resolution output DEM
% outSR = output spatial resolution (m), should match resolution of outDEM
% chunksize = size of spatial chunks (m), without buffer
% buffer = size of buffer around each chunk (m)
% outdir = directory for saving chunk data

% ncell on chunk side
window = chunksize/outSR;
bufcells = ceil(buffer/outSR);

% define chunks using finer resolution outDEM
m = matfile(outDEM);
lonfine = m.lon;
latfine = m.lat;
elevfine = m.elev;
clear m

%nchunk_horz = ceiling(size(lonfine,2)/window);
%nchunk_vert = ceiling(size(lonfine,1)/window);

% define coordinates for each chunk
st_row = 1:window:size(lonfine,1);
en_row = st_row - 1 + window; 
st_col = 1:window:size(lonfine,2);
en_col = st_col -1 + window; 

% add buffer
st_row_buf = st_row - bufcells;
en_row_buf = en_row + bufcells;
st_col_buf = st_col - bufcells;
en_col_buf = en_col + bufcells;

% keep st/en within bounds
en_row = min(en_row, size(lonfine,1));
en_col = min(en_col, size(lonfine,2));
en_row_buf = min(en_row_buf, size(lonfine,1));
en_col_buf = min(en_col_buf, size(lonfine,2)); 
st_row = max(st_row, 1);
st_col = max(st_row, 1);
st_row_buf = max(st_row_buf, 1);
st_col_buf = max(st_col_buf, 1);


nr = length(st_row);
nc = length(st_col);

% replicate so that their are coords for each chunk
st_row = repelem(st_row, nc);
en_row = repelem(en_row, nc);
st_row_buf = repelem(st_row_buf, nc);
en_row_buf = repelem(en_row_buf, nc);
st_col = repmat(st_col,1,nr);
en_col = repmat(en_col,1,nr);
st_col_buf = repmat(st_col_buf,1,nr);
en_col_buf = repmat(en_col_buf,1,nr);



% save outputs
chunk_coords = struct();
chunk_coords.st_row = st_row;
chunk_coords.en_row = en_row;
chunk_coords.st_col = st_col;
chunk_coords.en_col = en_col;
chunk_coords.st_row_buf = st_row_buf;
chunk_coords.en_row_buf = en_row_buf;
chunk_coords.st_col_buf = st_col_buf;
chunk_coords.en_col_buf = en_col_buf;
save([outdir,'chunk_coordinates.mat'], 'chunk_coords');

% alternatively, don't save as a structure
%save([outdir,'chunk_coordinates.mat'],...
%    'st_row','en_row','st_col','en_col',...
%    'st_row_buf','en_row_buf','st_col_buf','en_col_buf');






% DO WE WANT TO SAVE FINE AND COARSE MINI DEMS FOR EACH CHUNK? 
% OR MAYBE JUST INFO ABOUT HOW TO EXTRACT THEM?

%chunknum = 0;
%for rr = 1:nchunk_vert    
%    for cc = 1:nchunk_horz
%        chunknum = chunknum +1;
%        
%        % does a directory exist for this chunk? if not make it
%        
%        
%    end % end cols
%end % end rows




end