function[] = define_spatial_chunks(outDEM, outSR, chunksize, buffer, outdir, us_latlon)

% outDEM = matfile of lon, lat, elev for the finer resolution output DEM
% outSR = output spatial resolution (m), should match resolution of outDEM
% chunksize = size of spatial chunks (m), without buffer
% buffer = size of buffer around each chunk (m)
% outdir = directory for saving chunk data
% us_latlon = location of US latlon coordinates

% load US outline
us = matfile(us_latlon);

% ncell on chunk side
window = chunksize/outSR;
bufcells = ceil(buffer/outSR);

% define chunks using finer resolution outDEM
m = matfile(outDEM);
lonfine = m.lon;
latfine = m.lat;
elevfine = m.elev;
clear m

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
st_col = max(st_col, 1);
st_row_buf = max(st_row_buf, 1);
st_col_buf = max(st_col_buf, 1);

% Visualize chunks
imAlpha = ones(size(elevfine));
imAlpha(isnan(elevfine))=0;
figure(1);clf;
imagesc(elevfine, 'AlphaData',imAlpha); 
colorbar(); hold on;
plot(reshape(repelem(st_col,2),2,length(st_col)), [min(st_row), max(en_row)], 'k');
plot([min(st_col), max(en_col)], reshape(repelem(st_row,2),2,length(st_row)), 'k');
title([num2str(outSR),'m chunk map. chunksize = ',num2str(chunksize/1000),'km. # chunks = ',num2str(length(st_col)*length(st_row)),'.']);
print([outdir,'chunks/chunk_map_',num2str(outSR),'m.png'],'-dpng','-r300');


nr = length(st_row);
nc = length(st_col);

% replicate so that there are coords for each chunk
st_row = repelem(st_row, nc);
en_row = repelem(en_row, nc);
st_row_buf = repelem(st_row_buf, nc);
en_row_buf = repelem(en_row_buf, nc);
st_col = repmat(st_col,1,nr);
en_col = repmat(en_col,1,nr);
st_col_buf = repmat(st_col_buf,1,nr);
en_col_buf = repmat(en_col_buf,1,nr);


% Identify chunks that are completely outside of US
us = matfile(us_latlon);
in_us = ones(size(st_row,2),1)*NaN;
for ii = 1:size(st_row,2)
    ilat = latfine(st_row(ii):en_row(ii),st_col(ii):en_col(ii));
    ilon = lonfine(st_row(ii):en_row(ii),st_col(ii):en_col(ii));
    %s2 = size(ilat,2);
    %s1 = size(ilat,1);
    % test outer points first to save time
    % this should be sufficient, especially if we decrease chunk sizes
    %in = inpolygon(ilon([1,s1],[1,s2]),ilat([1,s1],[1,s2]),us.us_lon,us.us_lat);
    % test outer sides to save time
    %lons = [ilon(1,:):ilon(end,:),ilon(:,1)':ilon(:,end)'];
    %lats = [ilat(1,:),ilat(end,:),ilat(:,1)',ilat(:,end)'];
    % check at regular intervals along each side of the chunk
    lons=[ilon(1,1:100:end), ilon(end,1:100:end), ilon(1:100:end,1)', ilon(1:100:end,end)', ilon(end,end)];
    lats=[ilat(1,1:100:end), ilat(end,1:100:end), ilat(1:100:end,1)', ilat(1:100:end,end)', ilat(end,end)];
 
    in = inpolygon(lons,lats,us.us_lon,us.us_lat);
    %%%%if in == 1
       % in_us(ii) = true;
    %else
    %    in = inpolygon(ilon,ilat,us.us_lon,us.us_lat); % this takes
    %    forever
        in_us(ii) = sum((sum(in))) > 0;
    %end
    
    % also check if outside of bounds:
    if (max(max(ilat))<30.99 || max(max(ilon))>-102.99 || (min(min(ilon))>-105 && max(max(ilat))<31.4))
        in_us(ii) = 0;
    end
end
p = in_us ==0;
plot((st_col(p)+en_col(p))./2, (st_row(p)+en_row(p))./2, '.r','MarkerSize',20);
clear ilat ilon s2 s1 in


% Only save chunk information for chunks within US
in_us = logical(in_us);
%st_row = st_row(in_us);
%en_row = en_row(in_us);
%st_col = st_col(in_us);
%en_col = en_col(in_us);
%st_row_buf = st_row_buf(in_us);
%en_row_buf = en_row_buf(in_us);
%st_col_buf = st_col_buf(in_us);
%en_col_buf = en_col_buf(in_us);


% save outputs
%----------------------------------------
chunk_coords = struct();
chunk_coords.st_row = st_row;
chunk_coords.en_row = en_row;
chunk_coords.st_col = st_col;
chunk_coords.en_col = en_col;
chunk_coords.st_row_buf = st_row_buf;
chunk_coords.en_row_buf = en_row_buf;
chunk_coords.st_col_buf = st_col_buf;
chunk_coords.en_col_buf = en_col_buf;
chunk_coords.in_us = in_us;
save([outdir,'chunks/chunk_coordinates_',num2str(outSR),'m.mat'], 'chunk_coords');

end