mdir = '/home/abby/DATA/WRF/downscaled/WUS/chunks/';


finechunks = matfile([mdir,'chunk_coordinates_',num2str(outSRf),'m.mat']);
finechunks = finechunks.chunk_coords;
nchunk = size(finechunks.st_col,2);
    
lonf = ones(1,0);
latf = ones(1,0);
elevf = ones(1,0);
lonc = ones(1,0);
latc = ones(1,0);
elevc = ones(1,0);

for ii = 1:nchunk
    
   if finechunks.in_us(ii)==1 % if inside the us, skip this chunk
       m = matfile([mdir,'points_to_model_chunk_',num2str(ii),'.mat']);

      lonf = [lonf; m.outlonf];
      latf = [latf; m.outlatf];
      elevf = [elevf; m.outelevf];
      lonc = [lonc; m.outlonc];
      latc = [latc; m.outlatc]; 
      elevc = [elevc; m.outelevc];

   end
end


save([mdir,'coarse_fine_cell_locations.mat'],'lonf','latf','elevf','lonc','latc','elevc');