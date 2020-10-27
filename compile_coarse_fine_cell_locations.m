mdir = 'DATA/WRF/downscaled/WUS/chunks/';

for ii = 1:nchunk
   m = matfile([mdir,'points_to_model_chunk_',num2str(ii),'.mat']);
   if ii==1
       lonf = m.outlonf;
       latf = m.outlatf;
       lonc = m.outlonc;
       latc = m.outlatc;
   else
      lonf = [lonf; m.outlonf];
      latf = [latf; m.outlatf];
      lonc = [lonc; m.outlonc];
      latc = [latc; m.outlatc]; 
   end
end


save([mdir,'coarse_fine_cell_locations.mat'],'lonf','latf','lonc','latc');