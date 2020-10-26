mdir = '/DATA/WRF/downscaled/WUS/chunks/';

for ii = 1:nchunks
   m = matfile([mdir,'points_to_model_chunk_',num2str(ii),'.mat']);
   if ii==1
       lonf = m.outlonf;
       latf = m.outlatf;
       lonc = m.outlonc;
       latc = m.outlatc;
   else
      lonf = [lonf, outlonf];
      latf = [latf, outlatf];
      lonc = [lonc, outlonc];
      latc = [latc, outlatc]; 
   end
end


save([mdir,'coarse_fine_cell_locations.mat'],'lonf','latf','lonc','latc');