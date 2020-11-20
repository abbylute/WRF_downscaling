function[paramfilename, tcfilename] = write_solar_paramfile(ch, inDEM, outSR, outDEM, outlon, outlat, outTR, era, solarparamdir, outdir,reggmttz)

tcfilename = [outdir,era,'/ACSWDNB/solartc_',era,'_',num2str(outSR),'m_chunk',num2str(ch),'.mat'];

paramtext = struct();
paramtext.ch=ch; 
%paramtext.inDEM = inDEM; 
paramtext.outSR=outSR;
paramtext.outDEM=outDEM; 
paramtext.outlon=outlon; 
paramtext.outlat=outlat;
paramtext.outTR = outTR; 
paramtext.outfile = tcfilename;
paramtext.finalGMTtz = 0;
paramtext.regionGMTtz = reggmttz;

paramfilename = [solarparamdir,'chunk_',num2str(ch),'_',num2str(outSR),'m.mat'];
save(paramfilename,'-struct','paramtext','-v6');
        
end