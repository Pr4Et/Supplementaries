%Function requires AstraToolbox

for t=1:16
    clusteralign_astra_reconstruct_BP(0,90,0,sprintf('Z:\\shared\\ArinaData\\2023Jan30_phages\\2023Jan30_KE13_1_tom2_ring%d_reorder_ali.mrc',t),'Z:\shared\ArinaData\2023Jan30_phages\2023Jan30_KE13_1_tom2_ring1_reorder.rawtlt','',1,300);
    clusteralign_astra_reconstruct_BP(0,90,0,sprintf('Z:\\shared\\ArinaData\\2023Jan30_phages\\2023Jan30_KE13_1_tom2_sect%d_reorder_ali.mrc',t),'Z:\shared\ArinaData\2023Jan30_phages\2023Jan30_KE13_1_tom2_ring1_reorder.rawtlt','',1,300);
end

clusteralign_astra_reconstruct(0,90,0,'Z:\\shared\\ArinaData\\2023Jan30_phages\\2023Jan30_KE13_1_tom2_iCOM_reorder_ali.mrc','Z:\shared\ArinaData\2023Jan30_phages\2023Jan30_KE13_1_tom2_ring1_reorder.rawtlt','',1,300);
