### Make Networkin DB ####
Y_DIR<-"/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/data/NETWORKS/Y/Default/INDIVIDUAL_NETWORKS/"
Y_networks<-purrr::map(list.files(Y_DIR,full.names = T),function(x){
  data.table::fread(x)%>%as.data.frame()
},.progress = T)%>%do.call(rbind,.)%>%distinct(.)

map_Y_networks <- clusterProfiler::bitr(
  geneID = Y_networks$KSTAR_ACCESSION,
  fromType = "UNIPROT",
  toType = "SYMBOL",
  ,OrgDb = "org.Hs.eg.db"
)
Y_networks$gene<-plyr::mapvalues(Y_networks$KSTAR_ACCESSION,map_Y_networks$UNIPROT,map_Y_networks$SYMBOL)

colnames(ptm_dbs[[1]])
networkkinY<-data.frame("p_site"=paste0(Y_networks$gene,"_",Y_networks$KSTAR_SITE),
                        "enzyme_genesymbol"=Y_networks$KSTAR_KINASE,
                        "mor"=1,
                        likelihood=1)
saveRDS(networkkinY,"./data/ptm_databases/ptm.networKIN.rds")

### Make Networkin DB ####
Y_DIR<-"/mnt/NAS4/user_data/np-dennis/projects/Proteomics/ProteoLab/data/NETWORKS/Y/Default/INDIVIDUAL_NETWORKS/"
Y_networks<-purrr::map(list.files(Y_DIR,full.names = T),function(x){
  data.table::fread(x)%>%as.data.frame()
},.progress = T)%>%do.call(rbind,.)%>%distinct(.)

map_Y_networks <- clusterProfiler::bitr(
  geneID = Y_networks$KSTAR_ACCESSION,
  fromType = "UNIPROT",
  toType = "SYMBOL",
  ,OrgDb = "org.Hs.eg.db"
)
Y_networks$gene<-plyr::mapvalues(Y_networks$KSTAR_ACCESSION,map_Y_networks$UNIPROT,map_Y_networks$SYMBOL)

colnames(ptm_dbs[[1]])
networkkinY<-data.frame("p_site"=paste0(Y_networks$gene,"_",Y_networks$KSTAR_SITE),
                        "enzyme_genesymbol"=Y_networks$KSTAR_KINASE,
                        "mor"=1,
                        likelihood=1)%>%distinct()
saveRDS(networkkinY,"./data/ptm_databases/ptm.networkinY.rds")
