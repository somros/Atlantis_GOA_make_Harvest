# Alberto Rovellini
# 1/27/2023
# make mFC vectors based on following range of F:
# F = 0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1, 1.25, 1.5, 2

# Formula for mFC is mfc = 1-exp(-F / 365)

f_vals <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00, 1.25, 1.50, 2.00)
mfc_vals <- 1-exp(-f_vals / 365)
# 0.0000000000 0.0001369769 0.0002739351 0.0004108745 0.0005477951 0.0008215801 0.0013689252 0.0020526849 0.0027359764 0.0034188001 0.0041011562 0.0054644672

mfc_tab <- readLines('data/mfc_vector.prm')
mfc_names <- mfc_tab[grepl('mFC', mfc_tab)]
mfc_vec <- mfc_tab[!grepl('mFC', mfc_tab)]

for(mf in 1:length(f_vals)){
  
  file_label <- gsub('\\.', '', as.character(f_vals[mf]))
  
  newfile <- paste0('mFC_sensitivity/mFC_', file_label, '.prm')
  file.create(newfile)
  
  for(v in 1:length(mfc_vec)){
    
    this_v <- as.numeric(strsplit(mfc_vec[v], ' ')[[1]])
    if(this_v[1] > 0){
      this_v[1] <- mfc_vals[mf]
    }
    
    cat(mfc_names[v], file=newfile, append=TRUE,'\n')
    cat(this_v, file=newfile, append=TRUE, '\n') 
  }
}
