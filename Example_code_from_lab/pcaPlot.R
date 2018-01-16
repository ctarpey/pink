#PLOT INDIVIDUAL PCAs
library(adegenet)

#########PCA functions
get_pca=function(filename="WAKCombinedMasked_Filtered2.gen")
{
  genepop_file=read.genepop(file=filename)
  trans_data_no_missing=scaleGen(genepop_file,NA.method="mean")
  #trans_data_no_missing=scaleGen(genepop_file)
  ind_pca=dudi.pca(trans_data_no_missing,cent=FALSE,scale=FALSE,scannf=FALSE,nf=4)
  return(ind_pca)    
}

#function to plot pca
plot_pca=function(pca_obj=ind_pca_WR_all_snps,pop_names=c("Anvil B","Yako B","Little Togiak R","Agulowak R","Teal C"," Hansen C")
                  , colors=c("dodgerblue","blue","forestgreen","green","firebrick1","firebrick4","purple","pink","gold1"), 
                  pop_samp_sizes=c(45,47,47,48,44,45),pcs_to_plot=c(1,2),xlim=c(-30,20),ylim=c(-20,20),
                  title="",rev_pc1=1,rev_pc2=1,plot_leg="TRUE",cex_leg=1,leg_loc="topright")
{
  eigenvects_to_plot=pca_obj$li
  variance_explained=(pca_obj$eig)/sum(pca_obj$eig)
  
  col_vect=NULL
  for(i in 1:length(colors))
  {
    rep_vec=rep(colors[i],pop_samp_sizes[i])
    col_vect=c(col_vect,rep_vec)
  }
  xlabel=paste("PC",pcs_to_plot[1]," ",round(variance_explained[pcs_to_plot[1]],3)*100,"% of variance",sep="")
  ylabel=paste("PC",pcs_to_plot[2]," ",round(variance_explained[pcs_to_plot[2]],3)*100,"% of variance",sep="")
  
  plot(rev_pc1*eigenvects_to_plot[,pcs_to_plot[1]],rev_pc2*eigenvects_to_plot[,pcs_to_plot[2]],
       bg=col_vect, pch=21, col="black",xlab=xlabel,ylab=ylabel,
       xlim=xlim, ylim=ylim)
  abline(h=0)
  abline(v=0)
  title(title)
  if(plot_leg=="TRUE"){legend(leg_loc,pop_names,cex=cex_leg, fill=colors,bty="n")}
  return(eigenvects_to_plot)
}



#Chum PCA
ind_pca_chumSNPs=get_pca(filename="E:/CHUM/RAD/Filtering/filteredGenos_filteredSamples_MAF05_singletons_80PercGeno_genepop.gen")

pca_all_inds=plot_pca(pca_obj=ind_pca_chumSNPs,xlim=c(-50,50),ylim=c(-50,50),
                      pop_names=c("Eldorado", "Fish", "Holokuk","Kokwok","Nulato","Otter"),
                      pop_samp_sizes=c(48,39,48,46,38,48),
                      colors=c("blue","light blue","purple","yellow","green","light green"),
                      plot_leg="TRUE",cex_leg=0.7,leg_loc="bottomright",pcs_to_plot=c(1,2))