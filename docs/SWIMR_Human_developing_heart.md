load package

    library(SingleCellExperiment)
    library(SWIMR)

We obtained the human developing heart data collected using in situ
sequencing (ISS) technology. This dataset has both single-cell and
spot-level Spatial Transcriptomics (ST) data from similar biological
replicates (i.e., two human embryonic cardiac tissue sample at the same
stage, with justification on checking the correlation among expression
similarities in the original paper).

    load(system.file("extdata",'PCW6.5_1_ISS_st.RData',package = 'SWIMR'))
    #### it contains the true proportion stored in data$trueProp

    spatial_count = apply(data$counts,2,as.numeric)
    rownames(spatial_count) = rownames(data$counts)

    spatial_location = apply(data$spatial,2,as.numeric)
    rownames(spatial_location) = rownames(data$spatial)
    spatial_location = as.data.frame(spatial_location)

    sc_internal_count = read.table(system.file("extdata",'internal_PCW6.5_1_scRNA.cnt.tsv',package = 'SWIMR'))
    sc_internal_count = as.matrix(sc_internal_count)
    sc_internal_meta = read.table(system.file("extdata",'internal_PCW6.5_1_scRNA.mta.tsv',package = 'SWIMR'), sep="\t", header = TRUE)

    sc_ex_count = read.table(system.file("extdata",'external_development_heart.scRNA.processed.cnt.genexrow.tsv',package = 'SWIMR'))
    sc_ex_count = as.matrix(sc_ex_count)
    sc_ex_meta = read.table(system.file("extdata",'external_development_heart.scRNA.processed.mta.tsv',package = 'SWIMR'), sep="\t", header = TRUE)

    sce1 <- SingleCellExperiment(
      assays = list(counts = t(sc_internal_count)),
      colData = sc_internal_meta
    )
    sce2 <- SingleCellExperiment(
      assays = list(counts = sc_ex_count),
      colData = sc_ex_meta
    )
    scList = list(sce1, sce2)

We matched the cell types from the external data to our internal data by
looking at which cell types had the most similar gene expression
patterns in both internal and external scRNA-seq datasets. Since the
correlation between Fibroblast-like (related to larger vascular
development) and (8) Fibroblast-like (outflow tract & valve related) and
the correlation between Fibroblast-like (related to larger vascular
development) and (4) Fribroblast-like (coronary & mediastinal
vasculature related) are 1, we simple combined the two cell types into
one cell type, denoted as (15) Fibroblast-like (outflow tract & valve
related & coronary & mediastinal vasculature related).

    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Ventricular cardiomyocytes')] = "(1) Ventricular cardiomyocytes"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Atrial cardiomyocytes')] = "(7) Artrial cardiomyocytes"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Capillary endothelium')] = "(0) Capillary endothelium"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Epicardial cells')] = "(9) Epicardial cells"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Endothelium / pericytes / adventitia')] = "(10) Endothelium/ pericytes/ adventia"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Smooth muscle cells / fibroblast-like)')] =  "(5) Smooth muscle cells/ fibroblast-like"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Cardiac neural crest cells & Schwann progenitor cells')] = "(14) Cells related to cardiac neural crest"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Myoz2-enriched cardiomyocytes')] = "(12) Cardiomyocytes"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Epicardium-derived cells')] = "(3) Subepicardial cells"## not sure about this


    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Fibroblast-like (related to cardiac skeleton connective tissue)')] = "(2) Fribroblast-like (AV mesenchyme related)"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Fibroblast-like (related to larger vascular development)')] = "(8) Fibroblast-like (outflow tract & valve related)"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=='Fibroblast-like (related to smaller vascular development)')] = "(4) Fribroblast-like (coronary & mediastinal vasculature related)"

    ## we combine smaller vascular and larger vascular because the correlation = 1
    colData(scList[[1]])[,2][which(colData(scList[[1]])[,2]=="(8) Fibroblast-like (outflow tract & valve related)")] = "(15) Fibroblast-like (outflow tract & valve related & coronary & mediastinal vasculature related)"
    colData(scList[[1]])[,2][which(colData(scList[[1]])[,2]=="(4) Fribroblast-like (coronary & mediastinal vasculature related)")] = "(15) Fibroblast-like (outflow tract & valve related & coronary & mediastinal vasculature related)"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=="(8) Fibroblast-like (outflow tract & valve related)")] = "(15) Fibroblast-like (outflow tract & valve related & coronary & mediastinal vasculature related)"
    colData(scList[[2]])[,2][which(colData(scList[[2]])[,2]=="(4) Fribroblast-like (coronary & mediastinal vasculature related)")] = "(15) Fibroblast-like (outflow tract & valve related & coronary & mediastinal vasculature related)"

We futher constructed lists as input for SWIMR.

    colData(scList[[1]])$Sample = 1
    colData(scList[[2]])$Sample = 1

    #### start SWIMR
    scRef = list("internal"=list("countData"=counts(scList[[1]]),"metaData"=colData(scList[[1]])),
                 "external"=list("countData"=counts(scList[[2]]),"metaData"=colData(scList[[2]])))

We performed deconvolution using SWIMR.

    SWIMR_obj = createSWIMRObject(
        scRNAseq=scRef,
        spatial_count = spatial_count,
        spatial_location = spatial_location,
        ct.varname_list = list("bio_celltype","bio_celltype"),
        ct.select_list =  list(unique(colData(scList[[1]])$bio_celltype),unique(colData(scList[[2]])$bio_celltype)),
        sample.varname_list = list("Sample","Sample"),
        minCountGene = 100,
        minCountSpot = 5
    )

    SWIMR_obj = SWIMR_deconvolution(SWIMR_obj)

    weight = SWIMR_obj@Weight_est
    proportion = SWIMR_obj@Proportion_est

We checked the root mean squared error (RMSE) for SWIMR and visulized
the weight distribution among two references and the cell type
proportion obtained from SWIMR.

    ## evaluation matric was utilized from Chen te al.: https://github.com/JiawenChenn/St-review/blob/main/scripts/evaluation_metric.R
    library(energy)
    library(data.table)
    library(dplyr)
    library(funfun)
    evaluation_metric = function(truth,data){
      result_all=list()
      intersect_spot=intersect(rownames(truth),rownames(data))
      truth=truth[intersect_spot,]
      data=data[intersect_spot,]

      data=data.frame(data,check.names = F)
      truth=data.frame(truth,check.names=F)
      
      intersect_celltype=intersect(colnames(truth),colnames(data))
      truth=truth[,intersect_celltype]
      data=data[,intersect_celltype]
      
      mse=calc.mse(t(truth), t(data), rsq = FALSE)
      dcor=rep(0,ncol(data))
      diff_all = c()
      for(i in 1:ncol(data)){
        dcor[i]=as.numeric(dcor.test(truth[,i], data[,i], R=200)$statistic)
        diff_all=rbind(diff_all,data.frame(V1=rownames(truth),value=data[,i]-truth[,i],cell_type=colnames(data)[i]))
      }
      result_all[['mse']]=data.table(V1=rownames(truth),value=mse)
      result_all[['rmse']]=data.table(V1=rownames(truth),value=sqrt(mse))
      result_all[['dcor']]=data.table(cell_type=colnames(data),value=dcor)
      result_all[['diff']]=diff_all
      return(result_all)
    }

    trueDat = as.data.frame(data$trueProp)
    rownames(trueDat) = colnames(spatial_count)
    trueDat = sweep(trueDat,1,rowSums(trueDat),'/')
    trueDat$`(15) Fibroblast-like (outflow tract & valve related & coronary & mediastinal vasculature related)` = 
      trueDat$`(8) Fibroblast-like (outflow tract & valve related)`+
      trueDat$`(4) Fribroblast-like (coronary & mediastinal vasculature related)`
    trueDat = trueDat[,match(colnames(proportion),colnames(trueDat))]
    trueDat = trueDat[match(rownames(proportion),rownames(trueDat)),]

    SWIMR_compare = evaluation_metric(trueDat,proportion)
    print(summary(SWIMR_compare$mse$value))

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## 6.630e-07 1.960e-03 4.132e-03 6.807e-03 7.791e-03 1.098e-01

    ## CARD.visualize.pie was from CARD
    library(CARD)
    CARD.visualize.pie(
      proportion = weight,
      spatial_location = spatial_location[rownames(trueDat),], 
      colors = c("#994636","#DBD8AE"), 
      radius = 250)+ggtitle('Weight distribution')+
      scale_y_reverse() +
      theme(plot.title = element_text(size=20, face='bold',hjust = 0.5),
            legend.text = element_text(size = 12),
            legend.title = element_blank())+
      guides(fill = guide_legend(title = "Reference"))

![](/net/mulan/home/pyzhao/Project/Project1/code/SWIMR/vignettes/SWIMR_Human_developing_heart_files/figure-markdown_strict/my_plot-1.png)

    CARD.visualize.pie(
      proportion = proportion,
      spatial_location = spatial_location[rownames(proportion),], 
      colors = NULL, 
      radius = 220)+ggtitle('Cell type Proportion')+
      scale_y_reverse() +
      theme(plot.title = element_text(size=20, face='bold',hjust = 0.5),
            legend.text = element_text(size = 8),
            legend.title = element_blank())+
      guides(fill = guide_legend(nrow = 6))

![](/net/mulan/home/pyzhao/Project/Project1/code/SWIMR/vignettes/SWIMR_Human_developing_heart_files/figure-markdown_strict/my_plot-2.png)
