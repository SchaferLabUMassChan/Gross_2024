
library(Seurat)

library(dplyr)

library(Matrix)

library(reticulate)

library(scCustomize)

library(ggplot2)


########data#######
#read in data, transpose, remove blanks, remove volume <50, add metadata and normalize for volume

SC3M <- read.table(file='counts_and_metadata_3M.csv', sep=",", header=TRUE,row.names=1)
SC3M <- t(SC3M)
SC3M.Genes <- SC3M[1:400,]
SC3M[467,] <- colSums(SC3M.Genes)
Real.Cells <- which(SC3M[467,] > 10)
SC3M <- SC3M[,Real.Cells]
Big.Cells <- which(SC3M[464,]>50)
SC3M <- SC3M[,Big.Cells]
volume <- SC3M[464,]
head(volume)
center.x <- SC3M[465,]
head(center.x)
FOV <- SC3M[463,]
head(FOV)
center.y <- SC3M[466,]
head(center.y)
mean.RNA <-mean(SC3M[467,])
head(mean.RNA)
SC3M <- SC3M[1:400,]
SC3M <- SC3M/mean.RNA
SC3M <- SC3M/volume
SC3M_s <-CreateSeuratObject(counts=SC3M)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
SC3M_s@meta.data <- cbind(SC3M_s@meta.data,Volumes)
SC3M_s@meta.data <- cbind(SC3M_s@meta.data,Center.x)
SC3M_s@meta.data <- cbind(SC3M_s@meta.data,Center.y)



#create metadata for conditions
SC3M_s$sample<- "3M"

SC18M <- read.table(file='counts_and_metadata_3M.csv', sep=",", header=TRUE,row.names=1)
SC18M <- t(SC18M)
SC18M.Genes <- SC18M[1:400,]
SC18M[467,] <- colSums(SC18M.Genes)
Real.Cells <- which(SC18M[467,] > 10)
SC18M <- SC18M[,Real.Cells]
Big.Cells <- which(SC18M[464,]>50)
SC18M <- SC18M[,Big.Cells]
volume <- SC18M[464,]
head(volume)
center.x <- SC18M[465,]
head(center.x)
FOV <- SC18M[463,]
head(FOV)
center.y <- SC18M[466,]
head(center.y)
mean.RNA <-mean(SC18M[467,])
head(mean.RNA)
SC18M <- SC18M[1:400,]
SC18M <- SC18M/mean.RNA
SC18M <- SC18M/volume
SC18M_s <-CreateSeuratObject(counts=SC18M)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
SC18M_s@meta.data <- cbind(SC18M_s@meta.data,Volumes)
SC18M_s@meta.data <- cbind(SC18M_s@meta.data,Center.x)
SC18M_s@meta.data <- cbind(SC18M_s@meta.data,Center.y)



#create metadata for conditions
SC18M_s$sample<- "18M"


#######main clustering######

#merge into single object 
all.combined <- merge(SC3M_s, y = c(SC18M_s), add.cell.ids = c("1", "2"))


#######lesion vs no lesion#######
library(sp)

##### sample1  ######
#creates the subset, get 1 sample/region from the combined object
SC3M <- subset(x = all.combined, subset = (sample =="3M"))
#gets the x coordinate for every cell
SC3M.x <- SC3M@meta.data[,5]
head(SC3M.x)
#gets the y coordinate for every cell
SC3M.y <- SC3M@meta.data[,6]
head(SC3M.y)

#coordinates for the regions of interest
SC3M_SC_lesion_1_x <- c(5322,
                         5404,
                         5491,
                         5535,
                         5599,
                         5578,
                         5548,
                         5437,
                         5424,
                         5442,
                         5303,
                         5202,
                         5178,
                         5223,
                         5322
)
SC3M_SC_lesion_1_y <- c(3795,
                         3788,
                         3707,
                         3548,
                         3315,
                         3247,
                         3241,
                         3390,
                         3497,
                         3130,
                         3029,
                         3372,
                         3644,
                         3804,
                         3795
)

SC3M_SC_nonlesion_1_x <- c(4460,
                            4468,
                            4908,
                            5110,
                            5034,
                            4713,
                            4424,
                            4449
)
SC3M_SC_nonlesion_1_y <- c(3085,
                            3107,
                            3160,
                            3027,
                            2931,
                            2876,
                            2977,
                            3120
)


SC3M_SC_lesion_2_x <- c(9662,
                         9845,
                         9939,
                         9906,
                         9803,
                         9639,
                         9615,
                         9595,
                         9652
)
SC3M_SC_lesion_2_y <- c(4561,
                         4470,
                         4241,
                         3970,
                         3877,
                         4156,
                         4337,
                         4368,
                         4561
)

SC3M_SC_nonlesion_2_x <- c(10536,
                            10481,
                            10653,
                            10836,
                            10862,
                            10658,
                            10543,
                            10536
                            
)
SC3M_SC_nonlesion_2_y <- c(4409,
                            4548,
                            4980,
                            5003,
                            4867,
                            4496,
                            4406,
                            4409
                            
)


# find cells in your dataset that are inside your area of interest using point in polygon
SC3M_pol_lesion_1 <- point.in.polygon(SC3M.x, SC3M.y, SC3M_SC_lesion_1_x, SC3M_SC_lesion_1_y)
SC3M_pol_nonlesion_1 <- point.in.polygon(SC3M.x, SC3M.y, SC3M_SC_nonlesion_1_x, SC3M_SC_nonlesion_1_y)
SC3M_pol_lesion_2 <- point.in.polygon(SC3M.x, SC3M.y, SC3M_SC_lesion_2_x, SC3M_SC_lesion_2_y)
SC3M_pol_nonlesion_2 <- point.in.polygon(SC3M.x, SC3M.y, SC3M_SC_nonlesion_2_x, SC3M_SC_nonlesion_2_y)


SC3M.info <- cbind.data.frame(SC3M_pol_lesion_1, SC3M_pol_nonlesion_1, SC3M_pol_lesion_2, SC3M_pol_nonlesion_2)
head(SC3M.info, 50)

SC3M.info <- data.frame(SC3M.info)
SC3M.info
SC3M@meta.data <- cbind(SC3M@meta.data,SC3M.info)

SC3M_lesion_1 <- subset(SC3M,  SC3M_pol_lesion_1 == 1)
head(SC3M_lesion_1, 100)
SC3M_lesion_2 <- subset(SC3M,  SC3M_pol_lesion_2 == 1)
head(SC3M_lesion_1, 100)
SC3M_nonlesion_1 <- subset(SC3M,  SC3M_pol_nonlesion_1 == 1)
head(SC3M_nonlesion_1, 100)
SC3M_nonlesion_2 <- subset(SC3M,  SC3M_pol_nonlesion_2 == 1)
head(SC3M_nonlesion_2, 100)


SC3M_lesion_1$lesion <- "Lesion"
SC3M_lesion_2$lesion <- "Lesion"
SC3M_nonlesion_1$lesion <- "Non lesion"
SC3M_nonlesion_2$lesion <- "Non lesion"

SC3M_lesion_1$number <- "1_3M"
SC3M_lesion_2$number <- "2_3M"
SC3M_nonlesion_1$number <- "1_3M"
SC3M_nonlesion_2$number <- "2_3M"

SC3M_SC <- merge(SC3M_lesion_1, y=c(SC3M_lesion_2,SC3M_nonlesion_1,SC3M_nonlesion_2 ))

#Example of spatial plotting of region of interest
dfall<- data.frame(SC3M_lesion_1$center.x[WhichCells(object = subset(SC3M_lesion_1,subset = sample=="3M"))],SC3M_lesion_1$center.y[WhichCells(object = subset(SC3M_lesion_1,subset = sample=="3M"))],SC3M_lesion_1@active.ident[WhichCells(object = subset(SC3M_lesion_1,subset = sample=="3M"))])
x <- SC3M_lesion_1$center.x[WhichCells(object = subset(SC3M_lesion_1,subset = sample=="3M"))]
y <- SC3M_lesion_1$center.y[WhichCells(object = subset(SC3M_lesion_1,subset = sample=="3M"))]
group <- SC3M_lesion_1@active.ident[WhichCells(object = subset(SC3M_lesion_1,subset = sample=="3M"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=1)+ theme_linedraw()+theme_classic()+coord_fixed()
ggsave('full_xy_plot.png', height=10, width=13)+coord_fixed(ratio=1)


##### sample2 ######
#creates the subset, get 1 sample/region from the combined object
SC18M <- subset(x = all.combined, subset = (sample =="18M"))
#gets the x coordinate for every cell
SC18M.x <- SC18M@meta.data[,5]
head(SC18M.x)
#gets the y coordinate for every cell
SC18M.y <- SC18M@meta.data[,6]
head(SC18M.y)

#coordinates for the regions of interest
SC18M_SC_lesion_1_x <- c(655.511,
                         757,
                         759,
                         759,
                         703,
                         556,
                         486,
                         315,
                         234,
                         228,
                         293,
                         369,
                         398,
                         522,
                         655
                         
)
SC18M_SC_lesion_1_y <- c(2404.84,
                         2312,
                         2254,
                         2114,
                         2030,
                         1840,
                         1788,
                         1702,
                         1755,
                         1916,
                         2067,
                         2211,
                         2255,
                         2352,
                         2404
                         
)

SC18M_SC_nonlesion_1_x <- c(1032,
                            1262,
                            1210,
                            905,
                            779,
                            869,
                            1032
                            
)
SC18M_SC_nonlesion_1_y <- c(3419,
                            3428,
                            3200,
                            2889,
                            3039,
                            3222,
                            3419
                            
)
SC18M_SC_lesion_2_x <- c(3601,
                         3828,
                         3927,
                         4023,
                         4078,
                         4070,
                         4037,
                         3996,
                         3731,
                         3563,
                         3601,
                         3601,
                         3638,
                         3601
                         
                         
)
SC18M_SC_lesion_2_y <- c(4652,
                         4780,
                         4771,
                         4728,
                         4689,
                         4585,
                         4542,
                         4517,
                         4404,
                         4444,
                         4652,
                         4652,
                         4697,
                         4652
                         
)

SC18M_SC_nonlesion_2_x <- c(1421,
                            1528,
                            1753,
                            2020,
                            1865,
                            1620,
                            1446,
                            1421
                            
                            
)
SC18M_SC_nonlesion_2_y <- c(3507,
                            3491,
                            3350,
                            3081,
                            2904,
                            3038,
                            3312,
                            3507
                            
                            
)




# find cells in your dataset that are inside your area of interest using point in polygon
SC18M_pol_lesion_1 <- point.in.polygon(SC18M.x, SC18M.y, SC18M_SC_lesion_1_x, SC18M_SC_lesion_1_y)
SC18M_pol_nonlesion_1 <- point.in.polygon(SC18M.x, SC18M.y, SC18M_SC_nonlesion_1_x, SC18M_SC_nonlesion_1_y)
SC18M_pol_lesion_2 <- point.in.polygon(SC18M.x, SC18M.y, SC18M_SC_lesion_2_x, SC18M_SC_lesion_2_y)
SC18M_pol_nonlesion_2 <- point.in.polygon(SC18M.x, SC18M.y, SC18M_SC_nonlesion_2_x, SC18M_SC_nonlesion_2_y)


SC18M.info <- cbind.data.frame(SC18M_pol_lesion_1, SC18M_pol_nonlesion_1, SC18M_pol_lesion_2,SC18M_pol_nonlesion_2)
head(SC18M.info, 50)


SC18M.info <- data.frame(SC18M.info)
SC18M.info
SC18M@meta.data <- cbind(SC18M@meta.data,SC18M.info)

SC18M_lesion_1 <- subset(SC18M,  SC18M_pol_lesion_1 == 1)
head(SC18M_lesion_1, 100)

SC18M_nonlesion_1 <- subset(SC18M,  SC18M_pol_nonlesion_1 == 1)
head(SC18M_nonlesion_1, 100)

SC18M_lesion_1 <- subset(SC18M,  SC18M_pol_lesion_1 == 1)
head(SC18M_lesion_1, 100)
SC18M_lesion_2 <- subset(SC18M,  SC18M_pol_lesion_2 == 1)
head(SC18M_lesion_1, 100)
SC18M_nonlesion_1 <- subset(SC18M,  SC18M_pol_nonlesion_1 == 1)
head(SC18M_nonlesion_1, 100)
SC18M_nonlesion_2 <- subset(SC18M,  SC18M_pol_nonlesion_2 == 1)
head(SC18M_nonlesion_2, 100)


SC18M_lesion_1$lesion <- "Lesion"
SC18M_lesion_2$lesion <- "Lesion"
SC18M_nonlesion_1$lesion <- "Non lesion"
SC18M_nonlesion_2$lesion <- "Non lesion"

SC18M_lesion_1$number <- "1_18M"
SC18M_lesion_2$number <- "2_18M"
SC18M_nonlesion_1$number <- "1_18M"
SC18M_nonlesion_2$number <- "2_18M"

SC18M_SC <- merge(SC18M_lesion_1, y=c(SC18M_lesion_2,SC18M_nonlesion_1,SC18M_nonlesion_2 ))

dfall<- data.frame(SC18M_lesion_1$center.x[WhichCells(object = subset(SC18M_lesion_1,subset = sample=="18M"))],SC18M_lesion_1$center.y[WhichCells(object = subset(SC18M_lesion_1,subset = sample=="18M"))],SC18M_lesion_1@active.ident[WhichCells(object = subset(SC18M_lesion_1,subset = sample=="18M"))])
x <- SC18M_lesion_1$center.x[WhichCells(object = subset(SC18M_lesion_1,subset = sample=="18M"))]
y <- SC18M_lesion_1$center.y[WhichCells(object = subset(SC18M_lesion_1,subset = sample=="18M"))]
group <- SC18M_lesion_1@active.ident[WhichCells(object = subset(SC18M_lesion_1,subset = sample=="18M"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=1)+ theme_linedraw()+theme_classic()+coord_fixed()
ggsave('full_xy_plot.png', height=10, width=13)+coord_fixed(ratio=1)

#merge all regions of interest
all.combined <- merge(SC3M_SC, y=c(SC18M_SC))

Idents(SC3M_SC) <- SC3M_SC@meta.data$lesion
levels(SC3M_SC)

#spatial plot of lesion & non lesion 
dfall<- data.frame(SC3M_SC$center.x[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))],SC3M_SC$center.y[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))],SC3M_SC@active.ident[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))])
x <- SC3M_SC$center.x[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))]
y <- SC3M_SC$center.y[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))]
group <- SC3M_SC@active.ident[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=1)+ theme_linedraw()+theme_classic()+coord_fixed()

Idents(SC3M_SC) <- SC3M_SC@meta.data$number
levels(SC3M_SC)

#spatial plot of the different samples within one condition
dfall<- data.frame(SC3M_SC$center.x[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))],SC3M_SC$center.y[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))],SC3M_SC@active.ident[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))])
x <- SC3M_SC$center.x[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))]
y <- SC3M_SC$center.y[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))]
group <- SC3M_SC@active.ident[WhichCells(object = subset(SC3M_SC,subset = sample=="3M"))]
ggplot(dfall,aes(x,y,colour=group))+geom_point(size=1)+ theme_linedraw()+theme_classic()+coord_fixed()

Idents(all.combined) <- all.combined@meta.data$lesion
levels(all.combined)

#QC checks 
VlnPlot(all.combined, features=c('nFeature_RNA', 'nCount_RNA'), ncol=3, pt.size=0)

#remove non-cells
all.combined <- subset(all.combined,nFeature_RNA > 10)

VlnPlot(all.combined,features="nCount_RNA",group.by="sample",pt.size=0)
VlnPlot(all.combined,features="nCount_RNA",group.by="number",pt.size=0)

#normalize data
all.combined <- NormalizeData(all.combined)

#find variable features (reduce from default 2000 to 400 since 400 genes in the panel)

all.combined <- FindVariableFeatures(all.combined, selection.method = "vst", nfeatures = 400)
all.genes <- rownames(all.combined)
all.genes
all.combined <- ScaleData(all.combined, features = all.genes)
all.combined <- RunPCA(all.combined, features = VariableFeatures(object = all.combined),npcs = 50)

# test of number of significant Principal Components
all.combined <- JackStraw(all.combined, num.replicate = 100,dims=40)
all.combined <- ScoreJackStraw(all.combined, dims = 1:40)
JackStrawPlot(all.combined, dims = 1:40)
ElbowPlot(object = all.combined, ndims = 40)

#cluster
all.combined <- FindNeighbors(all.combined, dims = 1:20)
all.combined <- FindClusters(all.combined, resolution = 1.45)
all.combined <- RunUMAP(all.combined, dims = 1:20)
DimPlot(all.combined, reduction = "umap")
DimPlot(all.combined, reduction = "umap", split.by = "lesion")

Idents(all.combined) <- all.combined@meta.data$CellType
levels(all.combined)
Idents(all.combined) <- all.combined@meta.data$seurat_clusters
levels(all.combined)
new.cluster.ids <- c('Microglia', 
                     'OLs', 
                     'Macrophages', 
                     'OLs', 
                     'OLs', 
                     'Astrocytes',
                     'Microglia', 
                     'Endothelial Cells',
                     'Astrocytes', 
                     'OPCs')
names(new.cluster.ids) <- levels(all.combined)
all.combined <- RenameIdents(all.combined, new.cluster.ids)
all.combined$CellType <- Idents(all.combined)
table(Idents(all.combined))

levels(x = all.combined) <- c("Microglia", "Macrophages", "Astrocytes", "OLs", "OPCs", "Endothelial Cells")


markers <- c("P2ry12", "Cx3cr1",  "Lyz2","Aldh1l1", "Aqp4","Plp1", "Sox10", "Cspg4",   'Pdgfra', "Pecam1", "Cldn5" ) 
DotPlot(all.combined, features = markers) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  scale_colour_gradient2(low = "green", mid = "white", high = "darkred")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DimPlot(all.combined, reduction = "umap",  label = T, split.by = "lesion")


Idents(all.combined) <- all.combined@meta.data$CellType
levels(all.combined)

prop.table(table(Idents(all.combined)))
proportion <- prop.table(table(Idents(all.combined), all.combined$lesion))
proportion.df <- data.frame(proportion)
proportion.df

proportion.plot <- ggplot(proportion.df, aes(fill = Var1,
                                             y = Freq, x = Var2))+
  geom_bar(position = "fill", stat = "identity")+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_flip()
proportion.plot + theme_classic()

Idents(all.combined) <- all.combined@meta.data$lesion
levels(all.combined)

table(Idents(all.combined))

senmayo <- list(c('Axl',
                  'C3',
                  'Ccl20',
                  'Ccl3',
                  'Ccl8',
                  'Cd9',
                  'Csf1',
                  'Csf2',
                  'Ctsb',
                  'Cxcl1',
                  'Cxcl10',
                  'Cxcl2',
                  'Cxcl3',
                  'Dkk1',
                  'Edn1',
                  'Esm1',
                  'Fgf1',
                  'Fgf2',
                  'Fgf7',
                  'Gdf15',
                  'Hgf',
                  'Hmgb1',
                  'Igfbp2',
                  'Igfbp3',
                  'Igfbp5',
                  'Il10',
                  'Il13',
                  'Il15',
                  'Il1a',
                  'Il1b',
                  'Il6',
                  'Il7',
                  'Lcp1',
                  'Mmp13',
                  'Mmp12',
                  'Mmp2',
                  'Mmp3',
                  'Mmp9',
                  'Pecam1',
                  'Sema3f',
                  'Serpine1',
                  'Spp1',
                  'Spx',
                  'Timp2',
                  'Tnf',
                  'Vegfc',
                  'Vgf',
                  'Wnt16'
))

all.combined <- AddModuleScore(
  object = all.combined, 
  features = senmayo, 
  ctrl = 5, 
  name = 'senmayo'
)

senmayo <- c(     'Axl',
                  'C3',
                  'Ccl20',
                  'Ccl3',
                  'Ccl8',
                  'Cd9',
                  'Csf1',
                  'Csf2',
                  'Ctsb',
                  'Cxcl1',
                  'Cxcl10',
                  'Cxcl2',
                  'Cxcl3',
                  'Dkk1',
                  'Edn1',
                  'Esm1',
                  'Fgf1',
                  'Fgf2',
                  'Fgf7',
                  'Gdf15',
                  'Hgf',
                  'Hmgb1',
                  'Igfbp2',
                  'Igfbp3',
                  'Igfbp5',
                  'Il10',
                  'Il13',
                  'Il15',
                  'Il1a',
                  'Il1b',
                  'Il6',
                  'Il7',
                  'Lcp1',
                  'Mmp13',
                  'Mmp12',
                  'Mmp2',
                  'Mmp3',
                  'Mmp9',
                  'Pecam1',
                  'Sema3f',
                  'Serpine1',
                  'Spp1',
                  'Spx',
                  'Timp2',
                  'Tnf',
                  'Vegfc',
                  'Vgf',
                  'Wnt16'
)
Idents(all.combined) <- all.combined@meta.data$CellType
levels(all.combined)
DotPlot_scCustom(all.combined, features = senmayo)+coord_flip()


Idents(all.combined) <- all.combined@meta.data$sample
levels(all.combined)
DotPlot_scCustom(all.combined, features = senmayo)+coord_flip()


### subset microglia/macrophages

onlymicromacro <- subset(all.combined,idents = c("Microglia", "Macrophages"))

onlymicromacro <- AddModuleScore(
  object = onlymicromacro, 
  features = senmayo, 
  ctrl = 5, 
  name = 'senmayo_Mg'
)





