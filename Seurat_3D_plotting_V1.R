# 3D tSNE plotting of scRNAseq Data
# The following is a length of code generated to create nice 
# 3D tSNE plots of seurat v3.0.0 objects utilizing the visualization 
# scatterplot3d

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code :)

# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. Fr more information please see: https://satijalab.org/seurat/

# Contributors (by their Github handles):
# @Dragonmasterx87 (Dept. of Cell Biology, UM)

#Install packages and dependencies
install.packages('scatterplot3d')
install.packages('rgl')
install.packages('rmarkdown')
install.packages('magick')

# Load packages
library(magick)
library(rmarkdown)

yourseuratobject <- IAmASeruatObjectCreatedWith.Seurat3.0.0.9150

yourseuratobject <- RunTSNE(yourseuratobject,
                        reduction.use = "pca",
                        dims.use = 1:10,
                        dim.embed = 3)

tsne_1 <- yourseuratobject[["tsne"]]@cell.embeddings[,1]
tsne_2 <- yourseuratobject[["tsne"]]@cell.embeddings[,2]
tsne_3 <- yourseuratobject[["tsne"]]@cell.embeddings[,3]

# If you are using Seurat v2.3.4 the following code is what you need to extract information for cell embeddings
# tsne_2 <- yourseuratobject@reductions$tsne@cell.embeddings[,2]
# tsne_3 <- yourseuratobject@reductions$tsne@cell.embeddings[,3]

# if you get errors in the color scheme it’s because the x:y numbers are off, make sure y = number of clusters
# It’s nice to look at gene expression of the ‘orientation’ 2D tSNE so that you can identify clusters in the 3D map
# If you have run a analysis and know your color combinations use the same ones in the same order, they will auto-correspond to correct clusters (check to be sure)
# The following example is for a Seurat object which has 21 clusters (0-20)

TSNEPlot(yourseuratobject, label = FALSE, 
         cols = c("lightseagreen",
                  "gray50",
                  "darkgreen",
                  "red4",
                  "red",
                  "turquoise4",
                  "black",
                  "yellow4",
                  "royalblue1",
                  "lightcyan3",
                  "peachpuff3",
                  "khaki3",
                  "gray20",
                  "orange2",
                  "royalblue4",
                  "yellow3",
                  "gray80",
                  "darkorchid1",
                  "lawngreen",
                  "plum2",
                  "darkmagenta"),
         pt.size = 2)

FeaturePlot(object = yourseuratobject, features = c("PECAM1"), min.cutoff =0, max.cutoff = 1, label = FALSE, 
            cols = c("grey", "red"), pt.size = 2)

# 3D plotting
# Note how the color combinations remain the same

library(rgl) #interactive 3d plotting
plot3d(x = tsne_1, y = tsne_2, z = tsne_3,
       col = c("lightseagreen",
               "gray50",
               "darkgreen",
               "red4",
               "red",
               "turquoise4",
               "black",
               "yellow4",
               "royalblue1",
               "lightcyan3",
               "peachpuff3",
               "khaki3",
               "gray20",
               "orange2",
               "royalblue4",
               "yellow3",
               "gray80",
               "darkorchid1",
               "lawngreen",
               "plum2",
               "darkmagenta")[yourseuratobject@active.ident],
       type = "s", 
       size = 0.5, 
       box = FALSE)

# Run plot3d and while the rgl widget is open run the code below to generate a html file in the plots panel of RStudio
# You can save as a html file via export. These files work with optimal resolution and are most user friendly in Google Chrome 

rgl::rglwidget() #save as html
