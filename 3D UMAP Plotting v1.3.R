# Interacive multimodal 3D UMAP plotting of scRNA sequencing datasets
# The following is a length of code generated to create nice 
# 3D UMAP plots of seurat v3.0.0-v3.1.1 objects utilizing the visualization 
# package plot_ly

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code :)

# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/

# Contributors (by their Github handles):
# @Dragonmasterx87 (Dept. of Cell Biology, UM)
# @msaadsadiq (Dept. of Electrical and Computer Engineering, UM)

# Install plot_ly
install.packages('plotly')

# Load plot_ly
library(plotly)

# Construct a dataframe using data from your pre-clustered Seurat v3.1.1 object
# Here 'seurat_clusters' is list of numeric cluster identities, you can find it here: yourseuratobject[["seurat_cluster"]], 
# or yourseuratobject$seurat_clusters, where 'yourseuratobject' is a Seurat object created with Seurat v3.1.1 (works for v3.0.0 as well)
yourseuratobject <- ThisIsWhateverYourSeuratObjectIsEvenIfItsIntegrated

# Re-run UMAPs that you have accurate calculations for all UMAP(s)
yourseuratobject <- RunUMAP(yourseuratobject,
                            dims = 1:10,
                            n.components = 3L)

# Extract tSNE information from Seurat Object
umap_1 <- yourseuratobject[["umap"]]@cell.embeddings[,1]
umap_2 <- yourseuratobject[["umap"]]@cell.embeddings[,2]
umap_3 <- yourseuratobject[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = yourseuratobject, reduction = "umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~seurat_clusters, 
        colors = c("lightseagreen",
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
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 5, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


# Updates stemming from Issue #9 Having a fixed scale on axes while selecting particular clusters
# @rtoddler thanks for the suggestions!
# Before you plot, set the ranges of the axis you desire. This set axis range will be 
# present across all clusters, and plotly will not adjust for axis length anymore
# this axis length will persist even when selecting some clusters

# xaxis
axx <- list(
  nticks = 4,
  range = c(-10,10) #select range of xaxis
)

# yaxis
axy <- list(
  nticks = 4,
  range = c(-10,10) #select range of yaxis
)

#zaxis
axz <- list(
  nticks = 4,
  range = c(-10,10) #select range of zaxis
)

fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
fig
fig_cube

# Say you wanto make a gene-expression 3D plot, where you can plot gene expression against a color scale
# Here using the same seurat object as above, we extract gene expression information for beta-actin 'ACTB'
# Here we concentrate on SCT normalized data, or log normalized RNA NOT raw counts.
# In addition if you want, you may look at normalised-RNA, SCT or integrated slots, to look at gene expression
# Setting your DefaultAssay() will inform R which assay to pick up expression data from.
DefaultAssay(object = yourseuratobject)
DefaultAssay(object = yourseuratobject) <- "RNA"
DefaultAssay(object = yourseuratobject) <- "integrated"
DefaultAssay(object = yourseuratobject) <- "SCT"

# create a dataframe
plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "ACTB"), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plot.data$changed <- ifelse(test = plot.data$ACTB <1, yes = plot.data$ACTB, no = 1)

# Add the label column, so that now the column has 'cellname-its expression value'
plot.data$label <- paste(rownames(plot.data)," - ", plot.data$ACTB, sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text"
)

# On running this code the HTML output should appear in RStudio. You can save the output as a
# HTML file. Once you have saved, just open the HTML file in any web browser (double click on the html- file
# and if asked select to open with any web browser like google chrome/safari/mozilla/explorer etc).
# It should be have all of the integrated features you saw in the RStudio output file.

########## #
########## #

# Alternative method as designed by @vertesy (Thanks for the suggestions!)
# create a dataframe
goi <- "TOP2A"
plotting.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=goi), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'Expression' column of your dataframe
# Cutoff <- 2
Cutoff <- quantile(plotting.data[,goi], probs = .95)
plotting.data$"ExprCutoff" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plotting.data,
        # name = goi,
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~ExprCutoff, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgrey', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 1), 
        text=~label,
        hoverinfo="text"
) %>%layout(title=goi)

# Thank you for reading and using this code to further your scRNAseq analysis!
# If you liked it, dont forget to acknowledge, fork and star!
# Citation information is within the Readme, please dont forget to cite!
# Have a wonderful day!!
