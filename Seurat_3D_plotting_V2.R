# 3D tSNE plotting of scRNAseq Data
# The following is a length of code generated to create nice 
# 3D tSNE plots of seurat v3.0.0 objects utilizing the visualization 
# package plot_ly

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code :)

# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab please see: https://satijalab.org/seurat/

# Contributors (by their Github handles):
# @Dragonmasterx87 (Dept. of Cell Biology, UM)
# @msaadsadiq (Dept. of Electrical and Computer Engineering, UM)
# @andrewwbutler (Center for Genomics and Systems Biology, NYU)  

# Install plot_ly
install.packages('plotly')

# Load plot_ly
library(plotly)

# Construct a datasframe using data from your pre-clustered Seurat v3.0.0 object
# Here 'seurat_clusters' is list of numeric cluster identities, you can find it here: yourseuratobject[["seurat_cluster"]], 
# or yourseuratobject$seurat_clusters, where 'yourseuratobject' is a Seurat object created with Seurat v3.0.0
plot.data <- FetchData(object = yourseuratobject, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plotting.data, 
        x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
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
        text=~label, #This is that extra column we made earlier for which we will use
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

# Say you wanto make a gene-expression 3D plot, where you can plot gene expression against a color scale
# Here using the same seurat object as above, we extract gene expression information for beta-actin 'ACTB'
# creat a dataframe
plotting.data <- FetchData(object = yourseuratobject, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "ACTB"))

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression 
# will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plotting.data$changed <- ifelse(test = plotting.data$KRT19 <1, yes = plotting.data$KRT19, no = 1)

# Change the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data$ACTB, sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plotting.data, 
        x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
        color = ~changed,
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text"
)

# Thank you for reading and using this code to further your scRNAseq analysis!
# If you liked it, dont forget to acknowledge, fork and star!
# Have a wonderful day!!

