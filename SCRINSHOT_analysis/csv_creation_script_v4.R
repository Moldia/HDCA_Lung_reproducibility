
# Install (if it is necessary) and load there libraries
library(dplyr)
library(mefa)
library(Jmisc)
library(stringr)

# give the dimensions of the image (not of the tiles) :
x= 10250
y= 15850


# give number of images (Note: the number should be at least "2", since DAPI-nuclei is the 1st image).
images=33


# set the number of tiles per dimension, according to the settings of the "tiling_with_fiji_v2.ijm" file.
n_x=2
n_y=2


# write the input directory with the folders of the tiled images. Note that the *.csv 
# file with the coordinates will be saved in that folder
b = as.character("C:\\Users\\alexandros.sountoul\\Desktop\\w6_analysis\\input\\")
# set the name of the output file. 
name="hyb1_tiling.csv"


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# run the script
{
setwd(b)

x_size = x/n_x
y_size = y/n_y
tile_n = n_x*n_y

nuclei_path <- paste0(b,"hyb1_c1_ORG\\")


merge <- data.frame(matrix(data = 0, nrow = n_x*n_y,ncol = 3))

for (i in 1:tile_n){
  merge[i,1] <-str_pad(i, 4, pad = "0")
}

for (i in 1:n_y){
  for (j in (0:(n_x-1))){
    for (k in (i + j*n_y)){
      merge[k, 2] <- (i*x_size)-x_size
    }
  }
}

r <- c(1:tile_n)
splitted <- split(r,             
                  cut(seq_along(r),
                      n_x,
                      labels = FALSE))

for (i in 1:n_x){
  for (j in splitted[[i]])
    merge[j, 3] <- (i*y_size)-y_size
}
names(merge) <- c("Metadata_position", "Tile_xPos", "Tile_yPos")

merge <- addCol(merge, value=c(Hyb_step="hyb1"))


merge <- addCol(merge, value=c(Image_PathName_Nuclei=nuclei_path))

merge <- addCol(merge, value=c(Image_FileName_Nuclei="Stack-"))

tiles <- paste0(merge$Image_FileName_Nuclei, merge$Metadata_position, ".tif")

merge$Image_FileName_Nuclei<- tiles

for (i in 2:images) {
  gene_path <- paste0(b, "hyb1_c", i, "_ORG/")
  merge <- addCol(merge, value=c(Image_PathName_gene=gene_path))
  gene_path_name <- paste0("Image_PathName_gene", i-1)
  colnames(merge)[colnames(merge) == "Image_PathName_gene"] <- gene_path_name
  Image_FileName_gene <- merge$Image_FileName_Nuclei
  Image_FileName_gene_name <- paste0("Image_FileName_gene", i-1)
  merge <- cbind(merge, Image_FileName_gene)
  colnames(merge)[colnames(merge) == "Image_FileName_gene"] <- Image_FileName_gene_name
}


write.csv(merge, file=name, row.names=FALSE)

}

