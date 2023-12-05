library(ggtree)
library(ggplot2)
library(treeio)
library(stringr)
library(colorspace)

# treepath <- "/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/FINALISED/datasets/test/T2.raxml.bestTreeCollapsed"
# treepath <- "/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/FINALISED/datasets/trngg/T2.raxml.bestTreeCollapsed"
# treepath <- "/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/FINALISED/datasets/dengue/K1.raxml.bestTree"
input <- read.tree(treepath)
# clpath <- "/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/FINALISED/datasets/test/clusterlabels/"
# clpath <- "/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/FINALISED/datasets/trngg/clusterlabels/"
# clpath <- "/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/FINALISED/datasets/dengue/clusterlabels/"
# clusterlabelpath = "/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/FINALISED/datasets/test/clusterlabels/hier0.000301.csv"
for (measure in c("Leven", "Raw", "JC", "K2P", "Tamura", "Phylo")) {
  if (measure != "Tamura") {
    next
  }
  clusterdirpath <- paste0(clpath, measure, "/")
  files <- list.files(path = clusterdirpath)
  for (f in files) {
    fname <- paste0(clusterdirpath, f)
    # print(fname)
    colourTree(fname)
  }
}
# colourTree("/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/FINALISED/datasets/trngg/clusterlabels/grdtruth2.csv")
colourTree <- function(clusterlabelpath) {
  print(clusterlabelpath)
  clusterss <- read.csv(clusterlabelpath, header = FALSE)

  x <- clusterss[[2]]
  # x <- str_replace_all(x, "[\n]" , "")
  x <- str_replace_all(x, ",", "")
  x <- str_split(x, "\\s+")
  partition <- as.data.frame(x)
  partition <- as.numeric(unlist(x))
  print(max(partition))

  y <- clusterss[[1]] # 1st column
  # y <- str_replace_all(y, "'" , "")
  y <- str_replace_all(y, ",", "")
  # y <- str_replace_all(y, "[\n]" , "")
  y <- str_split(y, "\\s+")
  labs <- as.data.frame(y)
  everything <- rbind(partition, labs)
  everything <- as.data.frame(everything)
  everything <- t(everything) # transpose
  colnames(everything) <- c("partition", "label") # everything$part[3] is 3rd elem's part number and so on
  df <- as.data.frame(everything)
  tree <- full_join(input, df, by = "label")
  colorPalette <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40")
  cPallete <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40", "#895e4e", "#dfae3d", "#e8ded2", "#92efe1", "#575b9c")
  rainbowcolorPalette <- rainbow_hcl(max(partition) - min(partition) + 1, l = 69, c = 69)
  print(max(partition) - min(partition) + 1)
  if (max(partition) - min(partition) + 1 <= 4) {
    plot <- ggtree(tree) + geom_tiplab(aes(colour = factor(partition), data = df$labs), align = TRUE, size = 0.5) + scale_color_manual(values = colorPalette)
  } else {
    plot <- ggtree(tree) + geom_tiplab(aes(colour = factor(partition), data = df$labs), align = TRUE, size = 0.5) + scale_color_manual(values = rainbowcolorPalette)
  }


  treefname <- strsplit(fname, "clusterlabels")[[1]][2]
  treefname <- strsplit(treefname, ".csv")[[1]][1]
  outfilepath <- paste0(dirname(treepath), "/treeplots/tree", treefname, ".png")
  # outfilepath = "/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/FINALISED/datasets/trngg/treeplots/tree/grdtruth.png"
  ggsave(outfilepath, plot, dpi = 700)
}
