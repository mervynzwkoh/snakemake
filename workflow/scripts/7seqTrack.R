library(argparse)

parser <- ArgumentParser("This script generates a SeqTrack transmission nodal network based on the hellingerMatrix.Rda and the accessionList.xlsx file")
parser$add_argument("source", type = "character", help = "The source directory containing the seqtrack directory")
args <- parser$parse_args()

src <- args$source
num_partitions <- length(list.dirs(paste0(src, "/seqtrack"), recursive = FALSE))

# build_tree function not used
build_tree <- function(masterlist, csvfromvcfdir, samplenumber, outputloc) {
  library(adegenet)
  library(ape)
  library(stats)
  library(graphics)
  # load data matrix
  # load("/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/AY.122/AY.3/pt0/hellingerMatrix.Rda")
  load("/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/testDataset/lofreq/tinydemo/modification/hellingerMatrix.Rda")
  d <- as.dist(Dmat)
  tree <- hclust(d)
  # plot(tree, labels=rownames(d), col=tree$height, cex=0.5, hang=-1)
  hcd <- as.dendrogram(tree)
  # Default plot
  plot(hcd, type = "rectangle", ylab = "Height")
  return(NULL)
}

build_network <- function(masterlist, csvfromvcfdir, outputloc, x) {
  print("hi")
  library(readxl)

  data <- readxl::read_excel(masterlist)
  if (nrow(data) == 1) {
    return()
  }

  library(adegenet)
  library(ape)
  library(data.table)
  library(visNetwork)
  library(stringr)
  library(webshot)
  library(tidyverse)
  library(hash)
  library(igraph) # current library for graphing

  # Master List #To do: set sample no=number of csv file
  ref <- read_excel(masterlist)
  dates <- as.Date(ref$ReleaseDate) # store release dates
  # load("/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/AY.122/AY.3/pt0/hellingerMatrix.Rda")
  load(paste0(outputloc, "hellingerMatrix.Rda"))
  print(x) # status
  Dmat <- as.matrix(Dmat)
  # JUST TESTING TO ADD HELLINGER DISTANCES
  # Dmat[3,][9] = Dmat[3,][9] + 0.00243
  # Dmat[9,][3] = Dmat[9,][3] + 0.00243
  print(length(Dmat[1, ]))
  consensus <- seqTrack(Dmat, x.names = ref$Accession, x.dates = dates) # to fill in labels
  print(consensus)
  if (length(unique(consensus$date)) == 1) {
    print("breaking on...")
    print(x)
    return(NULL)
  }
  # consensustre <- plot(consensus, main = paste0('SeqTrack reconstruction',ref$Accession), col.pal=funky) #this works
  # newwd <- paste0(outputloc, "seqtrackDated.pdf")
  # print(newwd)
  # pdf(newwd)
  # Function to increase node separation (for explanatory details, see the link below)
  # Source: http://stackoverflow.com/a/28722680/496488
  # edge.color: bright red (more mut) -> dark red (less mut)
  # plot(consensus, main = paste0("SeqTrack Reconstruction of partition ", x), edge.color = "black", edge.arrow.size = 0.5, vertex.shape = "none", vertex.size = 1, vertex.label = ref$Accession, edge.label = NA, vertex.label.cex = 0.6, arrow.size = 0.3, margin = -0.08)

  # visNetwork nodal network
  samplenum <- length(Dmat[1, ])
  nodes <- data.frame(id = 1:samplenum, label = ref$Accession, shape = "circle", font.size = 20, color = "lightgrey")
  # nodes$color[52] <- "red"
  edges <- data.frame(from = consensus$id, to = consensus$ances, value = consensus$weight)
  network <- visNetwork(nodes, edges, main = paste0("SeqTrack Reconstruction of Partition", x)) %>%
    # the bottom part can omit
    visLegend() %>% visExport(
      type = "pdf", name = "export-network",
      float = "left", label = "Save network"
    )

  # vector for strength of each edge but dk how to implement in the form of thickness of the edge
  darkness <- c()
  for (i in 1:samplenum) {
    id <- ref$Accession[i]
    ances <- consensus$ances[i]
    if (is.na(ances)) {
      darkness <- c(darkness, NA)
      next
    }
    ancessaccession <- ref$Accession[ances]
    darkness <- c(darkness, Dmat[id, ancessaccession])
  }
  network <- network %>%
    visEdges(arrows = "from", width = "value")

  visSave(network, file = paste0("partition", x, ".html"))
  # dev.off()
  # consensustre2 <- plot(consensus, main = paste0('SeqTrack reconstruction'), edge.color='black', vertex.label=ref$Accession)
  return() # CHANGED
  # scft = max(dstmtx[FTM1 < FTM2]$helldst) * min(dstmtx[FTM1 < FTM2]$daydiff) # scale factor
  # scft = max(dstmtx[FTM1 < FTM2]$helldst) / max(dstmtx[FTM1 < FTM2]$spacediff / dstmtx[FTM1 < FTM2]$daydiff) # scale factor

  # code below is not used by me... (visnetwork)
  dscft <- c()
  for (h in dstmtx$helldst) {
    print(h) # DEBUG
    print(dscft)
    dscft <- c(dscft, as.integer(strsplit(as.character(h), "e-")[[1]][2])) # error ? hard coding
  }
  dscft <- max(dscft)
  print(dscft)
  mysqtk <- list()
  # barriosens = list()
  sqtkindex <- sort(which(month(ref$ReleaseDate) %in% c(1:12)))
  # dscft<-max(dstmtx$helldst)
  # dscft<-max(dstmtx[FTM1 < FTM2]$helldst) / max(dstmtx[FTM1 < FTM2]$spacediff / dstmtx[FTM1 < FTM2]$daydiff) # scale factor

  # seqTrack
  for (i in 1:101)
  {
    Dmtx <- Dmat + seq(0, 1, 0.01)[i] * dscft * D2
    mysqtk[[i]] <- seqTrack(Dmtx, x.names = ref$Accession, x.dates = as.Date(ref$ReleaseDate)[sqtkindex], best = "min")
    mysqtk[[i]] <- as.data.table(mysqtk[[i]])
    mysqtk[[i]][, daydiff := as.numeric(date - ances.date)]
    mysqtk[[i]][, ":="(i1 = sqtkindex[ances], i2 = sqtkindex[id])]
    cgmrt <- ref
    # for when geographical data is given:
    '
  mysqtk[[i]][, ":="(`House ID1` = cgmrt$`House ID`[i1], `House ID2` = cgmrt$`House ID`[i2],
                     `Participant ID1` = cgmrt$`Participant ID`[i1], `Participant ID2` = cgmrt$`Participant ID`[i2],
                     BARRIO1 = cgmrt$BARRIO[i1], BARRIO2 = cgmrt$BARRIO[i2],
                     LATITUD1 = cgmrt$LATITUD[i1], LATITUD2 = cgmrt$LATITUD[i2],
                     LONGITUD1 = cgmrt$LONGITUD[i1], LONGITUD2 = cgmrt$LONGITUD[i2],
                     FTM1 = cgmrt$FTM[i1], FTM2 = cgmrt$FTM[i2])]
  mysqtk[[i]] = mysqtk[[i]][!is.na(daydiff), .(i1, i2, `House ID1`, `House ID2`, `Participant ID1`, `Participant ID2`,
                                               `BARRIO1`, `BARRIO2`, `LONGITUD1`, `LONGITUD2`, `LATITUD1`, `LATITUD2`,
                                               FTM1, FTM2, daydiff, weight)]

  mysqtk[[i]][, i1i2 := paste(i1, i2, sep = "_")]
  barriosens[[i]] = mysqtk[[i]][daydiff %in% 4:28, .N, by = .(BARRIO1, BARRIO2)][barriodt0,
                                                                                 on = c(BARRIO1 = "BARRIO1", BARRIO2 = "BARRIO2")]
  barriosens[[i]][is.na(N), N := 0]

  barriosens[[i]][, ":="(Nsum = sum(N)), by = .(BARRIO2)]
  barriosens[[i]][, ":="(p = N/Nsum)]
                     '
  }

  'barriosensdt = do.call("rbind", lapply(1:101, function(i) barriosens[[i]][, .(BARRIO1, BARRIO2, N)]))
barriosensdt = barriosensdt[, .(N = sum(N)), by = .(BARRIO1, BARRIO2)]
barriosensdt[, ":="(Nsum = sum(N)), by = .(BARRIO2)]
barriosensdt[, ":="(p = N/Nsum)]

mysqtkdt = do.call("rbind", lapply(1:101, function(i) mysqtk[[i]]))
mysqtk[[1]][, robust := sapply(i1i2, function(x) nrow(mysqtkdt[i1i2 == x])/101)]'

  # SeqTrack:
  lofrq <- seqTrack(Dmat, x.names = modified, x.dates = dates, best = "min")
  print(lofrq)
  nodes <- data.frame(
    id = 1:samplenumber,
    label = unlist(modified), font.size = 20, # changed size
    shape = "circle",
    color = "#87CEEB", # changed color
    shadow = FALSE
  )
  edges <- data.frame(
    from = dstmtx$i1, to = dstmtx$i2, length = 800, # changed to default length
    width = unlist(dstmtx$helldst) # set width to correspond to hellinger distance
    , arrows = "to",
    dashes = FALSE,
    smooth = FALSE,
    shadow = FALSE
  )

  y <- edges[order(edges$width), ] # order hellinger distance in increasing order
  x <- edges[order(edges$width), ] # order hellinger distance in increasing order
  # we use increasing order since the larger the hellinger distance the smaller the probability
  # of transmission thus, the smaller widths would be preferred.
  store <- hash() # hash contains key value pairs. key: node values
  # Value: nodes accessible by the node in the key

  findpath <- function(start, end, ht) {
    # check if there is already a path from start node to end node
    # output true if there is already a path from start to end and false otherwise
    n <- start
    flag <- FALSE
    n <- start
    visitedkey <- c()
    while (!(end %in% n) & length(intersect(start, hash::values(ht))) >= 1) {
      ni <- c()
      for (k in hash::keys(ht)) {
        visitedkey <- c(visitedkey, k)
        if (length(intersect(ht[[k]], n)) >= 1) {
          ni <- c(ni, k)
          print(ni)
        }
      }
      n <- ni
      print("here") # DEBUG
      if (length(n) == 0) {
        break
      }
    }
    if (end %in% n) {
      flag <- TRUE
    }
    if (start %in% hash::keys(ht)) {
      if (end %in% ht[[toString(start)]]) {
        flag <- TRUE
      }
    }

    return(flag)
  }

  visited <- c()
  toremove <- c()
  for (r in 1:nrow(y)) {
    # loop through each row to find which rows are irrelevant for final tree
    # this ensures that each node only has one arrow pointing to it so that the final
    # output is a hierarchical tree

    front <- y[r, ]$from
    taill <- y[r, ]$to

    if (length(unique(visited)) == samplenumber) {
      # we break the loop when all the samples have already been accounted for
      break
    } else {
      # we ignore if there is already a path from front to taill or there is already an arrow pointing to taill
      if (findpath(front, taill, store) | (taill %in% visited)) {
        toremove <- c(toremove, r)
      } else {
        # if not we need to join front to taill and update the visited values and store the
        # fact that taill can be accessed by front
        visited <- c(visited, taill)
        print(paste0("visited:", visited))
        if (has.key(toString(front), store)) {
          store[[toString(front)]] <- c(store[[toString(front)]], taill)
        } else {
          store[[toString(front)]] <- c()
          store[[toString(front)]] <- c(store[[toString(front)]], taill)
        }
      }
    }
  }
  print(toremove) # DEBUG
  print(x)

  if (!is.null(toremove)) { # added check
    x <- x[-toremove, ]
  }

  # remove irrelevant rows and we should be left with the same number of rows as the number of samples
  # but for 2 samples, only 1 object / 1 arrow
  newwidths <- seq(from = 10, to = 10 / length(rownames(x)), by = -10 / length(rownames(x)))
  # create new widths for the final network, the larger the width, the higher the probability of trasmission
  x$width <- newwidths
  network <<- visNetwork(nodes, x, height = "1000px", width = "100%", main = list(text = "LoFreq", style = "font-family:Comic Sans MS;color:blue;font-size:25px;text-align:left;"))
  network <<- network %>%
    visHierarchicalLayout(direction = "LR") %>%
    visPhysics(enabled = FALSE)
  # store viNetwork in html file
  fname <- paste0("network", ".html")
  visSave(network, fname)
  # save html file as png file:
  webshot(fname,
    delay = 0.5, zoom = 2, file = paste0("faster", ".png"),
    vwidth = 900, vheight = 900
  )
  # move file to the outputloc folder:
  file.copy(paste0(getwd(), "/", "faster.png"), outputloc, overwrite = TRUE)
  file.remove(paste0(getwd(), "/", "faster.png"))
  return(dstmtx)
}

finalfn <- function(inputfiledir, vcfinputfilename, fastainputfilename, outputfiledir, x) {
  masterexcelfile <- paste0(inputfiledir, list.files(inputfiledir, ".xlsx"))
  build_network(masterexcelfile, outputfiledir, outputfiledir, x)
}

# 32 samples: 1hr 21 min
# 3 samples: 7 min

# for (x in 0:3) {
#  #files named as such: ERR1234567.lofreq.fasta, ERR1234567.lofreq.vcf
#  dr = paste0("/Users/windsorkoh/OneDrive - National University of Singapore/Y2 UROPs/MDS/partitions/pt",x,"/")
#  finalfn(dr,"vcf","fasta",dr,x)
# }

for (x in 0:(num_partitions - 1)) {
  setwd(paste0(src, "/seqtrack/pt", x))
  dr <- paste0(src, "/seqtrack/pt", x, "/")
  finalfn(dr, "vcf", "fasta", dr, x)
}
