## This script is a stripped-down version of the worksheet for Bowman et al., 2016, ISMEJ.  It should contain
## all essential analyses and plots, though not all plots will appear exactly as they do in the
## publication, except for the classification of flow cytometry events.


## the dataset contains some experimental data that is not relevant to this analysis.
## Read in a file that identifies these for later removal.

treat <- read.table('treatment_sample_data.txt', header = T, row.names = 1)

#### read in pathway abundance  and generate date means ####

paths <- read.csv('pal_timeseries.path_tally.csv', row.names = 1, header = T)
row.names(paths) <- gsub('.', '', row.names(paths), fixed = T)

paths <- paths[which(!row.names(paths) %in% row.names(treat)),]

paths$date <- strptime(sapply(strsplit(as.character(row.names(paths)), '_'), '[[', 3), "%Y%m%d")

path.rownames <- unique(paths$date)
path.matrix <- matrix(ncol = length(colnames(paths)) - 1, nrow = length(path.rownames))

colnames(path.matrix) <- colnames(paths)[1:length(colnames(paths)) - 1]
rownames(path.matrix) <- as.character(path.rownames)

for(date in as.character(paths$date)){
  temp <- paths[which(as.character(paths$date) == date),]
  temp$date <- NULL
  temp <- data.matrix(temp)
  temp_total <- colSums(temp)
  temp_norm <- temp_total / max(temp_total)
  path.matrix[which(rownames(path.matrix) == as.character(date)),] <- temp_norm
}

#### read in edge abundances ####

edge.tally <- read.table('pal_timeseries.edge_tally.csv', sep = ',', header = T, row.names = 1)
edge.data <- as.data.frame(t(read.table('pal_timeseries.edge_data.csv', sep = ',', header = T, row.names = 1)))

row.names(edge.tally) <- gsub('.', '', row.names(edge.tally), fixed = T)
row.names(edge.data) <- gsub('.', '', row.names(edge.data), fixed = T)

## remove experiments and a couple of bad libraries

edge.tally <- edge.tally[which(!row.names(edge.tally) %in% row.names(treat)),]
edge.data <- edge.data[which(!row.names(edge.data) %in% row.names(treat)),]

edge.tally <- edge.tally[which(!row.names(edge.tally) == 'PAL_243_20140121_F02_2844'),]
edge.data <- edge.data[which(!row.names(edge.data) == 'PAL_243_20140121_F02_2844'),]

## convert sample names to dates

edge.tally$date <- strptime(sapply(strsplit(as.character(row.names(edge.tally)), '_'), '[[', 3), "%Y%m%d")
edge.data$date <- strptime(sapply(strsplit(as.character(row.names(edge.data)), '_'), '[[', 3), "%Y%m%d")

## order both dataframes by date

edge.tally <- edge.tally[order(edge.tally$date),]
edge.data <- edge.data[order(edge.data$date),]

## convert all NA to 0

edge.tally[is.na(edge.tally)] <- 0

## get mean abundances for each date

row.names <- unique(edge.tally$date)
col.names <- colnames(edge.tally[,sapply(edge.tally, is.numeric)])

edge.tally.matrix <- matrix(ncol = length(col.names), nrow = length(row.names))
colnames(edge.tally.matrix) <- col.names
rownames(edge.tally.matrix) <- as.character(row.names)

edge.count.matrix <- matrix(ncol = 1, nrow = length(row.names))
rownames(edge.count.matrix) <- as.character(row.names)

edge.sd.matrix <- matrix(ncol = length(col.names), nrow = length(row.names))
colnames(edge.sd.matrix) <- col.names
rownames(edge.sd.matrix) <- as.character(row.names)

edge.sum.matrix <- matrix(ncol = 1, nrow = length(row.names))
rownames(edge.sum.matrix) <- as.character(row.names)

for(date in unique(as.character(edge.tally$date))){
  temp <- edge.tally[which(as.character(edge.tally$date) == date),]
  temp$date <- NULL
  temp$X <- NULL
  temp <- data.matrix(temp)
  temp.total <- mean(rowSums(temp))
  temp.sum <- colSums(temp)
  temp.norm <- temp.sum / sum(temp.sum)
  temp.sd <- apply(temp, 2, sd) / temp.sum
  temp.count <- length(temp[,1])
  edge.tally.matrix[which(rownames(edge.tally.matrix) == as.character(date)),] <- temp.norm
  edge.count.matrix[which(rownames(edge.count.matrix) == as.character(date)),1] <- temp.count
  edge.sum.matrix[which(rownames(edge.sum.matrix) == as.character(date)),1] <- temp.total
  edge.sd.matrix[which(rownames(edge.sd.matrix) == as.character(date)),] <- temp.sd
}

## eliminate edges that have abundance of 0 in all samples

edge.norm <- edge.tally.matrix[,colSums(edge.tally.matrix) > 0]

## get mean data for each date

col.names <- colnames(edge.data)

edge.data.matrix <- matrix(ncol = length(col.names) - 1, nrow = length(row.names))
colnames(edge.data.matrix) <- col.names[which(!col.names == 'date')]
rownames(edge.data.matrix) <- as.character(row.names)

for(date in unique(as.character(edge.data$date))){
  temp <- edge.data[which(as.character(edge.data$date) == date),]
  temp$date <- NULL
  temp <- data.matrix(temp)
  temp.mean <- apply(temp, 2, mean)
  edge.data.matrix[which(rownames(edge.data.matrix) == as.character(date)),] <- temp.mean
}

edge.data.df <- as.data.frame(edge.data.matrix) # will cause you headaches if you leave as matrix
colnames(edge.data.df) <- colnames(edge.data.matrix)
edge.data.df$date <- as.POSIXct(row.names(edge.data.df))

#### calculate ESOM and clusters ####

library('kohonen')

### build ESOM - don't do this unless you explicitly want to rebuild the map!! ###

#som.grid <- somgrid(xdim = 5, ydim=5, topo="hexagonal")

## ESOM based on pathway abundance

# som.model.paths <- som(path.matrix, 
#                  grid = som.grid, 
#                  rlen = 100,
#                  alpha = c(0.05,0.01),
#                  keep.data = TRUE,
#                  n.hood = "circular",
#                  toroidal = F)
# 
# path.k <- 8
# 
# som.cluster.paths <- kmeans(som.model.paths$codes, centers = path.k)

## ESOM based on edge abundance

# som.model.edges <- som(edge.norm, 
#                  grid = som.grid, 
#                  rlen = 100,
#                  alpha = c(0.05,0.01),
#                  keep.data = TRUE,
#                  n.hood = "circular",
#                  toroidal = T)
# 
# som.cluster.edges <- kmeans(som.model.edges$codes, centers = k)
# 
# plot(som.model.edges,
#      main = '',
#      type = "property",
#      property = som.cluster.edges$cluster,
#      palette.name = topo.colors)
# add.cluster.boundaries(som.model.edges, som.cluster.edges$cluster)

### Cluster the map nodes, edge-based map ###

## Evaluate how many clusters might be in there using simprof

# library(clustsig)
# simprof.codes.edges <- simprof(som.model.edges$codes, method.distance = 'braycurtis', silent = F)

## Evaluate visually using within-clusters sum of squares

# wss.edges <- (nrow(som.model.edges$codes)-1)*sum(apply(som.model.edges$codes,2,var)) 
# for (i in 2:15) {
#   wss.edges[i] <- sum(kmeans(som.model.edges$codes, centers=i)$withinss)
# }
# 
# plot(wss.edges)
# 
# pdf('16S_wss.pdf',
#     width = 5,
#     height = 5)
# 
# plot(wss.edges,
#      pch = 19,
#      ylab = 'Within-clusters sum of squares',
#      xlab = 'K')
# 
# lines(c(8,8),
#       c(0,10000),
#       lty = 2)
# 
# lines(c(0,20),
#       c(wss.edges[8],wss.edges[8]),
#       lty = 2)
# 
# dev.off()

## Pick a reasonable value close to the inflection point and cluster with kmeans

# k <- 8
# som.cluster.edges <- kmeans(som.model.edges$codes, centers = k)

## write out model and clusters so they can be used in future analyses
## don't do this unless you want to overwrite existing!

#save(list = c('som.model.paths', 'som.cluster.paths', 'som.model.edges', 'som.cluster.edges'), file = 'pal_timeseries_som.Rdata')

#### load existing ESOMs, clusters, and make some basic plots ####

load('pal_timeseries_som.Rdata')

## add mode assignments to edge data ##

edge.data.df$mode <- som.cluster.edges$cluster[som.model.edges$unit.classif]
edge.data.df$f.mode <- som.cluster.paths$cluster[som.model.paths$unit.classif]

## Make some descriptive plots

pdf('som_property_map_clusters.pdf', width = 4, height = 4)

plot(som.model.edges,
     main = '',
     type = "property",
     property = som.cluster.edges$cluster,
     palette.name = topo.colors)
add.cluster.boundaries(som.model.edges, som.cluster.edges$cluster)

dev.off()

pdf('som_quality_map.pdf', width = 4, height = 4)
plot(som.model.edges, type = 'quality', pch = 19, palette.name = topo.colors, main = '')
add.cluster.boundaries(som.model.edges, som.cluster.edges$cluster)
dev.off()

pdf('som_mapping_map.pdf', width = 4, height = 4)
plot(som.model.edges, type = 'mapping', pch = 19, palette.name = topo.colors, main = '')
add.cluster.boundaries(som.model.edges, som.cluster.edges$cluster)
dev.off()

## Make some maps of abundance of different edges

property.map <- function(property, model, cluster){
  pdf(paste0('som_property_map_', property, '.pdf'), width = 4, height = 4)
  plot(model,
       main = '',
       type = "property",
       property = model$codes[,which(colnames(model$codes) == property)],
       palette.name = topo.colors)
  add.cluster.boundaries(model, cluster$cluster)
  dev.off()
}

property.map('X424', som.model.edges, som.cluster.edges)
property.map('X93', som.model.edges, som.cluster.edges)
property.map('X210', som.model.edges, som.cluster.edges)
property.map('X500', som.model.edges, som.cluster.edges)
property.map('X82', som.model.edges, som.cluster.edges)
property.map('X4571', som.model.edges, som.cluster.edges)
property.map('X73', som.model.edges, som.cluster.edges)

## Conduct a PCA of codebook vector edge abundance

edge.node.pca <- prcomp(som.model.edges$codes)

edge.node.frac <- colSums(abs(edge.node.pca$x)) / sum(colSums(abs(edge.node.pca$x)))

pdf('16S_node_pca.pdf',
    width = 5,
    height = 5)

plot(edge.node.pca$x[,1], edge.node.pca$x[,2], type = 'n',
     xlab = paste0('PC1: ', {round(edge.node.frac[1] * 100, 2)},' %'),
     ylab = paste0('PC2: ', {round(edge.node.frac[2] * 100, 2)},' %'),
     ylim = c(-0.3, 0.3),
     xlim = c(-0.5, 0.3))

points(edge.node.pca$x[which(som.cluster.edges$cluster == 1),1], edge.node.pca$x[which(som.cluster.edges$cluster == 1),2],
       pch = 19, col = 'lightblue')

points(edge.node.pca$x[which(som.cluster.edges$cluster == 2),1], edge.node.pca$x[which(som.cluster.edges$cluster == 2),2],
       pch = 19, col = 'darkblue')

points(edge.node.pca$x[which(som.cluster.edges$cluster == 3),1], edge.node.pca$x[which(som.cluster.edges$cluster == 3),2],
       pch = 19, col = 'lightgreen')

points(edge.node.pca$x[which(som.cluster.edges$cluster == 4),1], edge.node.pca$x[which(som.cluster.edges$cluster == 4),2],
       pch = 19, col = 'darkgreen')

points(edge.node.pca$x[which(som.cluster.edges$cluster == 5),1], edge.node.pca$x[which(som.cluster.edges$cluster == 5),2],
       pch = 19, col = 'yellow')

points(edge.node.pca$x[which(som.cluster.edges$cluster == 6),1], edge.node.pca$x[which(som.cluster.edges$cluster == 6),2],
       pch = 19, col = 'orange')

points(edge.node.pca$x[which(som.cluster.edges$cluster == 7),1], edge.node.pca$x[which(som.cluster.edges$cluster == 7),2],
       pch = 19, col = 'red')

points(edge.node.pca$x[which(som.cluster.edges$cluster == 8),1], edge.node.pca$x[which(som.cluster.edges$cluster == 8),2],
       pch = 19, col = 'grey')

major.edges <- order(rowSums(abs(edge.node.pca$rotation[,1:2])), decreasing = T)[1:7]

scaling.factor = 0.4

arrows(0,
       0,
       edge.node.pca$rotation[major.edges,1] * scaling.factor,
       edge.node.pca$rotation[major.edges,2] * scaling.factor)

legend('bottomleft',
       bg = 'white',
       legend = c('Mode 1',
                  'Mode 2',
                  'Mode 3',
                  'Mode 4',
                  'Mode 5',
                  'Mode 6',
                  'Mode 7',
                  'Mode 8'),
       pch = 19,
       col = c('lightblue', 'darkblue', 'lightgreen', 'darkgreen', 'yellow', 'orange', 'red', 'grey'))

dev.off()

plot(sort(rowSums(abs(edge.node.pca$rotation[,1:5])), decreasing = T)[1:20],
     type = 'b',
     pch = 19,
     ylab = 'abs(PC1 + PC2)',
     xlab = 'Rank')

row.names(edge.node.pca$rotation)[major.edges]
edge.node.pca$rotation[major.edges,1:2]

colSums(abs(edge.node.pca$rotation[major.edges,1:2])) / colSums(abs(edge.node.pca$rotation[,1:2]))

#### Work with the most abundant edges ####

## Get the most abundant edges, and convert to edge numbers to taxon names

edge.names <- read.csv('edge_taxa.csv', as.is = T)

edge.norm.1314 <- subset(edge.norm, strptime(row.names(edge.norm), '%Y-%m-%d') >= as.POSIXct("2013-07-01"))
edge.norm.1314.abund <- t(edge.norm.1314[,order(colSums(edge.norm), decreasing = T)[0:30]])
row.names(edge.norm.1314.abund) <- sub('X', '', row.names(edge.norm.1314.abund))

known.names <- edge.names$organism_name[match(as.numeric(row.names(edge.norm.1314.abund)), edge.names$clade)]
row.names(edge.norm.1314.abund)[is.na(known.names)]

## These are the CEG, which cannot be determined automatically

known.names[which(row.names(edge.norm.1314.abund) == '212')] <- 'Candidatus Pelagibacter'
known.names[which(row.names(edge.norm.1314.abund) == '93')] <- 'Rhodobacteraceae'
known.names[which(row.names(edge.norm.1314.abund) == '4619')] <- 'Francisella'
known.names[which(row.names(edge.norm.1314.abund) == '412')] <- 'Flavobacteriaceae'
known.names[which(row.names(edge.norm.1314.abund) == '499')] <- 'Flavobacteriales'
known.names[which(row.names(edge.norm.1314.abund) == '4676')] <- 'Colwellia'

row.names(edge.norm.1314.abund) <- known.names

## Make a heatmap of most abundant taxa for 13-14 season, with colorbar showing mode

col.color.palette <- c('lightblue', 'darkblue', 'lightgreen', 'darkgreen', 'yellow', 'orange', 'red', 'grey')
col.color <- col.color.palette[edge.data.df$mode[match(colnames(edge.norm.1314.abund), row.names(edge.data.df))]]

heatcolors <- colorRampPalette(c('white', 'red'))(100)

pdf('mean_edge_abund_1314.pdf',
    width = 8,
    height = 8)

heatmap(edge.norm.1314.abund,
        Colv = NA,
        Rowv = NA,
        col = heatcolors,
        margins = c(10, 23),
        scale = 'n',
        ColSideColors = col.color,
        cexRow = 1.2)

dev.off()

#### Environmental parameters ####

library(zoo)

### BP data ###

pal.st.bp <- read.table('pal_bp.txt', header = T, stringsAsFactors = F, sep = '\t', na.strings = '-999')
pal.st.bp$Date.GMT <- strptime(pal.st.bp$Date.GMT, format = "%m/%d/%Y")
colnames(pal.st.bp) <- c('year', 'date', 'event', 'station', 'depth', 'abundance', 'thy.pmol.l.hr', 'leu.pmol.l.hr')

pal.st.bp.b.10 <- subset(pal.st.bp, station == 'B' & depth == 10)
pal.st.bp.swi <- subset(pal.st.bp, station == 'SWI')
pal.st.bp.swi.b <- rbind(pal.st.bp.b.10, pal.st.bp.swi)

pal.st.bp.swi.b$date <- as.POSIXct(pal.st.bp.swi.b$date)
pal.st.bp.swi.b <- pal.st.bp.swi.b[!duplicated(as.POSIXct(pal.st.bp.swi.b$date)),]

### Flow cytometer data ###

sg.13.14 <- read.csv('SG_PAL1314_esom.cluster.tally.csv', header = T, row.names = 1)
sg.12.13 <- read.csv('SG_PAL1213_esom.cluster.tally.csv', header = T, row.names = 1)
sg.11.12 <- read.csv('SG_PAL1112_esom.cluster.tally.csv', header = T, row.names = 1)
sg.10.11 <- read.csv('SG_PAL1011_esom.cluster.tally.csv', header = T, row.names = 1)

sg <- rbind(sg.13.14, sg.12.13, sg.11.12, sg.10.11)

sg <- sg[grep('PH|B_10', row.names(sg)),]
sg$X4 <- NULL
sg.abund <- rowSums(sg) * (1000/14) # bacterial abundance in ml
sg$abund <- sg.abund
sg$fhna <- sg$X1 / (sg$X1 + sg$X2 + sg$X3) 

sg$date <- sub('SG_PAL1314_|SGr_PAL1314_|SG_PAL1213_|SG_PAL1112_|SG_PAL1011_', '', row.names(sg))
sg$date <- as.POSIXct(strptime(sg$date, '%y%m%d'))

sg <- sg[!duplicated(sg$date),]
sg <- sg[which(!sg$date == as.POSIXct('2012-10-31')),] # remove bad date

### Convert to time series and interpolate ###
## unfortunately need to interpolate seperately for each season, to constrain within seasonal limits

## 2010-2011

sg.10.11.select <- sg[which(sg$date >= as.POSIXct('2010-07-01') & sg$date < as.POSIXct('2011-07-01')),]
sg.10.11.fhna <- zoo(sg.10.11.select$fhna, sg.10.11.select$date)
sg.10.11.abund <- zoo(sg.10.11.select$abund, sg.10.11.select$date)
sg.10.11.time <- merge(sg.10.11.fhna, sg.10.11.abund)
sg.10.11.expand <- merge(sg.10.11.time, zoo(, seq(start(sg.10.11.time), end(sg.10.11.time), by = 'DSTday')), all = T)
sg.10.11.interp <- na.approx(sg.10.11.expand)
colnames(sg.10.11.interp) <- c('fhna', 'abund')

## 2011-2012

sg.11.12.select <- sg[which(sg$date >= as.POSIXct('2011-07-01') & sg$date < as.POSIXct('2012-07-01')),]
sg.11.12.fhna <- zoo(sg.11.12.select$fhna, sg.11.12.select$date)
sg.11.12.abund <- zoo(sg.11.12.select$abund, sg.11.12.select$date)
sg.11.12.time <- merge(sg.11.12.fhna, sg.11.12.abund)
sg.11.12.expand <- merge(sg.11.12.time, zoo(, seq(start(sg.11.12.time), end(sg.11.12.time), by = 'DSTday')), all = T)
sg.11.12.interp <- na.approx(sg.11.12.expand)
colnames(sg.11.12.interp) <- c('fhna', 'abund')

## 2012-2013

sg.12.13.select <- sg[which(sg$date >= as.POSIXct('2012-07-01') & sg$date < as.POSIXct('2013-07-01')),]
sg.12.13.fhna <- zoo(sg.12.13.select$fhna, sg.12.13.select$date)
sg.12.13.abund <- zoo(sg.12.13.select$abund, sg.12.13.select$date)
sg.12.13.time <- merge(sg.12.13.fhna, sg.12.13.abund)
sg.12.13.expand <- merge(sg.12.13.time, zoo(, seq(start(sg.12.13.time), end(sg.12.13.time), by = 'DSTday')), all = T)
sg.12.13.interp <- na.approx(sg.12.13.expand)
colnames(sg.12.13.interp) <- c('fhna', 'abund')

## 2013-2014

sg.13.14.select <- sg[which(sg$date >= as.POSIXct('2013-07-01') & sg$date < as.POSIXct('2014-07-01')),]
sg.13.14.fhna <- zoo(sg.13.14.select$fhna, sg.13.14.select$date)
sg.13.14.abund <- zoo(sg.13.14.select$abund, sg.13.14.select$date)
sg.13.14.time <- merge(sg.13.14.fhna, sg.13.14.abund)
sg.13.14.expand <- merge(sg.13.14.time, zoo(, seq(start(sg.13.14.time), end(sg.13.14.time), by = 'DSTday')), all = T)
sg.13.14.interp <- na.approx(sg.13.14.expand)
colnames(sg.13.14.interp) <- c('fhna', 'abund')

## combine

sg.time <- rbind(sg.10.11.interp, sg.11.12.interp, sg.12.13.interp, sg.13.14.interp)

## Do the same with bp ##

## 2013-2014

pal.st.bp.swi.b.13.14 <- pal.st.bp.swi.b[which(pal.st.bp.swi.b$date >= as.POSIXct("2013-07-01") & pal.st.bp.swi.b$date < as.POSIXct("2014-07-01")),]
bp.13.14.time <- zoo(pal.st.bp.swi.b.13.14$leu.pmol.l.hr, pal.st.bp.swi.b.13.14$date)
bp.13.14.expand <- merge(bp.13.14.time, zoo(, seq(start(bp.13.14.time), end(bp.13.14.time), by = 'DSTday')), all = T)
bp.13.14.interp <- na.approx(bp.13.14.expand)

## 2012-2013

pal.st.bp.swi.b.12.13 <- pal.st.bp.swi.b[which(pal.st.bp.swi.b$date >= as.POSIXct("2012-07-01") & pal.st.bp.swi.b$date < as.POSIXct("2013-07-01")),]
bp.12.13.time <- zoo(pal.st.bp.swi.b.12.13$leu.pmol.l.hr, pal.st.bp.swi.b.12.13$date)
bp.12.13.expand <- merge(bp.12.13.time, zoo(, seq(start(bp.12.13.time), end(bp.12.13.time), by = 'DSTday')), all = T)
bp.12.13.interp <- na.approx(bp.12.13.expand)

## 2011-2012

pal.st.bp.swi.b.11.12 <- pal.st.bp.swi.b[which(pal.st.bp.swi.b$date >= as.POSIXct("2011-07-01") & pal.st.bp.swi.b$date < as.POSIXct("2012-07-01")),]
bp.11.12.time <- zoo(pal.st.bp.swi.b.11.12$leu.pmol.l.hr, pal.st.bp.swi.b.11.12$date)
bp.11.12.expand <- merge(bp.11.12.time, zoo(, seq(start(bp.11.12.time), end(bp.11.12.time), by = 'DSTday')), all = T)
bp.11.12.interp <- na.approx(bp.11.12.expand)

## 2010-2011

pal.st.bp.swi.b.10.11 <- pal.st.bp.swi.b[which(pal.st.bp.swi.b$date >= as.POSIXct("2010-07-01") & pal.st.bp.swi.b$date < as.POSIXct("2011-07-01")),]
bp.10.11.time <- zoo(pal.st.bp.swi.b.10.11$leu.pmol.l.hr, pal.st.bp.swi.b.10.11$date)
bp.10.11.expand <- merge(bp.10.11.time, zoo(, seq(start(bp.10.11.time), end(bp.10.11.time), by = 'DSTday')), all = T)
bp.10.11.interp <- na.approx(bp.10.11.expand)

## 2009-2010

pal.st.bp.swi.b.9.10 <- pal.st.bp.swi.b[which(pal.st.bp.swi.b$date >= as.POSIXct("2009-07-01") & pal.st.bp.swi.b$date < as.POSIXct("2010-07-01")),]
bp.9.10.time <- zoo(pal.st.bp.swi.b.9.10$leu.pmol.l.hr, pal.st.bp.swi.b.9.10$date)
bp.9.10.expand <- merge(bp.9.10.time, zoo(, seq(start(bp.9.10.time), end(bp.9.10.time), by = 'DSTday')), all = T)
bp.9.10.interp <- na.approx(bp.9.10.expand)

## combine

bp.time <- rbind(bp.9.10.interp, bp.10.11.interp, bp.11.12.interp, bp.12.13.interp, bp.13.14.interp)

## Convert mode and functional mode to timeseries

mode.time <- zoo(edge.data.df$mode, order.by = as.POSIXct(edge.data.df$date))
f.mode.time <- zoo(edge.data.df$f.mode, order.by = as.POSIXct(edge.data.df$date))

## Combine all timeseries data

mode.params <- merge(chl.interp, bp.time, mode.time, f.mode.time, sg.time)

### Make a plot of bp and mode for all seasons ###

pdf('bp_mode_all_seasons.pdf',
    width = 6,
    height = 5)

plot(mode.params$bp,
     xlim = c(as.POSIXct('2009-07-01'), as.POSIXct('2014-07-01')),
     lwd = 2,
     col = 'orange',
     xaxt = 'n',
     ylab = expression(BP ~ (pmol ~ leu ~ L^{-1} ~ h^{-1})),
     xlab = 'Date')

p.sf <- max(mode.params$bp, na.rm = T) / 8

for(l in seq(p.sf, 8 * p.sf, p.sf)){
  lines(c(as.POSIXct('2008-07-01'), as.POSIXct('2015-07-01')),
        c(l, l),
        lty = 2,
        col = 'grey')
}

points(mode.params$mode * p.sf,
       pch = 21,
       bg = 'white')

points(na.omit(mode.params)$mode * p.sf,
       pch = 19)

axis(4,
     at = seq(p.sf, 8 * p.sf, p.sf),
     labels = 1:8)

axis.POSIXct(1,
     at = as.POSIXct(c('2010-01-01', '2011-01-01', '2012-01-01', '2013-01-01', '2014-01-01')),
     labels = 2010:2014)

dev.off()

#### Nested linear models ####

## Mode.params.df limited to those observations with bp, chl, and mode.
## Use of a df is necessary to allow factors.

mode.params.df <- as.data.frame(na.omit(mode.params))
colnames(mode.params.df) <- c('chl', 'bp', 'mode', 'f.mode', 'fhna', 'abund')

## Calculate specific bacterial production

mode.params.df$spb <- (mode.params.df$bp / 1000) / mode.params.df$abund

## Make sure modes are factors!

mode.params.df$mode <- as.factor(mode.params.df$mode)
mode.params.df$f.mode <- as.factor(mode.params.df$f.mode)

### stepwise regression for samples that have fhna, bp, mode ###

select.bp.f.mode.lm <- lm(mode.params.df$bp ~ mode.params.df$f.mode)
select.bp.mode.lm <- lm(mode.params.df$bp ~ mode.params.df$mode)
select.bp.hna.lm <- lm(mode.params.df$bp ~ mode.params.df$fhna)
select.bp.mode.hna.lm <- lm(mode.params.df$bp ~ mode.params.df$mode + mode.params.df$fhna)

select.bp.abund.lm <- lm(mode.params.df$bp ~ mode.params.df$abund)
select.bp.abund.mode.lm <- lm(mode.params.df$bp ~ mode.params.df$abund + mode.params.df$mode)
select.bp.abund.mode.fhna.lm <- lm(mode.params.df$bp ~ mode.params.df$abund + mode.params.df$mode + mode.params.df$fhna)
select.bp.abund.f.mode.fhna.lm <- lm(mode.params.df$bp ~ mode.params.df$abund + mode.params.df$f.mode + mode.params.df$fhna)


anova(select.bp.f.mode.lm, select.bp.hna.lm)
anova(select.bp.mode.lm, select.bp.hna.lm)
anova(select.bp.mode.hna.lm, select.mode.hna.lm)
anova(select.bp.abund.lm, select.bp.mode.hna.lm)
anova(select.bp.abund.mode.lm, select.bp.abund.lm)
anova(select.bp.abund.mode.lm, select.bp.abund.mode.fhna.lm)
anova(select.bp.abund.mode.fhna.lm, select.bp.abund.f.mode.fhna.lm)

summary(select.bp.mode.hna.lm)

## Calculate variance inflation factors

library(cars)

vif(select.bp.abund.mode.fhna.lm)

summary(select.bp.abund.lm)
summary(select.bp.abund.mode.lm)
summary(select.bp.abund.mode.fhna.lm) ## top scoring model

bp.aic <- AIC(select.bp.f.mode.lm, select.bp.mode.lm, select.bp.hna.lm, select.bp.mode.hna.lm, select.bp.abund.lm, select.bp.abund.mode.lm, select.bp.abund.mode.fhna.lm)
bp.aic.min <- min(bp.aic[,2])

for(aic in bp.aic[,2]){
  print(exp((bp.aic.min - aic) / 2))
}

### Nested models for specific bacterial production ###

select.sbp.f.mode.lm <- lm(mode.params.df$spb ~ mode.params.df$f.mode)
select.sbp.mode.lm <- lm(mode.params.df$spb ~ mode.params.df$mode)
select.sbp.hna.lm <- lm(mode.params.df$spb ~ mode.params.df$fhna)
select.sbp.mode.hna.lm <- lm(mode.params.df$spb ~ mode.params.df$mode + mode.params.df$fhna)

select.sbp.abund.lm <- lm(mode.params.df$spb ~ mode.params.df$abund)
select.sbp.abund.mode.lm <- lm(mode.params.df$spb ~ mode.params.df$abund + mode.params.df$mode)
select.sbp.abund.mode.fhna.lm <- lm(mode.params.df$spb ~ mode.params.df$abund + mode.params.df$mode + mode.params.df$fhna)
select.sbp.abund.f.mode.fhna.lm <- lm(mode.params.df$spb ~ mode.params.df$abund + mode.params.df$f.mode + mode.params.df$fhna)

summary(select.sbp.f.mode.lm)
anova(select.sbp.hna.lm, select.sbp.f.mode.lm)

sbp.aic <- AIC(select.sbp.f.mode.lm, select.sbp.mode.lm, select.sbp.hna.lm, select.sbp.mode.hna.lm, select.sbp.abund.lm, select.sbp.abund.mode.lm, select.sbp.abund.mode.fhna.lm)
sbp.aic.min <- min(sbp.aic[,2])

for(aic in sbp.aic[,2]){
  print(exp((sbp.aic.min - aic) / 2))
}


### Plots showing correlation for parameters used in stepwise regression ###

plot.colors <- c('lightblue', 'darkblue', 'lightgreen', 'darkgreen', 'yellow', 'orange', 'red', 'grey')
plot.pch <- c(0,1,2,3,4,5,6,8)

## Specific bacterial production

pdf('SBP_mode_select.pdf',
    width = 6,
    height = 5)

plot(mode.params.df$spb ~ mode.params.df$mode,
     xlab = 'Mode',
     ylab = NA,
     yaxt = 'n')

axis(2, at = c(0, 1e-7, 2e-7),
     labels = c(0, 1e-7, 2e-7))

#ylab = expression(BP ~ (pmol ~ leu ~ L^{-1} ~ h^{-1})))

dev.off()

pdf('SBP_fhna_select.pdf',
    width = 6,
    height = 5)

plot(mode.params.df$spb ~ mode.params.df$fhna,
     #ylab = expression(BP ~ (pmol ~ leu ~ L^{-1} ~ h^{-1}))
     xlab = 'fHNA',
     ylab = NA,
     pch = 19,
     type = 'n',
     yaxt = 'n')

axis(2, at = c(0, 1e-7, 2e-7),
     labels = c(0, 1e-7, 2e-7))

m = 0
for(mode in 1:8){
  m = m + 1
  try({points(mode.params.df$spb[which(as.character(mode.params.df$mode) == mode)] ~ 
         mode.params.df$fhna[which(as.character(mode.params.df$mode) == mode)],
         #col = plot.colors[m],
         pch = plot.pch[m])}, silent = T)
}

legend('topright',
       bg = 'white',
       legend = c('Mode 1',
                  'Mode 2',
                  'Mode 3',
                  'Mode 4',
                  'Mode 5',
                  'Mode 6',
                  'Mode 7',
                  'Mode 8'),
       #col = plot.colors,
       pch = plot.pch)

dev.off()

##

pdf('SBP_abund_select.pdf',
    width = 6,
    height = 5)

plot(mode.params.df$spb ~ mode.params.df$abund,
     #ylab = expression(BP ~ (pmol ~ leu ~ L^{-1} ~ h^{-1}))
     xlab = expression(Bacteria ~ mL^{-1}),
     ylab = NA,
     pch = 19,
     type = 'n',
     xaxt = 'n',
     yaxt = 'n',
     xlim = c(0, 1.5e6))

axis(1, at = c(0, 500000, 1e6, 1.5e6),
     labels = c(0, '5e5', '1e6', '1.5e6'))

axis(2, at = c(0, 1e-7, 2e-7),
     labels = c(0, 1e-7, 2e-7))

m = 0
for(mode in 1:8){
  m = m + 1
  try({points(mode.params.df$spb[which(as.character(mode.params.df$mode) == mode)] ~ 
           mode.params.df$abund[which(as.character(mode.params.df$mode) == mode)],
         #col = plot.colors[m],
         pch = plot.pch[m])}, silent = T)
}

legend('topright',
       bg = 'white',
       legend = c('Mode 1',
                  'Mode 2',
                  'Mode 3',
                  'Mode 4',
                  'Mode 5',
                  'Mode 6',
                  'Mode 7',
                  'Mode 8'),
       #col = plot.colors,
       pch = plot.pch)

dev.off()

## bacterial production

pdf('BP_mode_select.pdf',
    width = 6,
    height = 5)

plot(mode.params.df$bp ~ mode.params.df$mode,
     xlab = 'Mode',
     ylab = NA)
#ylab = expression(BP ~ (pmol ~ leu ~ L^{-1} ~ h^{-1})))

dev.off()

pdf('BP_abund_select.pdf',
    width = 6,
    height = 5)

plot(mode.params.df$bp ~ mode.params.df$abund,
     xlab = expression(Bacteria ~ mL^{-1}),
     ylab = NA,
     pch = 19,
     xaxt = 'n',
     xlim = c(0,1.5e6),
     type = 'n')

axis(1, at = c(0, 500000, 1e6, 1.5e6),
     labels = c(0, '5e5', '1e6', '1.5e6'))

abline(lm(mode.params.df$bp ~ mode.params.df$abund))

m = 0
for(mode in 1:8){
  m = m + 1
  try({points(mode.params.df$bp[which(as.character(mode.params.df$mode) == mode)] ~ 
                mode.params.df$abund[which(as.character(mode.params.df$mode) == mode)],
              #col = plot.colors[m],
              pch = plot.pch[m])}, silent = T)
}

legend('topleft',
       bg = 'white',
       legend = c('Mode 1',
                  'Mode 2',
                  'Mode 3',
                  'Mode 4',
                  'Mode 5',
                  'Mode 6',
                  'Mode 7',
                  'Mode 8'),
       pch = plot.pch)

dev.off()

##

pdf('BP_fhna_select.pdf',
    width = 6,
    height = 5)

plot(mode.params.df$bp ~ mode.params.df$fhna,
     xlab = 'fHNA',
     ylab = NA,
     pch = 19,
     type = 'n')

abline(lm(mode.params.df$bp ~ mode.params.df$fhna))

m = 0
for(mode in 1:8){
  m = m + 1
  try({points(mode.params.df$bp[which(as.character(mode.params.df$mode) == mode)] ~ 
           mode.params.df$fhna[which(as.character(mode.params.df$mode) == mode)],
         #col = plot.colors[m],
         pch = plot.pch[m])}, silent = T)
}

legend('topleft',
       bg = 'white',
       legend = c('Mode 1',
                  'Mode 2',
                  'Mode 3',
                  'Mode 4',
                  'Mode 5',
                  'Mode 6',
                  'Mode 7',
                  'Mode 8'),
       pch = plot.pch)

dev.off()

#### Plots for hna/lna ####

### absolute Pelagibacter abundance ###

pelagibacter.1314.time <- zoo(edge.norm.1314.abund[1,], order.by = as.POSIXct(colnames(edge.norm.1314.abund)))
pelagibacter.1314.abund.time <- merge.zoo(pelagibacter.1314.time, sg.13.14.expand)
pelagibacter.1314.abund.time <- na.approx(pelagibacter.1314.abund.time)
plot(pelagibacter.1314.abund.time$pelagibacter.1314.time * pelagibacter.1314.abund.time$sg.13.14.abund)

### absolute Rhodo abundance ###

rhodo.1314.time <- zoo(edge.norm.1314.abund[5,], order.by = as.POSIXct(colnames(edge.norm.1314.abund)))
rhodo.1314.abund.time <- merge.zoo(rhodo.1314.time, sg.13.14.expand)
rhodo.1314.abund.time <- na.approx(rhodo.1314.abund.time)
plot(rhodo.1314.abund.time$rhodo.1314.time * rhodo.1314.abund.time$sg.13.14.abund)

pdf('pal_1314_abund.pdf',
    width = 8,
    height = 5)

plot(sg.13.14.interp$abund,
     ylab = expression(Bacterial ~ abundance ~ (ml^{-1})),
     #ylab = ''
     xlab = NA,
     xaxt = 'n',
     yaxt = 'n')

axis(2,
     labels = c('5e5', '1e6', '1.5e6'),
     at = c('5e5', '1e6', '1.5e6'))

points(sg.13.14.time$sg.13.14.abund,
       pch = 19)

points((pelagibacter.1314.abund.time$pelagibacter.1314.time * pelagibacter.1314.abund.time$sg.13.14.abund),
       lwd = 2,
       lty = 1,
       type = 'l',
       col = 'red')

points((rhodo.1314.abund.time$rhodo.1314.time * rhodo.1314.abund.time$sg.13.14.abund),
       lwd = 2,
       lty = 1,
       type = 'l',
       col = 'green')

points((dokdonia.1314.abund.time$dokdonia.1314.time * dokdonia.1314.abund.time$sg.13.14.abund),
       lwd = 2,
       lty = 1,
       type = 'l',
       col = 'yellow')

for(x in c(as.POSIXct(c('2013-11-01', '2013-12-01', '2014-01-01', '2014-02-01', '2014-03-01')))){
  lines(c(x, x), c(-2, 2),
        lwd = 1,
        lty = 'dotted',
        col = 'lightgrey')
}

box()

legend('topleft',
       legend = c('All bacteria', 'Candidatus P. ubique HTCC1062', 'Dokdonia sp. MED134', 'Rhodobacteraceae'),
       lty = 1,
       lwd = c(1, 2, 2, 2),
       col = c('black', 'red', 'yellow', 'green'),
       bg = 'white')

dev.off()

pdf('pal_1314_fhna.pdf',
    width = 8,
    height = 5)

plot(sg.13.14.interp$fhna,
     ylab = 'fHNA',
     xlab = NA,
     xaxt = 'n')

points(sg.13.14.time$sg.13.14.fhna,
       pch = 19)

for(x in c(as.POSIXct(c('2013-11-01', '2013-12-01', '2014-01-01', '2014-02-01', '2014-03-01')))){
  lines(c(x, x), c(-2, 2),
        lwd = 1,
        lty = 'dotted',
        col = 'lightgrey')
}

dev.off()

### Calculate abundance and fHNA anomalies, then delta anomaly ###

sg.13.14.abund.anom <- (sg.13.14.interp$abund - mean(sg.13.14.interp$abund)) / mean(sg.13.14.interp$abund)
sg.13.14.fhna.anom <- (sg.13.14.interp$fhna - mean(sg.13.14.interp$fhna)) / mean(sg.13.14.interp$fhna)
sg.13.14.delta.anomaly <- sg.13.14.abund.anom - sg.13.14.fhna.anom
sg.13.14.delta.anomaly.inst <- sg.13.14.delta.anomaly[index(sg.13.14.delta.anomaly) %in% index(sg.13.14.time),]

pdf('fhna_abund_anomaly.pdf',
    width = 8,
    height = 5)

plot(sg.13.14.delta.anomaly,
     ylab = expression(paste(Delta, " anomaly")),
     xlab = 'Date',
     ylim = c(-1, 1.5),
     col = 'blue',
     type = 'n',
     lwd = 4,
     xaxt = 'n')

points(sg.13.14.delta.anomaly[which(sg.13.14.delta.anomaly > 0)],
       col = 'blue',
       type = 'h',
       lwd = 4)

points(sg.13.14.delta.anomaly[which(sg.13.14.delta.anomaly < 0)],
     col = 'orange',
     type = 'h',
     lwd = 4)

points(sg.13.14.delta.anomaly.inst[which(sg.13.14.delta.anomaly.inst > 0)],
       type = 'h',
       lwd = 1)

points(sg.13.14.delta.anomaly.inst[which(sg.13.14.delta.anomaly.inst < 0)],
       type = 'h',
       lwd = 1)

axis.POSIXct(1, at = c(as.POSIXct(c('2013-11-01', '2013-12-01', '2014-01-01', '2014-02-01', '2014-03-01'))),
                       labels = c('Nov', 'Dec', 'Jan', 'Feb', 'March'))

abline(0,0)

for(x in c(as.POSIXct(c('2013-11-01', '2013-12-01', '2014-01-01', '2014-02-01', '2014-03-01')))){
  lines(c(x, x), c(-2, 2),
        lwd = 1,
        lty = 'dotted',
        col = 'lightgrey')
}

dev.off()

#### Community succession and persistence ####

## Now that we've defined modes is there evidence of succession?  tally which mode normally follows a given mode
## row is leading mode (m), column is following mode (m + 1)

succession.matrix <- matrix(ncol = k, nrow = k, data = 0)

for(i in 1:length(edge.data.df$mode) - 1){
  m <- edge.data.df$mode[i]
  m.1 <- edge.data.df$mode[i + 1]
  print(paste(m, m.1))
  succession.matrix[m, m.1] <- succession.matrix[m, m.1] + 1
}

succession.matrix.rel <- succession.matrix / rowSums(succession.matrix, na.rm = T)

## Evaluate the probability of this distribution occurring randomly with permutations.

permutation.matrix <- matrix(ncol = k, nrow = k, data = 0)

for(r in 1:100){
  print(paste(r))
  temp.matrix <- matrix(ncol = k, nrow = k, data = 0)

  for(i in 1:length(edge.data.df$mode) - 1){
    m <- edge.data.df$mode[i]
    m.1 <- sample(1:k, 1)
    temp.matrix[m, m.1] <- temp.matrix[m, m.1] + 1
  }
  
  ## for each random transition, if the number of observed transitions exceeded the number of random transitions
  ## add "1" to the number of that transition in the permutation matrix
  
  permutation.matrix[which(temp.matrix >= succession.matrix)] <- permutation.matrix[which(temp.matrix >= succession.matrix)] + 1
}

permutation.matrix[succession.matrix == 0] <- NA
succession.matrix[succession.matrix == 0] <- NA

### Make a nice heatmap of succession, with borders around significant cells ###

heat.color <- colorRampPalette(c('blue', 'red'))(100)

## if the number of transitions in the data exceeded the number of random transitions in 95% or greater of
## events, this is deemed signficant

sig.permutations <- matrix(ncol = k, nrow = k)

for(i in 1:k){
  for(j in 1:k){
    try({
      if(permutation.matrix[i,j] <= 5){
        sig.permutations[i,j] <- permutation.matrix[i,j]
      }}, silent = T)
    }
  }

library(gplots)

pdf('succession_matrix.pdf', width = 4, height = 4)

heatmap.2(succession.matrix,
        Rowv = NA,
        Colv=NA,
        xlab = 'Following mode',
        ylab = 'Leading mode',
        col = heat.color,
        na.rm = T,
        bg = 'gray',
        scale = 'none',
        key = F,
        na.color = ('gray'),
        colsep = 1:8,
        rowsep = 1:8,
        trace = 'none',
        cellnote = succession.matrix,
        notecol = 'white',
        notecex = 1.2)

dev.off()

pdf('permutation_matrix.pdf', width = 4, height = 4)

heatmap.2(sig.permutations,
          Rowv = NA,
          Colv=NA,
          xlab = 'Following mode',
          ylab = 'Leading mode',
          col = heat.color,
          na.rm = T,
          bg = 'gray',
          scale = 'none',
          key = F,
          na.color = ('gray'),
          colsep = 1:8,
          rowsep = 1:8,
          trace = 'none',
          cellnote = sig.permutations,
          notecol = 'white',
          notecex = 1.2)

dev.off()

#### Calculate and plot mean genomic characteristics for all samples by mode ####

## Mean predicted genome length

pdf('genome_length_mode.pdf',
    width = 6,
    height = 4)

boxplot(edge.data.df$genome_size.mean ~ edge.data.df$mode,
        yaxt = 'n',
        ylab = 'Predicted mean genome length (Mbp)',
        xlab = 'Mode')
order(tapply(edge.data.df$genome_size.mean, edge.data.df$mode, mean))
axis(2,
     at = c(2e6, 3e6, 4e6, 5e6),
     labels = c('2', '3', '4', '5'))

dev.off()

## Mean 16S rRNA gene copy number

pdf('n16S_mode.pdf',
    width = 6,
    height = 4)

boxplot(edge.data.df$n16S.mean ~ edge.data.df$mode,
        ylab = 'Predicted mean 16S rRNA gene copy number',
        xlab = 'Mode')

dev.off()


boxplot(edge.data.df$ncds.mean ~ edge.data.df$mode)

## n16S ~ genome size

pdf('n16S_genome_length.pdf',
    width = 6,
    height = 4)

gc.color <- colorRampPalette(c('white', 'red'))(100)

plot(edge.data.df$n16S.mean ~ edge.data.df$genome_size.mean,
     xlab = 'Predicted mean genome length (bp)',
     ylab = 'Predicted mean 16S rRNA gene copy number',
     xaxt = 'n',
     type = 'n')

abline(lm(edge.data.df$n16S.mean ~ edge.data.df$genome_size.mean))

axis(1,
     at = c(2e6, 3e6, 4e6, 5e6),
     labels = c('2e6', '3e6', '4e6', '5e6'))

points(edge.data.df$n16S.mean ~ edge.data.df$genome_size.mean,
       pch = 21,
       bg = gc.color[cut(edge.data.df$GC.mean, breaks = length(edge.data.df$GC.mean))],
       cex = 1.2)

# for(i in 1:length(edge.data.df$genome_size.mean)){
#   lines(rbind(c(edge.data.df$genome_size.mean[i], edge.data.df$n16S.mean[i]),
#               c(edge.data.df$genome_size.mean[i], edge.data.df$n16S.mean[i] + edge.data.df$n16S.sd[i])))
# }

dev.off()