## Before using this script you need to download all the fcs files from
## the palmer datazoo at http://oceaninformatics.ucsd.edu/datazoo/data/pallter/datasets

#### parameters ####

f.list <- list.files(pattern = "SG.*fcs") # all fcs files matching this pattern will be analyzed
output <- 'SG_PAL_esom'
model <- 'SG_pal_timeseries_esom.Rdata'

FL1.H.limit <- 3500 # lower bound for FL3.H, below 1e4 is usually noise

col = c('grey', 'red', 'green', 'blue', 'orange', 'yellow', 'black', 'purple', 'lightblue', 'lightgreen', 'pink')

#### aggregation and QC ####

fcs.results <- matrix(ncol = 9, nrow = 0)
colnames(fcs.results) <- c('file', 'cluster', "mean.FSC.H", "mean.SSC.H", "mean.FL1.H", "mean.FL2.H", "mean.FL3.H", "mean.FL4.H", 'size')

all.events <- data.frame(FSC.H = numeric(),
                         SSC.H = numeric(),
                         FL1.H = numeric(),
                         FL2.H = numeric(),
                         FL3.H = numeric(),
                         FL4.H = numeric(),
                         cluster = character(),
                         sample = character())

library('flowCore')

f.name <- f.list[1]

for(f.name in f.list){
  print(f.name)
  
  fcm <- read.FCS(f.name)
  colnames(fcm) <- sub('-', '.', colnames(fcm))
  fcm <- as.data.frame((exprs(fcm)))
  
  fcm$FSC.A <- NULL
  fcm$SSC.A <- NULL
  fcm$FL1.A <- NULL
  fcm$FL2.A <- NULL
  fcm$FL3.A <- NULL
  fcm$FL4.A <- NULL
  fcm$Width <- NULL
  fcm$Time <- NULL
  
  fcm[fcm == 0] <- NA
  fcm <- na.omit(fcm)
  
  fcm <- fcm[fcm$FL1.H > FL1.H.limit,]    ## remove noise with low SG fluorescence
  fcm <- fcm[fcm$FL1.H > 5 * fcm$FL2.H,] ## if FL1 >> FL2 not likely to be biological
  
  plot(fcm$FL1.H ~ fcm$SSC.H,
       log = 'xy')
  
  lines(c(1, 1e6),
        c(5, 5e6))


  fcm$sample <- sub('.fcs', '', f.name)
  all.events <- rbind(all.events, fcm)
  write.csv(fcm, sub('fcs', 'csv', f.name), quote = F)
}

write.csv(all.events, paste0(output, '.events.csv'), quote = F)

#### esom ####

library(kohonen)

events <- read.csv(paste0(output, '.events.csv'), stringsAsFactors = F)
events$cluster <- NA

#### training - don't retrain model unless you really want to ####

# sample.size <- 50000 # size of the training set
# 
# events.sample <- events[sample(1:length(events[,1]), sample.size),]
# sample.mat <- cbind(events.sample$FSC.H, events.sample$SSC.H, events.sample$FL1.H, events.sample$FL2.H, events.sample$FL3.H, events.sample$FL4.H)
# sample.mat <- log10(sample.mat)
# 
# grid.size <- ceiling(sample.size ^ (1/2.5)) # in practice this grid size works well, produces 40 x 40 grid for 10000 events
# 
# som_grid <- somgrid(xdim = grid.size, ydim = grid.size, topo="hexagonal")
# 
# som_model <- som(sample.mat, 
#                  grid = som_grid, 
#                  rlen = 100,
#                  alpha = c(0.05,0.01),
#                  keep.data = TRUE,
#                  n.hood = "circular",
#                  toroidal = F)
# 
# ## use kmeans and determine best number of clusters
# 
# som.events <- som_model$codes 
# wss <- (nrow(som.events) - 1) * sum(apply(som.events, 2, var)) 
# 
# for (i in 2:15) {
#   wss[i] <- sum(kmeans(som.events, centers=i)$withinss)
# }
# 
# pdf('fcm_wss.pdf',
#     width = 5,
#     height = 5)
# 
# plot(wss,
#      pch = 19,
#      ylab = 'Within-clusters sum of squares',
#      xlab = 'K')
# 
# lines(c(4,4),
#       c(0,10000),
#       lty = 2)
# 
# lines(c(0,20),
#       c(wss[4],wss[4]),
#       lty = 2)
# 
# dev.off()
# 
# ## Pick reasonable value near inflection point
# 
# k <- 4
# som_cluster <- kmeans(som.events, centers = k)
# 
# ## Descriptive plots of ESOM
# 
# pdf(paste0(output, '.properties.pdf'),
#     width = 6,
#     height = 5)
# 
# plot(som_model, type = 'counts')
# add.cluster.boundaries(som_model, som_cluster$cluster)
# 
# plot(som_model, type = 'property', property = (som.events[,1]), main = 'log10(FSC.H)')
# add.cluster.boundaries(som_model, som_cluster$cluster, lwd = 2)
# 
# plot(som_model, type = 'property', property = (som.events[,2]), main = 'log10(SSC.H)')
# add.cluster.boundaries(som_model, som_cluster$cluster)
# 
# plot(som_model, type = 'property', property = (som.events[,3]), main = 'log10(FL1.H)')
# add.cluster.boundaries(som_model, som_cluster$cluster)
# 
# plot(som_model, type = 'property', property = (som.events[,4]), main = 'log10(FL2.H)')
# add.cluster.boundaries(som_model, som_cluster$cluster)
# 
# plot(som_model, type = 'property', property = (som.events[,5]), main = 'log10(FL3.H)')
# add.cluster.boundaries(som_model, som_cluster$cluster)
# 
# plot(som_model, type = 'property', property = (som.events[,6]), main = 'log10(FL4.H)')
# add.cluster.boundaries(som_model, som_cluster$cluster)
# 
# dev.off()
# 
# ## property plot for pub
# 
# pdf('flh1_property.pdf',
#     width = 5,
#     height = 5)
# 
# plot(som_model, type = 'property', property = (som.events[,3]),
#      main = "",
#      heatkeywidth = 1.2,
#      palette.name = terrain.colors)
# 
# add.cluster.boundaries(som_model, som_cluster$cluster, lwd = 2)
# 
# dev.off()
# 
# ## pca on nodes
# 
# node.pca <- prcomp(som_model$codes)
# 
# pdf('fcm_node_pca.pdf',
#     width = 5,
#     height = 5)
# 
# plot(node.pca$x[,1], node.pca$x[,2], type = 'n',
#      xlab = 'PC1',
#      ylab = 'PC2')
# 
# points(node.pca$x[which(som_cluster$cluster == 1),1], node.pca$x[which(som_cluster$cluster == 1),2],
#        pch = 19, col = 'blue')
# 
# points(node.pca$x[which(som_cluster$cluster == 2),1], node.pca$x[which(som_cluster$cluster == 2),2],
#        pch = 19, col = 'red')
# 
# points(node.pca$x[which(som_cluster$cluster == 3),1], node.pca$x[which(som_cluster$cluster == 3),2],
#        pch = 19, col = 'green')
# 
# points(node.pca$x[which(som_cluster$cluster == 4),1], node.pca$x[which(som_cluster$cluster == 4),2],
#        pch = 19, col = 'grey')
# 
# scaling.factor = 3
# 
# arrows(0,
#        0,
#        node.pca$rotation[,1] * scaling.factor,
#        node.pca$rotation[,2] * scaling.factor)
# 
# legend('bottomleft',
#        bg = 'white',
#        legend = c('HNA',
#                   'MNA',
#                   'LNA',
#                   'Other'),
#        pch = 19,
#        col = c('blue', 'red', 'green', 'grey'))
# 
# dev.off()

### don't save unless you want to override the existing mapping! ###

#save(list = c('som_model', 'som_cluster', 'k'), file = model)

#### Classify all events ####

load(model)

cluster.tally <- matrix(nrow = length(unique(events$sample)), ncol = k)
colnames(cluster.tally) <- 1:k
row.names(cluster.tally) <- unique(events$sample)

pdf(paste0(output, '.clusters.pdf'))

for(sample in unique(events$sample)){
  print(sample)
  events.sample <- events[which(events$sample == sample),]
  sample.mat <- cbind(events.sample$FSC.H, events.sample$SSC.H, events.sample$FL1.H, events.sample$FL2.H, events.sample$FL3.H, events.sample$FL4.H)
  sample.mat <- log10(sample.mat)
  
  sample.predict <- predict(som_model, newdata = sample.mat, trainY = som_cluster$cluster[som_model$unit.classif])
  events$cluster[which(events$sample == sample)] <- sample.predict$prediction
  
  ylim <- c(5e3, 1.25 * max(events$FL1.H))
  xlim <- c(10, 1.25 * max(events$SSC.H))
  
  plot(events.sample$FL1.H ~ events.sample$FSC.H,
       type = 'n',
       log = 'xy',
       main = sample,
       xlab = 'FSC.H',
       ylab = 'FL1.H',
       ylim = ylim,
       xlim = xlim)
  
  for(cluster in 1:k){
    i = which(events$sample == sample & events$cluster == cluster)
    cluster.tally[sample, cluster] <- length(i)
    points(events$FL1.H[i] ~ events$FSC.H[i],
           pch = 19,
           col = col[cluster])
  }
  
  legend('topleft',
         c('cluster1', 'cluster2', 'cluster3', 'cluster4'),
         pch = 19,
         col = c(col[1], col[2], col[3], col[4]))
}

dev.off()

write.csv(cluster.tally, paste0(output, '.cluster.tally.csv'), quote = F)