library(tidyverse)
library(ggpubr)
library(ggtranscript)

#Get and format data
plotdata <- read.table('../Results/HP202.plotdata.2.3.tsv', sep = '\t', header = F)
colnames(plotdata) <- c('Rname', 'Position', 'Genotype', 'DelType')

#Subset by readnames for plotting
rname_df <- read.table('../Results/HP202.2.3.pos.rnames.txt', sep = '\t', header = F)
colnames(rname_df) <- c('Rname', 'Status')
rnames <-  c(rname_df$Rname)
plotdata <- subset(plotdata, Rname %in% rnames)
plotdata <- left_join(plotdata, rname_df, by = "Rname")


#Calculate values for path plotting
startsite <- function(incoord, x) {
  startpos <- incoord[x] - ((incoord[x] - incoord[x -1])/2)
  return(startpos)
}

plotdata$Start <- 0

for (x in 1:nrow(plotdata)) {
  if (x > 1 && x < nrow(plotdata)) {
    if (plotdata$Rname[x] == plotdata$Rname[x-1]) {
      plotdata$Start[x] <- startsite(plotdata$Position, x)
    }
    else {
      plotdata$Start[x] <- plotdata$Position[x]
    }
  }
  else {
    plotdata$Start[x] <- plotdata$Position[x]
  }
}

#replace last start value, for the end of the path plotting
plotdata$Start[plotdata$Position == 16276] <- 16276


#The plot
posdata <- plotdata[plotdata$Status == "True", ]

p <- ggplot(posdata, aes(Position,Rname,group=Rname,col=Genotype)) + 
  geom_point(size = 1.5) +
  geom_path(aes(Start,Rname,group=Rname,col=Genotype), linewidth = 1) +
  theme_minimal() +
  scale_color_manual(values=c("Orange", "Blue")) +
  theme(axis.text.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_continuous(name = 'Position',
                     breaks = c(unique(plotdata$Position)),
                     labels = c(unique(plotdata$Position)),
                     guide = guide_axis(check.overlap = T),
                     limits = c(1, 16276)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Reads") +
  theme(legend.position = "top")
p

annot <- read.table('../Results/NZB.summary.tsv', sep = '\t', header = F)
annot$V10 <- (annot$V2 + annot$V3)/2

a <- ggplot(annot) +
  geom_range(aes(xstart = V2, xend = V3, y = V1, fill = V6, )) +
  theme_minimal() +
  scale_x_continuous(name = 'Gene',
                     breaks = c(annot$V10),
                     labels = c(annot$V7),
                     #guide = guide_axis(check.overlap = T),
                     limits = c(1, 16276)) +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill=guide_legend(title="Gene type")) +
  ylab("Annotation") +
  theme(axis.text.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'none')

fig <- ggarrange(p, a,
                 labels = c("Recombinant reads", ""),
                 ncol = 1, nrow = 2,
                 heights = c(8, 1))
fig

ggsave(file="HP202_positive.annotated.pdf", plot = p, width = 210, height = 297, units = "mm")
