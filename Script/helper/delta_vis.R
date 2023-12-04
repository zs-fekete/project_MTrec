library('tidyverse')
library('ggpubr')
library(ggtranscript)


#Get and format data
plotdata <- read.table('../Results/BH10.full.23-10-31.plotdata.tsv', sep = '\t', header = F)
colnames(plotdata) <- c('Rname', 'Position', 'Genotype', 'DelType')

#Subset by readnames for plotting
rname_df <- read.table('reads.txt', sep = '\t', header = F)
colnames(rname_df) <- c('Rname', 'Status')
rnames <-  c(rname_df$Rname)
plotdata <- subset(plotdata, Rname %in% rnames)
plotdata <- left_join(plotdata, rname_df, by = "Rname")

#Get data for delta reads
deltanames <- subset(plotdata, DelType == 'delta')
deltanames <- unique(deltanames$Rname)

deltaloc <- c(10522, 15487)
Rname <- rep(c(deltanames), each = 2)
loc <- rep(deltaloc, length(Rname))
deltable <- data.frame(Rname, loc)
deltable <- left_join(deltable, rname_df, by = "Rname")
#deltable$DelType <- 'delta'

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
plotdata$Start[plotdata$Position == 16207] <- 16207 


#The plot
posdata <- plotdata[plotdata$Status == "True", ]
posdel <- deltable[deltable$Status == "True", ]
negdata <- plotdata[plotdata$Status == "False", ]
negdel <- deltable[deltable$Status == "False", ]

p <- ggplot(posdata, aes(Position,Rname,group=Rname,col=Genotype)) + 
  geom_point(size = 2.5) +
  geom_path(aes(Start,Rname,group=Rname,col=Genotype), linewidth = 1.5) +
  theme_minimal() +
  scale_color_manual(values=c("Orange", "Blue")) +
  geom_path(data = posdel, aes(loc, Rname, group = Rname),
            colour = 'lightgrey', linewidth = 1.5) +
  theme(axis.text.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_continuous(name = 'Position',
                     breaks = c(unique(plotdata$Position), deltaloc),
                     labels = c(unique(plotdata$Position), deltaloc),
                     guide = guide_axis(check.overlap = T),
                     limits = c(1, 16207)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Reads") +
  theme(legend.position = "top")

n <- ggplot(negdata, aes(Position,Rname,group=Rname,col=Genotype)) + 
  geom_point(size = 2.5) +
  geom_path(aes(Start,Rname,group=Rname,col=Genotype), linewidth = 1.5) +
  theme_minimal() +
  scale_color_manual(values=c("Orange", "Blue"),
                     guide = guide_legend(override.aes = list(color = "white"))) +
  geom_path(data = negdel, aes(loc, Rname, group = Rname),
            colour = 'lightgrey', linewidth = 1.5) +
  theme(axis.text.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
scale_x_continuous(name = 'Position',
                     breaks = c(unique(plotdata$Position), deltaloc),
                     labels = c(unique(plotdata$Position), deltaloc),
                     guide = guide_axis(check.overlap = T),
                     limits = c(1, 16207)) +
  theme(axis.text.x = element_text(angle = 90),
        legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent")) +
  ylab("Reads") +
  theme(legend.position = "top")


annot <- read.table('../Results/BH10.summary.tsv', sep = '\t', header = F)
annot$V10 <- (annot$V2 + annot$V3)/2

a <- ggplot(annot) +
  geom_range(aes(xstart = V2, xend = V3, y = V1, fill = V6, )) +
  theme_minimal() +
  scale_x_continuous(name = 'Gene',
                     breaks = c(annot$V10),
                     labels = c(annot$V7),
                     guide = guide_axis(check.overlap = T),
                     limits = c(1, 16207)) +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(fill=guide_legend(title="Gene type")) +
  ylab("Annotation") +
  theme(axis.text.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = 'none')
  
    



fig <- ggarrange(p, n,
                labels = c("Recombinant", "Non-recombinant"),
                ncol = 1, nrow = 2)
#fig
f2 <- ggarrange(fig, a,
                ncol = 1, nrow = 2,
                heights = c(3, 1))
f2
