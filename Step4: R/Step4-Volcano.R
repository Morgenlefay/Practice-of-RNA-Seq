###--------------------------
rm(list=ls())
options(stringsAsFactors = F)

###火山图-1
rm(list=ls())
library(ggplot2)
volcano_data <-  read.csv("A315T_DEG.csv",header = TRUE)
loc_up <- intersect(which(volcano_data$padj<0.05),which(volcano_data$log2FoldChange>=1))
loc_down <- intersect(which(volcano_data$padj<0.05),which(volcano_data$log2FoldChange<=(-1)))
significant <- rep("Normal",times=nrow(volcano_data))
significant[loc_up] <- "Up"
significant[loc_down] <- "Down"
significant <- factor(significant,levels=c("Up","Down","Normal"))
p <- qplot(x=volcano_data$log2FoldChange,y=-log10(volcano_data$padj),xlab="log2(FC)",ylab="-log10(FDR)",colour=significant)
p <- p+ scale_color_manual(values=c("Up"="red","Normal"="black","Down"="green"))
xline=c(-log2(2),log2(2))
p <- p+geom_vline(xintercept=xline,lty=2,size=I(0.2),colour="grey11")
yline=-log(0.05,10)
p <- p+geom_hline(yintercept=yline,lty=2,size=I(0.2),colour="grey11")
p <- p+theme_bw()
p

#火山图-2
rm(list = ls())
library(ggplot2)
data <- read.csv("A315T_DEG.csv",header = TRUE)
data$color <- ifelse(data$padj<0.05 & abs(data$log2FoldChange)>= 1,ifelse(data$log2FoldChange > 1,'red','blue'),'black')
color <- c(red = "red",black = "black",blue = "blue")
p <- ggplot(data, aes(log2FoldChange, -log10(padj), col = color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (q-value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
p

###火山图-3
rm(list = ls())
library(EnhancedVolcano)
data <- read.csv("A315T_DEG.csv",header = TRUE, row.names = 1)
EnhancedVolcano(data,
                lab = rownames(data),
                x = "log2FoldChange",
                y = "padj",
                selectLab = c("ENSG00000182752","ENSG00000137558"),
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 0.00000000000000001,
                FCcutoff = 2.0,
                xlim = c(-8,8),
                transcriptLabSize = 3.0,
                colAlpha = 1,
                legend=c("NS","Log2 FC","Adjusted p-value",
                         "Adjusted p-value & Log2 FC"),
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0,
                DrawConnectors = FALSE,
                border = "full",
                borderWidth = 1.5,
                borderColour = "black",
                gridlines.major = FALSE,
                gridlines.minor = FALSE)
