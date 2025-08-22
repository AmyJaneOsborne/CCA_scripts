library("dplyr")
library("ggplot2")
library(grid)
`%notin%` <- Negate(`%in%`)

datam <- read.csv("gene_matrix.csv", header=T)
datam$name_study  = ifelse(datam$name_study %in% "PBMC", "Peripheral blood\nmononuclear cell", datam$name_study)

names(datam) <- c("Gene_symbol","study_ID","name_study","MN_vs")
datam.noempty <- datam[datam$study_ID %notin% c("",NA),]

plot <- ggplot(datam.noempty, aes(x = MN_vs, y = Gene_symbol, fill = name_study)) +
  geom_tile(color = "white")+
  #coord_fixed()+
  scale_fill_manual(values = c("grey80", "grey65","grey35","grey8"))+
  labs(x = "Membranous nephropathy vs.", y = "Gene", fill = "Cell type")+
  theme(
	panel.background = element_blank(),
	axis.line = element_line(colour = "black"))
	#panel.grid.minor = element_line(color = "grey"),
	#panel.grid.major = element_line(color = "grey"))
	#panel.background = element_rect(fill = "white", colour = "white"))
  
  # breaks = levels(datam.noempty$name_study),
  
jpeg("gene_matrix_out.jpeg", pointsize=4, width=250, height=250,units="mm",res = 294)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,heights=c(1))))
print(plot, vp=viewport(layout.pos.row=1,layout.pos.col=1))


  theme(
    legend.text=element_text(size=10),
	legend.key.size = unit(0.5, "cm"),
	legend.spacing.y = unit(0.5, 'cm'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill="white"),
    panel.border = element_rect(fill=NA, colour = "black", linetype="solid", size = 1),
    axis.text.x = element_text(family="Arial", face="bold", size = 7, angle=45,vjust = 1, hjust=1,color="black"),
	  axis.text.y = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.x = element_text(family="Arial",face="bold", size = 8,color="black"),
	  axis.title.y = element_text(family="Arial",face="bold", size = 8,color="black",margin = margin(t = 0, r = 12, b = 0, l = 0))
    #plot.title = element_text(family="Arial",face="bold", size = 8,color="black",hjust=0.5),
    #plot.tag = element_text(face="bold", size = 20,color="black")
)

