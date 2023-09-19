library(ggplot2)

data <- read.csv("table4_for_plot.csv", header=T)
names(data) <- c("analysis","term","cat","adj_pval","genes","count_gene")
data <- data[!is.na(data$adj_pval),]
data$term <- as.factor(data$term)
label_both = c("Expression quantitative trait loci\ngenes in CKDGen, NURTURE-CKD\nand Salford Kidney Study",
"Genes for 43 novel missense variants\nin CKDGen and BioBank Japan")
names(label_both) <- unique(data$analysis)
data$analysis <- factor(data$analysis, levels = unique(data$analysis), labels = c("Expression quantitative trait loci\ngenes in CKDGen, NURTURE-CKD\nand Salford Kidney Study",
"Genes for 43 novel missense variants\nin CKDGen and BioBank Japan"))
				  
plot <- ggplot(data, aes(y = reorder(stringr::str_wrap(term,35),count_gene), x = count_gene, fill = adj_pval, size = adj_pval)) +
		geom_dotplot(binaxis='y', stackdir='center',  stackratio=1.5, dotsize=1.2)+
		xlab("Gene count")+ 
		#xlab("Enrichment P-value (adjusted)")+ 
        ylab("Gene function/pathway term")+
		guides(fill=guide_legend(title="Enrichment\nP-value\n(adjusted)"))+	
		facet_wrap(vars(analysis))+
		theme(
			strip.background = element_rect( strip.background=element_rect(colour="black", fill="white"))
		)
