# for homozygous SNPs only:

`%notin%` <- Negate(`%in%`)

data <- read.table("data.vcf",	comment.char = "",skip=30,header=T)
	
newnames <- gsub("^.{0,141}", "", names(data[,c(10:2562)]))
data <- data[,c(1:5,10:2562)]
names(data) <- c("chr","pos","RSID","ref","alt",newnames)
GTs <- data
GTs$no.hom.0 <-  rowSums(GTs[,c(6:2558)] == "0/0")  
GTs$no.het.01 <- rowSums(GTs[,c(6:2558)] == "0/1") 
GTs$no.het.10 <- rowSums(GTs[,c(6:2558)] == "1/0") 
GTs$no.het <- GTs$no.het.01+GTs$no.het.10
GTs$no.hom.1 <- rowSums(GTs[,c(6:2558)] == "1/1") 
GTs$no.missing <- rowSums(GTs[,c(6:2558)] == "./.")  
GTs$total <- GTs$no.hom.0 + GTs$no.het + GTs$no.hom.1 + GTs$no.missing 

GTs$alt.af <- (GTs$no.het + (GTs$no.hom.1*2))/(GTs$total*2)
GTs$alt.ac <- GTs$no.het + (GTs$no.hom.1*2)
GTs$an <- GTs$total*2

GTs[,c(1:5,2559:2568)]
	
CKD_PTs <- read.csv("data.csv", header=T)

membr <- CKD_PTs[CKD_PTs$primary_diagnosis_Membranous.nephropathy...idiopathic ==1,]
GTs.membr <- GTs[,names(GTs) %in% membr$SentrixPosition]

GTs.membr$no.hom.0 <-  rowSums(GTs.membr[,c(1:ncol(GTs.membr))] == "0/0")  
GTs.membr$no.het.01 <- rowSums(GTs.membr[,c(1:ncol(GTs.membr))] == "0/1")  
GTs.membr$no.het.10 <- rowSums(GTs.membr[,c(1:ncol(GTs.membr))] == "1/0") 
GTs.membr$no.het <- GTs.membr$no.het.01+GTs.membr$no.het.10
GTs.membr$no.hom.1 <- rowSums(GTs.membr[,c(1:ncol(GTs.membr))] == "1/1")  
GTs.membr$no.missing <- rowSums(GTs.membr[,c(1:ncol(GTs.membr))] == "./.") 
GTs.membr$total <- GTs.membr$no.hom.0 + GTs.membr$no.het + GTs.membr$no.hom.1 + GTs.membr$no.missing 

GTs.membr$alt.af <- (GTs.membr$no.het + (GTs.membr$no.hom.1*2))/(GTs.membr$total*2)
GTs.membr$alt.ac <- GTs.membr$no.het + (GTs.membr$no.hom.1*2)
GTs.membr$an <- GTs.membr$total*2
GTs.membr$PT <- "GTs.membr"

# membr - combinations of up to 16 SNPs:
results_combo_Membr <- data.frame(no_SNPs = numeric(), SNPs = character(), overlap = numeric(), no_overlap = numeric(), stringsAsFactors=FALSE)

for(m in 1:6){
	print(paste0("combo length: ",m))
	combos <- combinations(n=16, r=m, v=c(1:16), repeats.allowed = FALSE)
	# n Size of the source vector, r=Size of the target vectors, v=Source vector. Defaults to 1:n
	for(x in 1:nrow(combos)){
		print(paste0("x: ",x))
		rowno = nrow(results_combo_Membr)+1
		idSnp = combos[x,]
		if(m==1){
			out = names(GTs.membr[idSnp[1],])[which(GTs.membr[idSnp[1],] %in% c("1/1"))]
			
		}else if(m==2){
			out = intersect(names(GTs.membr[idSnp[1],])[which(GTs.membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.membr[idSnp[2],])[which(GTs.membr[idSnp[2],] %in% c("1/1"))])
		}else if(m==3){
			out = intersect(intersect(names(GTs.membr[idSnp[1],])[which(GTs.membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.membr[idSnp[2],])[which(GTs.membr[idSnp[2],] %in% c("1/1"))]),
					names(GTs.membr[idSnp[3],])[which(GTs.membr[idSnp[3],] %in% c("1/1"))])
			
		}else if(m==4){
			out = intersect(intersect(intersect(names(GTs.membr[idSnp[1],])[which(GTs.membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.membr[idSnp[2],])[which(GTs.membr[idSnp[2],] %in% c("1/1"))]),
					names(GTs.membr[idSnp[3],])[which(GTs.membr[idSnp[3],] %in% c("1/1"))]),
						names(GTs.membr[idSnp[4],])[which(GTs.membr[idSnp[4],] %in% c("1/1"))])
			
		}else if (m==5){
			out = intersect(intersect(intersect(intersect(names(GTs.membr[idSnp[1],])[which(GTs.membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.membr[idSnp[2],])[which(GTs.membr[idSnp[2],] %in% c("1/1"))]),
					names(GTs.membr[idSnp[3],])[which(GTs.membr[idSnp[3],] %in% c("1/1"))]),
						names(GTs.membr[idSnp[4],])[which(GTs.membr[idSnp[4],] %in% c("1/1"))]),
							names(GTs.membr[idSnp[5],])[which(GTs.membr[idSnp[5],] %in% c("1/1"))])
				
		}else if(m==6){
			out = intersect(intersect(intersect(intersect(intersect(names(GTs.membr[idSnp[1],])[which(GTs.membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.membr[idSnp[2],])[which(GTs.membr[idSnp[2],] %in% c("1/1"))]),
					names(GTs.membr[idSnp[3],])[which(GTs.membr[idSnp[3],] %in% c("1/1"))]),
						names(GTs.membr[idSnp[4],])[which(GTs.membr[idSnp[4],] %in% c("1/1"))]),
							names(GTs.membr[idSnp[5],])[which(GTs.membr[idSnp[5],] %in% c("1/1"))]),
								names(GTs.membr[idSnp[6],])[which(GTs.membr[idSnp[6],] %in% c("1/1"))])
			
		}
		results_combo_Membr[rowno,1] = length(idSnp)
		results_combo_Membr[rowno,2] = paste(data$RSID[idSnp], collapse="_")
		results_combo_Membr[rowno,3] = paste0(out, collapse="; ")
		results_combo_Membr[rowno,4] = length(out)

	}
	x=0
}

results_combo_Membr$freq_MN = results_combo_Membr$no_overlap/87

### =======
# GTs - non-membr:
non_membr <- CKD_PTs[CKD_PTs$primary_diagnosis_Membranous.nephropathy...idiopathic ==0,]
GTs.non_membr <- GTs[,names(GTs) %in% non_membr$SentrixPosition]
# 87
GTs.non_membr$no.hom.0 <-  rowSums(GTs.non_membr[,c(1:ncol(GTs.non_membr))] == "0/0")   
GTs.non_membr$no.het.01 <- rowSums(GTs.non_membr[,c(1:ncol(GTs.non_membr))] == "0/1")  
GTs.non_membr$no.het.10 <- rowSums(GTs.non_membr[,c(1:ncol(GTs.non_membr))] == "1/0")  
GTs.non_membr$no.het <- GTs.non_membr$no.het.01+GTs.non_membr$no.het.10
GTs.non_membr$no.hom.1 <- rowSums(GTs.non_membr[,c(1:ncol(GTs.non_membr))] == "1/1")  
GTs.non_membr$no.missing <- rowSums(GTs.non_membr[,c(1:ncol(GTs.non_membr))] == "./.")  
GTs.non_membr$total <- GTs.non_membr$no.hom.0 + GTs.non_membr$no.het + GTs.non_membr$no.hom.1 + GTs.non_membr$no.missing  

GTs.non_membr$alt.af <- (GTs.non_membr$no.het + (GTs.non_membr$no.hom.1*2))/(GTs.non_membr$total*2)
GTs.non_membr$alt.ac <- GTs.non_membr$no.het + (GTs.non_membr$no.hom.1*2)
GTs.non_membr$an <- GTs.non_membr$total*2
GTs.non_membr$PT <- "GTs.non_membr"

# non-membr - combinations of up to 6 SNPs:
results_combo_nonMembr <- data.frame(no_SNPs = numeric(), SNPs = character(), overlap = numeric(), no_overlap = numeric(), stringsAsFactors=FALSE)

for(m in 1:6){
	print(paste0("combo length: ",m))
	combos <- combinations(n=16, r=m, v=c(1:16), repeats.allowed = FALSE)
	# n Size of the source vector, r=Size of the target vectors, v=Source vector. Defaults to 1:n
	for(x in 1:nrow(combos)){
		print(paste0("x: ",x))
		rowno = nrow(results_combo_nonMembr)+1
		idSnp = combos[x,]
		if(m==1){
			out = names(GTs.non_membr[idSnp[1],])[which(GTs.non_membr[idSnp[1],] %in% c("1/1"))]
			
		}else if(m==2){
			out = intersect(names(GTs.non_membr[idSnp[1],])[which(GTs.non_membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr[idSnp[2],])[which(GTs.non_membr[idSnp[2],] %in% c("1/1"))])
		}else if(m==3){
			out = intersect(intersect(names(GTs.non_membr[idSnp[1],])[which(GTs.non_membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr[idSnp[2],])[which(GTs.non_membr[idSnp[2],] %in% c("1/1"))]),
					names(GTs.non_membr[idSnp[3],])[which(GTs.non_membr[idSnp[3],] %in% c("1/1"))])
			
		}else if(m==4){
			out = intersect(intersect(intersect(names(GTs.non_membr[idSnp[1],])[which(GTs.non_membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr[idSnp[2],])[which(GTs.non_membr[idSnp[2],] %in% c("1/1"))]),
					names(GTs.non_membr[idSnp[3],])[which(GTs.non_membr[idSnp[3],] %in% c("1/1"))]),
						names(GTs.non_membr[idSnp[4],])[which(GTs.non_membr[idSnp[4],] %in% c("1/1"))])
			
		}else if (m==5){
			out = intersect(intersect(intersect(intersect(names(GTs.non_membr[idSnp[1],])[which(GTs.non_membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr[idSnp[2],])[which(GTs.non_membr[idSnp[2],] %in% c("1/1"))]),
					names(GTs.non_membr[idSnp[3],])[which(GTs.non_membr[idSnp[3],] %in% c("1/1"))]),
						names(GTs.non_membr[idSnp[4],])[which(GTs.non_membr[idSnp[4],] %in% c("1/1"))]),
							names(GTs.non_membr[idSnp[5],])[which(GTs.non_membr[idSnp[5],] %in% c("1/1"))])
				
		}else if(m==6){
			out = intersect(intersect(intersect(intersect(intersect(names(GTs.non_membr[idSnp[1],])[which(GTs.non_membr[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr[idSnp[2],])[which(GTs.non_membr[idSnp[2],] %in% c("1/1"))]),
					names(GTs.non_membr[idSnp[3],])[which(GTs.non_membr[idSnp[3],] %in% c("1/1"))]),
						names(GTs.non_membr[idSnp[4],])[which(GTs.non_membr[idSnp[4],] %in% c("1/1"))]),
							names(GTs.non_membr[idSnp[5],])[which(GTs.non_membr[idSnp[5],] %in% c("1/1"))]),
								names(GTs.non_membr[idSnp[6],])[which(GTs.non_membr[idSnp[6],] %in% c("1/1"))])
			
		}
		results_combo_nonMembr[rowno,1] = length(idSnp)
		results_combo_nonMembr[rowno,2] = paste(data$RSID[idSnp], collapse="_")
		results_combo_nonMembr[rowno,3] = paste0(out, collapse="; ")
		results_combo_nonMembr[rowno,4] = length(out)

	}
	x=0
}

results_combo_nonMembr$freq_nonMN = results_combo_nonMembr$no_overlap/817


### =======
# GTs - non-membr and non T1DN:
non_membr_nonT1DN <- CKD_PTs[CKD_PTs$primary_diagnosis_Membranous.nephropathy...idiopathic ==0 & CKD_PTs$primary_diagnosis_Diabetic.nephropathy.in.type.I.diabetes...no.histology ==0,]
GTs.non_membr_nonT1DN <- GTs[,names(GTs) %in% non_membr_nonT1DN$SentrixPosition]

GTs.non_membr_nonT1DN$no.hom.0 <-  rowSums(GTs.non_membr_nonT1DN[,c(1:ncol(GTs.non_membr_nonT1DN))] == "0/0") 
GTs.non_membr_nonT1DN$no.het.01 <- rowSums(GTs.non_membr_nonT1DN[,c(1:ncol(GTs.non_membr_nonT1DN))] == "0/1")  
GTs.non_membr_nonT1DN$no.het.10 <- rowSums(GTs.non_membr_nonT1DN[,c(1:ncol(GTs.non_membr_nonT1DN))] == "1/0") 
GTs.non_membr_nonT1DN$no.het <- GTs.non_membr_nonT1DN$no.het.01+GTs.non_membr_nonT1DN$no.het.10
GTs.non_membr_nonT1DN$no.hom.1 <- rowSums(GTs.non_membr_nonT1DN[,c(1:ncol(GTs.non_membr_nonT1DN))] == "1/1") 
GTs.non_membr_nonT1DN$no.missing <- rowSums(GTs.non_membr_nonT1DN[,c(1:ncol(GTs.non_membr_nonT1DN))] == "./.") 
GTs.non_membr_nonT1DN$total <- GTs.non_membr_nonT1DN$no.hom.0 + GTs.non_membr_nonT1DN$no.het + GTs.non_membr_nonT1DN$no.hom.1 + GTs.non_membr_nonT1DN$no.missing  #353

GTs.non_membr_nonT1DN$alt.af <- (GTs.non_membr_nonT1DN$no.het + (GTs.non_membr_nonT1DN$no.hom.1*2))/(GTs.non_membr_nonT1DN$total*2)
GTs.non_membr_nonT1DN$alt.ac <- GTs.non_membr_nonT1DN$no.het + (GTs.non_membr_nonT1DN$no.hom.1*2)
GTs.non_membr_nonT1DN$an <- GTs.non_membr_nonT1DN$total*2
GTs.non_membr_nonT1DN$PT <- "GTs.non_membr_nonT1DN"

# non-membr - combinations of up to 6 SNPs:
results_combo_nonMembr_nonT1DN <- data.frame(no_SNPs = numeric(), SNPs = character(), overlap = numeric(), no_overlap = numeric(), stringsAsFactors=FALSE)

for(m in 1:6){
	print(paste0("combo length: ",m))
	combos <- combinations(n=16, r=m, v=c(1:16), repeats.allowed = FALSE)
	# n Size of the source vector, r=Size of the target vectors, v=Source vector. Defaults to 1:n
	for(x in 1:nrow(combos)){
		print(paste0("x: ",x))
		rowno = nrow(results_combo_nonMembr_nonT1DN)+1
		idSnp = combos[x,]
		if(m==1){
			out = names(GTs.non_membr_nonT1DN[idSnp[1],])[which(GTs.non_membr_nonT1DN[idSnp[1],] %in% c("1/1"))]
			
		}else if(m==2){
			out = intersect(names(GTs.non_membr_nonT1DN[idSnp[1],])[which(GTs.non_membr_nonT1DN[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr_nonT1DN[idSnp[2],])[which(GTs.non_membr_nonT1DN[idSnp[2],] %in% c("1/1"))])
		}else if(m==3){
			out = intersect(intersect(names(GTs.non_membr_nonT1DN[idSnp[1],])[which(GTs.non_membr_nonT1DN[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr_nonT1DN[idSnp[2],])[which(GTs.non_membr_nonT1DN[idSnp[2],] %in% c("1/1"))]),
					names(GTs.non_membr_nonT1DN[idSnp[3],])[which(GTs.non_membr_nonT1DN[idSnp[3],] %in% c("1/1"))])
			
		}else if(m==4){
			out = intersect(intersect(intersect(names(GTs.non_membr_nonT1DN[idSnp[1],])[which(GTs.non_membr_nonT1DN[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr_nonT1DN[idSnp[2],])[which(GTs.non_membr_nonT1DN[idSnp[2],] %in% c("1/1"))]),
					names(GTs.non_membr_nonT1DN[idSnp[3],])[which(GTs.non_membr_nonT1DN[idSnp[3],] %in% c("1/1"))]),
						names(GTs.non_membr_nonT1DN[idSnp[4],])[which(GTs.non_membr_nonT1DN[idSnp[4],] %in% c("1/1"))])
			
		}else if (m==5){
			out = intersect(intersect(intersect(intersect(names(GTs.non_membr_nonT1DN[idSnp[1],])[which(GTs.non_membr_nonT1DN[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr_nonT1DN[idSnp[2],])[which(GTs.non_membr_nonT1DN[idSnp[2],] %in% c("1/1"))]),
					names(GTs.non_membr_nonT1DN[idSnp[3],])[which(GTs.non_membr_nonT1DN[idSnp[3],] %in% c("1/1"))]),
						names(GTs.non_membr_nonT1DN[idSnp[4],])[which(GTs.non_membr_nonT1DN[idSnp[4],] %in% c("1/1"))]),
							names(GTs.non_membr_nonT1DN[idSnp[5],])[which(GTs.non_membr_nonT1DN[idSnp[5],] %in% c("1/1"))])
				
		}else if(m==6){
			out = intersect(intersect(intersect(intersect(intersect(names(GTs.non_membr_nonT1DN[idSnp[1],])[which(GTs.non_membr_nonT1DN[idSnp[1],] %in% c("1/1"))],
				names(GTs.non_membr_nonT1DN[idSnp[2],])[which(GTs.non_membr_nonT1DN[idSnp[2],] %in% c("1/1"))]),
					names(GTs.non_membr_nonT1DN[idSnp[3],])[which(GTs.non_membr_nonT1DN[idSnp[3],] %in% c("1/1"))]),
						names(GTs.non_membr_nonT1DN[idSnp[4],])[which(GTs.non_membr_nonT1DN[idSnp[4],] %in% c("1/1"))]),
							names(GTs.non_membr_nonT1DN[idSnp[5],])[which(GTs.non_membr_nonT1DN[idSnp[5],] %in% c("1/1"))]),
								names(GTs.non_membr_nonT1DN[idSnp[6],])[which(GTs.non_membr_nonT1DN[idSnp[6],] %in% c("1/1"))])
			
		}
		results_combo_nonMembr_nonT1DN[rowno,1] = length(idSnp)
		results_combo_nonMembr_nonT1DN[rowno,2] = paste(data$RSID[idSnp], collapse="_")
		results_combo_nonMembr_nonT1DN[rowno,3] = paste0(out, collapse="; ")
		results_combo_nonMembr_nonT1DN[rowno,4] = length(out)

	}
	x=0
}

results_combo_nonMembr_nonT1DN$freq_nonMN = results_combo_nonMembr_nonT1DN$no_overlap/770

##########

# add in PLA2r data:
antipla2r <- read.csv("antiPLA2R_status.csv", header=T)
names(antipla2r) <- c("patient_id", "apla2r_status","Any_other_info.")

phenD_pla2r <- merge(membr, antipla2r, by = c("patient_id"))

GTs.antiPLA2rPositive <- GTs.membr[,names(GTs.membr) %in% phenD_pla2r[phenD_pla2r$apla2r_status %in% 1,3],]

GTs.antiPLA2rPositive$no.hom.0 <-  rowSums(GTs.antiPLA2rPositive[,c(1:ncol(GTs.antiPLA2rPositive))] == "0/0")
GTs.antiPLA2rPositive$no.het.01 <- rowSums(GTs.antiPLA2rPositive[,c(1:ncol(GTs.antiPLA2rPositive))] == "0/1")
GTs.antiPLA2rPositive$no.het.10 <- rowSums(GTs.antiPLA2rPositive[,c(1:ncol(GTs.antiPLA2rPositive))] == "1/0")
GTs.antiPLA2rPositive$no.het <- GTs.antiPLA2rPositive$no.het.01+GTs.antiPLA2rPositive$no.het.10
GTs.antiPLA2rPositive$no.hom.1 <- rowSums(GTs.antiPLA2rPositive[,c(1:ncol(GTs.antiPLA2rPositive))] == "1/1")
GTs.antiPLA2rPositive$no.missing <- rowSums(GTs.antiPLA2rPositive[,c(1:ncol(GTs.antiPLA2rPositive))] == "./.")
GTs.antiPLA2rPositive$total <- GTs.antiPLA2rPositive$no.hom.0 + GTs.antiPLA2rPositive$no.het + GTs.antiPLA2rPositive$no.hom.1 + GTs.antiPLA2rPositive$no.missing  #353

GTs.antiPLA2rPositive$alt.af <- (GTs.antiPLA2rPositive$no.het + (GTs.antiPLA2rPositive$no.hom.1*2))/(GTs.antiPLA2rPositive$total*2)
GTs.antiPLA2rPositive$alt.ac <- GTs.antiPLA2rPositive$no.het + (GTs.antiPLA2rPositive$no.hom.1*2)
GTs.antiPLA2rPositive$an <- GTs.antiPLA2rPositive$total*2
GTs.antiPLA2rPositive$PT <- "GTs.antiPLA2rPositive"
GTs.antiPLA2rPositive$SNP <- data$RSID

GTs.antiPLA2rNegative <- GTs.membr[,names(GTs.membr) %in% phenD_pla2r[phenD_pla2r$apla2r_status %in% 0,3],]
# 18
GTs.antiPLA2rNegative$no.hom.0 <-  rowSums(GTs.antiPLA2rNegative[,c(1:ncol(GTs.antiPLA2rNegative))] == "0/0") 
GTs.antiPLA2rNegative$no.het.01 <- rowSums(GTs.antiPLA2rNegative[,c(1:ncol(GTs.antiPLA2rNegative))] == "0/1") 
GTs.antiPLA2rNegative$no.het.10 <- rowSums(GTs.antiPLA2rNegative[,c(1:ncol(GTs.antiPLA2rNegative))] == "1/0")  
GTs.antiPLA2rNegative$no.het <- GTs.antiPLA2rNegative$no.het.01+GTs.antiPLA2rNegative$no.het.10
GTs.antiPLA2rNegative$no.hom.1 <- rowSums(GTs.antiPLA2rNegative[,c(1:ncol(GTs.antiPLA2rNegative))] == "1/1")  
GTs.antiPLA2rNegative$no.missing <- rowSums(GTs.antiPLA2rNegative[,c(1:ncol(GTs.antiPLA2rNegative))] == "./.")  
GTs.antiPLA2rNegative$total <- GTs.antiPLA2rNegative$no.hom.0 + GTs.antiPLA2rNegative$no.het + GTs.antiPLA2rNegative$no.hom.1 + GTs.antiPLA2rNegative$no.missing  #353

GTs.antiPLA2rNegative$alt.af <- (GTs.antiPLA2rNegative$no.het + (GTs.antiPLA2rNegative$no.hom.1*2))/(GTs.antiPLA2rNegative$total*2)
GTs.antiPLA2rNegative$alt.ac <- GTs.antiPLA2rNegative$no.het + (GTs.antiPLA2rNegative$no.hom.1*2)
GTs.antiPLA2rNegative$an <- GTs.antiPLA2rNegative$total*2
GTs.antiPLA2rNegative$PT <- "GTs.antiPLA2rNegative"
GTs.antiPLA2rNegative$SNP <- data$RSID

#write.table(data.frame(GTs.antiPLA2rNegative), "antipla2rNeg_1pla2r_and_5HLA_lead_ALLPredmissense_SNPs_AFs_updated.txt", col.names=F)
#write.table(data.frame(GTs.antiPLA2rPositive), "antipla2rPos_1pla2r_and_5HLA_lead_ALLPredmissense_SNPs_AFs_updated.txt", col.names=F)
#write.table(data.frame(GTs.membr), "1pla2r_and_5HLA_lead_ALLPredmissense_SNPs_AFs_updated.txt", col.names=F)

results_combo <- data.frame(no_SNPs = numeric(), SNPs = character(), overlap_positive = numeric(), no_overlap_positive = numeric(), 
	overlap_negative = numeric(), no_overlap_negative = numeric(), stringsAsFactors=FALSE)

for(m in 1:6){
	print(paste0("combo length: ",m))
	combos <- combinations(n=16, r=m, v=c(1:16), repeats.allowed = FALSE)
	# n Size of the source vector, r=Size of the target vectors, v=Source vector. Defaults to 1:n
	for(x in 1:nrow(combos)){
		print(paste0("x: ",x))
		rowno = nrow(results_combo)+1
		idSnp = combos[x,]
		if(m==1){
			out = names(GTs.antiPLA2rPositive[idSnp[1],])[which(GTs.antiPLA2rPositive[idSnp[1],] %in% c("1/1"))]
			out_n = names(GTs.antiPLA2rNegative[idSnp[1],])[which(GTs.antiPLA2rNegative[idSnp[1],] %in% c("1/1"))]
		}else if(m==2){
			out = intersect(names(GTs.antiPLA2rPositive[idSnp[1],])[which(GTs.antiPLA2rPositive[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rPositive[idSnp[2],])[which(GTs.antiPLA2rPositive[idSnp[2],] %in% c("1/1"))])
			out_n = intersect(names(GTs.antiPLA2rNegative[idSnp[1],])[which(GTs.antiPLA2rNegative[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rNegative[idSnp[2],])[which(GTs.antiPLA2rNegative[idSnp[2],] %in% c("1/1"))])
		}else if(m==3){
			out = intersect(intersect(names(GTs.antiPLA2rPositive[idSnp[1],])[which(GTs.antiPLA2rPositive[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rPositive[idSnp[2],])[which(GTs.antiPLA2rPositive[idSnp[2],] %in% c("1/1"))]),
					names(GTs.antiPLA2rPositive[idSnp[3],])[which(GTs.antiPLA2rPositive[idSnp[3],] %in% c("1/1"))])
			out_n = intersect(intersect(names(GTs.antiPLA2rNegative[idSnp[1],])[which(GTs.antiPLA2rNegative[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rNegative[idSnp[2],])[which(GTs.antiPLA2rNegative[idSnp[2],] %in% c("1/1"))]),
					names(GTs.antiPLA2rNegative[idSnp[3],])[which(GTs.antiPLA2rNegative[idSnp[3],] %in% c("1/1"))])
		}else if(m==4)		{
			out = intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp[1],])[which(GTs.antiPLA2rPositive[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rPositive[idSnp[2],])[which(GTs.antiPLA2rPositive[idSnp[2],] %in% c("1/1"))]),
					names(GTs.antiPLA2rPositive[idSnp[3],])[which(GTs.antiPLA2rPositive[idSnp[3],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp[4],])[which(GTs.antiPLA2rPositive[idSnp[4],] %in% c("1/1"))])
			out_n = intersect(intersect(intersect(names(GTs.antiPLA2rNegative[idSnp[1],])[which(GTs.antiPLA2rNegative[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rNegative[idSnp[2],])[which(GTs.antiPLA2rNegative[idSnp[2],] %in% c("1/1"))]),
					names(GTs.antiPLA2rNegative[idSnp[3],])[which(GTs.antiPLA2rNegative[idSnp[3],] %in% c("1/1"))]),
						names(GTs.antiPLA2rNegative[idSnp[4],])[which(GTs.antiPLA2rNegative[idSnp[4],] %in% c("1/1"))])
		}else if (m==5){
			out = intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp[1],])[which(GTs.antiPLA2rPositive[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rPositive[idSnp[2],])[which(GTs.antiPLA2rPositive[idSnp[2],] %in% c("1/1"))]),
					names(GTs.antiPLA2rPositive[idSnp[3],])[which(GTs.antiPLA2rPositive[idSnp[3],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp[4],])[which(GTs.antiPLA2rPositive[idSnp[4],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp[5],])[which(GTs.antiPLA2rPositive[idSnp[5],] %in% c("1/1"))])
			out_n = intersect(intersect(intersect(intersect(names(GTs.antiPLA2rNegative[idSnp[1],])[which(GTs.antiPLA2rNegative[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rNegative[idSnp[2],])[which(GTs.antiPLA2rNegative[idSnp[2],] %in% c("1/1"))]),
					names(GTs.antiPLA2rNegative[idSnp[3],])[which(GTs.antiPLA2rNegative[idSnp[3],] %in% c("1/1"))]),
						names(GTs.antiPLA2rNegative[idSnp[4],])[which(GTs.antiPLA2rNegative[idSnp[4],] %in% c("1/1"))]),
							names(GTs.antiPLA2rNegative[idSnp[5],])[which(GTs.antiPLA2rNegative[idSnp[5],] %in% c("1/1"))])			
		}else if(m==6){
			out = intersect(intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp[1],])[which(GTs.antiPLA2rPositive[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rPositive[idSnp[2],])[which(GTs.antiPLA2rPositive[idSnp[2],] %in% c("1/1"))]),
					names(GTs.antiPLA2rPositive[idSnp[3],])[which(GTs.antiPLA2rPositive[idSnp[3],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp[4],])[which(GTs.antiPLA2rPositive[idSnp[4],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp[5],])[which(GTs.antiPLA2rPositive[idSnp[5],] %in% c("1/1"))]),
								names(GTs.antiPLA2rPositive[idSnp[6],])[which(GTs.antiPLA2rPositive[idSnp[6],] %in% c("1/1"))])
			out_n = intersect(intersect(intersect(intersect(intersect(names(GTs.antiPLA2rNegative[idSnp[1],])[which(GTs.antiPLA2rNegative[idSnp[1],] %in% c("1/1"))],
				names(GTs.antiPLA2rNegative[idSnp[2],])[which(GTs.antiPLA2rNegative[idSnp[2],] %in% c("1/1"))]),
					names(GTs.antiPLA2rNegative[idSnp[3],])[which(GTs.antiPLA2rNegative[idSnp[3],] %in% c("1/1"))]),
						names(GTs.antiPLA2rNegative[idSnp[4],])[which(GTs.antiPLA2rNegative[idSnp[4],] %in% c("1/1"))]),
							names(GTs.antiPLA2rNegative[idSnp[5],])[which(GTs.antiPLA2rNegative[idSnp[5],] %in% c("1/1"))]),
								names(GTs.antiPLA2rNegative[idSnp[6],])[which(GTs.antiPLA2rNegative[idSnp[6],] %in% c("1/1"))])
		}
		results_combo[rowno,1] = length(idSnp)
		results_combo[rowno,2] = paste(data$RSID[idSnp], collapse="_")
		results_combo[rowno,3] = paste0(out, collapse="; ")
		results_combo[rowno,4] = length(out)
		results_combo[rowno,5] = paste0(out_n, collapse="; ")
		results_combo[rowno,6] = length(out_n)
	}
	x=0
}

# 11 in total pos., 18 in total neg.
results_combo$pla2r_pos_freq = results_combo$no_overlap_positive/11
results_combo$pla2r_neg_freq = results_combo$no_overlap_negative/18
results_combo$max_diff = results_combo$pla2r_pos_freq - results_combo$pla2r_neg_freq


# Controls (healthy) - with T1D cases included
data <- read.table("data.vcf",	comment.char = "",skip=30,header=T)

newnames <- gsub("^.{0,224}", "", names(data[,c(10:45)]))
data <- data[,c(1:5,10:45)]
names(data) <- c("chr","pos","RSID","ref","alt",newnames)

GTs$no.hom.0 <-  rowSums(GTs[,c(6:ncol(GTs))] == "0/0")
GTs$no.het.01 <- rowSums(GTs[,c(6:ncol(GTs))] == "0/1")
GTs$no.het.10 <- rowSums(GTs[,c(6:ncol(GTs))] == "1/0")
GTs$no.het <- GTs$no.het.01+GTs$no.het.10
GTs$no.hom.1 <- rowSums(GTs[,c(6:ncol(GTs))] == "1/1")
GTs$no.missing <- rowSums(GTs[,c(6:ncol(GTs))] == "./.") 
GTs$total <- GTs$no.hom.0 + GTs$no.het + GTs$no.hom.1 + GTs$no.missing

GTs$alt.af <- (GTs$no.het + (GTs$no.hom.1*2))/(GTs$total*2)
GTs$alt.ac <- GTs$no.het + (GTs$no.hom.1*2)
GTs$an <- GTs$total*2

results_combo_ctrl <- data.frame(no_SNPs = numeric(), SNPs = character(), overlap = numeric(), no_overlap = numeric(), stringsAsFactors=FALSE)

for(m in 1:6){
	print(paste0("combo length: ",m))
	combos <- combinations(n=16, r=m, v=c(1:16), repeats.allowed = FALSE)
	# n Size of the source vector, r=Size of the target vectors, v=Source vector. Defaults to 1:n
	for(x in 1:nrow(combos)){
		print(paste0("x: ",x))
		rowno = nrow(results_combo_ctrl)+1
		idSnp = combos[x,]
		if(m==1){
			out = names(GTs[idSnp[1],])[which(GTs[idSnp[1],] %in% c("1/1"))]
			
		}else if(m==2){
			out = intersect(names(GTs[idSnp[1],])[which(GTs[idSnp[1],] %in% c("1/1"))],
				names(GTs[idSnp[2],])[which(GTs[idSnp[2],] %in% c("1/1"))])
		}else if(m==3){
			out = intersect(intersect(names(GTs[idSnp[1],])[which(GTs[idSnp[1],] %in% c("1/1"))],
				names(GTs[idSnp[2],])[which(GTs[idSnp[2],] %in% c("1/1"))]),
					names(GTs[idSnp[3],])[which(GTs[idSnp[3],] %in% c("1/1"))])
			
		}else if(m==4)		{
			out = intersect(intersect(intersect(names(GTs[idSnp[1],])[which(GTs[idSnp[1],] %in% c("1/1"))],
				names(GTs[idSnp[2],])[which(GTs[idSnp[2],] %in% c("1/1"))]),
					names(GTs[idSnp[3],])[which(GTs[idSnp[3],] %in% c("1/1"))]),
						names(GTs[idSnp[4],])[which(GTs[idSnp[4],] %in% c("1/1"))])
			
		}else if (m==5){
			out = intersect(intersect(intersect(intersect(names(GTs[idSnp[1],])[which(GTs[idSnp[1],] %in% c("1/1"))],
				names(GTs[idSnp[2],])[which(GTs[idSnp[2],] %in% c("1/1"))]),
					names(GTs[idSnp[3],])[which(GTs[idSnp[3],] %in% c("1/1"))]),
						names(GTs[idSnp[4],])[which(GTs[idSnp[4],] %in% c("1/1"))]),
							names(GTs[idSnp[5],])[which(GTs[idSnp[5],] %in% c("1/1"))])
				
		}else if(m==6){
			out = intersect(intersect(intersect(intersect(intersect(names(GTs[idSnp[1],])[which(GTs[idSnp[1],] %in% c("1/1"))],
				names(GTs[idSnp[2],])[which(GTs[idSnp[2],] %in% c("1/1"))]),
					names(GTs[idSnp[3],])[which(GTs[idSnp[3],] %in% c("1/1"))]),
						names(GTs[idSnp[4],])[which(GTs[idSnp[4],] %in% c("1/1"))]),
							names(GTs[idSnp[5],])[which(GTs[idSnp[5],] %in% c("1/1"))]),
								names(GTs[idSnp[6],])[which(GTs[idSnp[6],] %in% c("1/1"))])
			
		}
		results_combo_ctrl[rowno,1] = length(idSnp)
		results_combo_ctrl[rowno,2] = paste(data$RSID[idSnp], collapse="_")
		results_combo_ctrl[rowno,3] = paste0(out, collapse="; ")
		results_combo_ctrl[rowno,4] = length(out)

	}
	x=0
}

results_combo_ctrl$freq = results_combo_ctrl$no_overlap/unique(GTs$total)
nrow(results_combo_ctrl)
# 14892
names(results_combo_ctrl) <- c("no_SNPs","SNPs","overlap_ctrls","no_overlap_ctrls","freq_ctrls")
mn_pla2r_vs_ctrl <- merge(results_combo, results_combo_ctrl, by = c("no_SNPs","SNPs"))
mn_pla2r_vs_ctrl$diff_pla2rPos_ctrl_freq = mn_pla2r_vs_ctrl$pla2r_pos_freq - mn_pla2r_vs_ctrl$freq

# =========== MERGE:
names(mn_pla2r_vs_ctrl) <- c("no_SNPs","SNPs","overlap_positive", "no_overlap_positive","overlap_negative","no_overlap_negative", "pla2r_pos_freq","pla2r_neg_freq","max_diff",
	"overlap_ctrl","freq_ctrl","diff_pla2rPos_ctrl_freq")
names(results_combo_nonMembr) <- c("no_SNPs","SNPs","overlap_nonMN","no_overlap_nonMN","freq_nonMN")
mn_pla2r_vs_ctrl.nonMN <- merge(mn_pla2r_vs_ctrl, results_combo_nonMembr, by = c("no_SNPs","SNPs"))
mn_pla2r_vs_ctrl.nonMN.MN <- merge(mn_pla2r_vs_ctrl.nonMN, results_combo_Membr, by = c("no_SNPs","SNPs"))
mn_pla2r_vs_ctrl.nonMN.MN.nonT1DN <- merge(mn_pla2r_vs_ctrl.nonMN.MN, results_combo_nonMembr_nonT1DN, by = c("no_SNPs","SNPs"))
#mn_pla2r_vs_ctrl.nonMN.MN[,c(2,7,8,12,17,20)]

#============

# Healthy controls - no T1D
data <- read.table("data.vcf", comment.char = "",skip=30,header=T)

newnames <- gsub("^.{0,224}", "", names(data[,c(10:45)]))
data <- data[,c(1:5,10:45)]
names(data) <- c("chr","pos","RSID","ref","alt",newnames)

# Remove T1Ds:
t1d_hc <- read.csv("controls_wdiabetes.csv", header=T)
names(t1d_hc) <- c("PPID","comorbidities","has_diabetes_1","has_diabetes_2")
hc <- read.csv("Control_GSA_samples.csv", header=T)
t1d_hc.hc <- merge(t1d_hc, hc[,c(3,4)], by = c("PPID"))
t1d_hc.hc.NoT1D <- t1d_hc.hc[t1d_hc.hc$has_diabetes_1 %in% FALSE,]
GTs <- data[,colnames(data) %in% c(colnames(data[,c(1:5)]), t1d_hc.hc.NoT1D$SampleID)]
#################

GTs$no.hom.0 <-  rowSums(GTs[,c(6:ncol(GTs))] == "0/0")
GTs$no.het.01 <- rowSums(GTs[,c(6:ncol(GTs))] == "0/1")
GTs$no.het.10 <- rowSums(GTs[,c(6:ncol(GTs))] == "1/0")
GTs$no.het <- GTs$no.het.01+GTs$no.het.10
GTs$no.hom.1 <- rowSums(GTs[,c(6:ncol(GTs))] == "1/1")
GTs$no.missing <- rowSums(GTs[,c(6:ncol(GTs))] == "./.")
GTs$total <- GTs$no.hom.0 + GTs$no.het + GTs$no.hom.1 + GTs$no.missing

GTs$alt.af <- (GTs$no.het + (GTs$no.hom.1*2))/(GTs$total*2)
GTs$alt.ac <- GTs$no.het + (GTs$no.hom.1*2)
GTs$an <- GTs$total*2

# Find SNP combination with max overlap?
#combos <- sapply(2:6, function(x) combinations(n=6, r=x, v=c(snp1, snp2, snp3, snp4, snp5, snp6), repeats.allowed = FALSE))
het_or_hom <- function(x){
	zyg = ifelse(x<17, "het", "hom")
}

# ctrl SNP order: rs1003878, rs1059535, rs1064994, rs10885, rs2233974, rs2395231, rs28630740, rs3130618, rs3132450, rs3134900, rs35771982, rs3749117, rs66484345, rs79350521, rs9260157, rs9268142 \
# pla2rpos SNP  : rs1003878, rs1059535, rs1064994, rs10885, rs2233974, rs2395231, rs28630740, rs3130618, rs3132450, rs3134900, rs35771982, rs3749117, rs66484345, rs79350521, rs9260157, rs9268142 \


results_combo_ctrl_v3 <- data.frame(no_SNPs = numeric(), SNPs = character(), HetHom=character(), overlap_ctrlExclT1D = numeric(), no_overlap_ctrlExclT1D = numeric(), 
overlap_PLA2RposMN = numeric(), no_overlap_PLA2RposMN = numeric(), stringsAsFactors=FALSE)

#for(m in 1:6){
m=6
	print(paste0("combo length: ",m))
	#combos <- combinations(n=32, r=m, v=c(1:32), repeats.allowed = FALSE)
	# n Size of the source vector, r=Size of the target vectors, v=Source vector. Defaults to 1:n
	for(x in 1:nrow(combos)){
		print(paste0("x:",x, " and m:", m))
		rowno = nrow(results_combo_ctrl_v3)+1
		idSnp = combos[x,]
		idSnp1 = idSnp[idSnp<17]
		idSnp2 = idSnp[idSnp>=17]
		idSnp2 = idSnp2 - 16
		if(m==1){		
			if(length(idSnp1)==1){
				out = names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))]
				out2 = names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))]
			}else if(length(idSnp2)==1){
				out = names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))]
				out2 = names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))]
			}
			
		}else if(m==2){
			if(length(idSnp2)==2){
				out = intersect(names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))],
					names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))])
				out2 = intersect(names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))],
					names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))])
			}else if(length(idSnp1)==2){
				out = intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))])
				out2 = intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))])
			} else{
				out = intersect(names(GTs[idSnp1,])[which(GTs[idSnp1,] %in% c("1/0","0/1"))],
					names(GTs[idSnp2,])[which(GTs[idSnp2,] %in% c("1/1"))])
				out2 = intersect(names(GTs.antiPLA2rPositive[idSnp1,])[which(GTs.antiPLA2rPositive[idSnp1,] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp2,])[which(GTs.antiPLA2rPositive[idSnp2,] %in% c("1/1"))])
			}
			
		}else if(m==3){

			if(length(idSnp1)==2){
				out = intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))])
				out2 = intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))])
			} else if(length(idSnp2)==2){
				out = intersect(intersect(names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))],
					names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
						names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))])
				out2 = intersect(intersect(names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))],
					names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))])
			} else if(length(idSnp1)==3){
				out = intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))])
				out2 = intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))])
			} else if(length(idSnp2)==3){
				out = intersect(intersect(names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))],
					names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
						names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))])
				out2 = intersect(intersect(names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))],
					names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))])
			}
			
		}else if(m==4){
			if(length(idSnp1)==2){
				out = intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))]),
							names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))])
			} else if(length(idSnp1)==4){
				out = intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs[idSnp1[4],])[which(GTs[idSnp1[4],] %in% c("1/0","0/1"))])
				out2 = intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs.antiPLA2rPositive[idSnp1[4],])[which(GTs.antiPLA2rPositive[idSnp1[4],] %in% c("1/0","0/1"))])
			} else if(length(idSnp1)==3){
				out = intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))])
			} else if(length(idSnp2)==3){
				out = intersect(intersect(intersect(names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))],
					names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
						names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))]),
							names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))])
				out2 = intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))],
					names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))])
			} else if(length(idSnp2)==4){
				out = intersect(intersect(intersect(names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))],
					names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
						names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))]),
							names(GTs[idSnp2[4],])[which(GTs[idSnp2[4],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))],
					names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[4],])[which(GTs.antiPLA2rPositive[idSnp2[4],] %in% c("1/1"))])
			}
			
		}else if (m==5){
			if(length(idSnp1)==1){
				out = intersect(intersect(intersect(intersect(names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))],
					names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
						names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))]),
							names(GTs[idSnp2[4],])[which(GTs[idSnp2[4],] %in% c("1/1"))]),
								names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))])
				out2 = intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))],
					names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[4],])[which(GTs.antiPLA2rPositive[idSnp2[4],] %in% c("1/1"))]),
								names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))])
			} else if(length(idSnp1)==2){
				out = intersect(intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))]),
							names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
								names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
								names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))])
			} else if(length(idSnp1)==3){
				out = intersect(intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))]),
								names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))]),
								names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))])
			} else if(length(idSnp1)==4){
				out = intersect(intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs[idSnp1[4],])[which(GTs[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs.antiPLA2rPositive[idSnp1[4],])[which(GTs.antiPLA2rPositive[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))])
			} else if(length(idSnp1)==5){
				out = intersect(intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs[idSnp1[4],])[which(GTs[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs[idSnp1[5],])[which(GTs[idSnp1[5],] %in% c("1/0","0/1"))])
				out2 = intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs.antiPLA2rPositive[idSnp1[4],])[which(GTs.antiPLA2rPositive[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs.antiPLA2rPositive[idSnp1[5],])[which(GTs.antiPLA2rPositive[idSnp1[5],] %in% c("1/0","0/1"))])
			} else if(length(idSnp2)==5){
				out = intersect(intersect(intersect(intersect(names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))],
					names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
						names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))]),
							names(GTs[idSnp2[4],])[which(GTs[idSnp2[4],] %in% c("1/1"))]),
								names(GTs[idSnp2[5],])[which(GTs[idSnp2[5],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))],
					names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[4],])[which(GTs.antiPLA2rPositive[idSnp2[4],] %in% c("1/1"))]),
								names(GTs.antiPLA2rPositive[idSnp2[5],])[which(GTs.antiPLA2rPositive[idSnp2[5],] %in% c("1/1"))])
			} 
				
		}else if(m==6){
			if(length(idSnp1)==1){
				out = intersect(intersect(intersect(intersect(intersect(names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))],
					names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
						names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))]),
							names(GTs[idSnp2[4],])[which(GTs[idSnp2[4],] %in% c("1/1"))]),
								names(GTs[idSnp2[5],])[which(GTs[idSnp2[5],] %in% c("1/1"))]),
									names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))])
				out2 = intersect(intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))],
					names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[4],])[which(GTs.antiPLA2rPositive[idSnp2[4],] %in% c("1/1"))]),
								names(GTs.antiPLA2rPositive[idSnp2[5],])[which(GTs.antiPLA2rPositive[idSnp2[5],] %in% c("1/1"))]),
									names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))])
			} else if(length(idSnp1)==2){
				out = intersect(intersect(intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))]),
							names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
								names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))]),
									names(GTs[idSnp2[4],])[which(GTs[idSnp2[4],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
								names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))]),
									names(GTs.antiPLA2rPositive[idSnp2[4],])[which(GTs.antiPLA2rPositive[idSnp2[4],] %in% c("1/1"))])
			} else if(length(idSnp1)==3){
				out = intersect(intersect(intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))]),
								names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
									names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))]),
								names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
									names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))])
			} else if(length(idSnp1)==4){
				out = intersect(intersect(intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs[idSnp1[4],])[which(GTs[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))]),
									names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs.antiPLA2rPositive[idSnp1[4],])[which(GTs.antiPLA2rPositive[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))]),
									names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))])
			} else if(length(idSnp1)==5){
				out = intersect(intersect(intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs[idSnp1[4],])[which(GTs[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs[idSnp1[5],])[which(GTs[idSnp1[5],] %in% c("1/0","0/1"))]),
									names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs.antiPLA2rPositive[idSnp1[4],])[which(GTs.antiPLA2rPositive[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs.antiPLA2rPositive[idSnp1[5],])[which(GTs.antiPLA2rPositive[idSnp1[5],] %in% c("1/0","0/1"))]),
									names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))])
			} else if(length(idSnp1)==6){
				out = intersect(intersect(intersect(intersect(intersect(names(GTs[idSnp1[1],])[which(GTs[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs[idSnp1[2],])[which(GTs[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs[idSnp1[3],])[which(GTs[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs[idSnp1[4],])[which(GTs[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs[idSnp1[5],])[which(GTs[idSnp1[5],] %in% c("1/0","0/1"))]),
									names(GTs[idSnp1[6],])[which(GTs[idSnp1[6],] %in% c("1/0","0/1"))])
				out2 = intersect(intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp1[1],])[which(GTs.antiPLA2rPositive[idSnp1[1],] %in% c("1/0","0/1"))],
					names(GTs.antiPLA2rPositive[idSnp1[2],])[which(GTs.antiPLA2rPositive[idSnp1[2],] %in% c("1/0","0/1"))]),
						names(GTs.antiPLA2rPositive[idSnp1[3],])[which(GTs.antiPLA2rPositive[idSnp1[3],] %in% c("1/0","0/1"))]),
							names(GTs.antiPLA2rPositive[idSnp1[4],])[which(GTs.antiPLA2rPositive[idSnp1[4],] %in% c("1/0","0/1"))]),
								names(GTs.antiPLA2rPositive[idSnp1[5],])[which(GTs.antiPLA2rPositive[idSnp1[5],] %in% c("1/0","0/1"))]),
									names(GTs.antiPLA2rPositive[idSnp1[6],])[which(GTs.antiPLA2rPositive[idSnp1[6],] %in% c("1/0","0/1"))])
			} else if(length(idSnp2)==6){
				out = intersect(intersect(intersect(intersect(intersect(names(GTs[idSnp2[1],])[which(GTs[idSnp2[1],] %in% c("1/1"))],
					names(GTs[idSnp2[2],])[which(GTs[idSnp2[2],] %in% c("1/1"))]),
						names(GTs[idSnp2[3],])[which(GTs[idSnp2[3],] %in% c("1/1"))]),
							names(GTs[idSnp2[4],])[which(GTs[idSnp2[4],] %in% c("1/1"))]),
								names(GTs[idSnp2[5],])[which(GTs[idSnp2[5],] %in% c("1/1"))]),
									names(GTs[idSnp2[6],])[which(GTs[idSnp2[6],] %in% c("1/1"))])
				out2 = intersect(intersect(intersect(intersect(intersect(names(GTs.antiPLA2rPositive[idSnp2[1],])[which(GTs.antiPLA2rPositive[idSnp2[1],] %in% c("1/1"))],
					names(GTs.antiPLA2rPositive[idSnp2[2],])[which(GTs.antiPLA2rPositive[idSnp2[2],] %in% c("1/1"))]),
						names(GTs.antiPLA2rPositive[idSnp2[3],])[which(GTs.antiPLA2rPositive[idSnp2[3],] %in% c("1/1"))]),
							names(GTs.antiPLA2rPositive[idSnp2[4],])[which(GTs.antiPLA2rPositive[idSnp2[4],] %in% c("1/1"))]),
								names(GTs.antiPLA2rPositive[idSnp2[5],])[which(GTs.antiPLA2rPositive[idSnp2[5],] %in% c("1/1"))]),
									names(GTs.antiPLA2rPositive[idSnp2[6],])[which(GTs.antiPLA2rPositive[idSnp2[6],] %in% c("1/1"))])
			}	
		}

		results_combo_ctrl_v3[rowno,1] = length(idSnp)
		results_combo_ctrl_v3[rowno,2] = paste0(c(data$RSID[idSnp1], data$RSID[idSnp2]), collapse="_")
		results_combo_ctrl_v3[rowno,3] = paste0(sapply(idSnp, het_or_hom), collapse="_")
		results_combo_ctrl_v3[rowno,4] = paste0(out, collapse="; ")
		results_combo_ctrl_v3[rowno,5] = length(out)
		results_combo_ctrl_v3[rowno,6] = paste0(out2, collapse="; ")
		results_combo_ctrl_v3[rowno,7] = length(out2)

	}
	x=0
}

write.csv(results_combo_ctrl_v3, "AF_combos_5lead_plus_PLA2RLead_and_AllPredmissense_16total_allGroups_incnonT1D_Hom_Het_v3.csv", quote=F, row.names=F)

results_combo_ctrl_v3$freq_ctrlExclT1D = results_combo_ctrl_v3$no_overlap_ctrlExclT1D/unique(GTs$total)
results_combo_ctrl_v3$freq_PLA2RposMN = results_combo_ctrl_v3$no_overlap_PLA2RposMN/unique(GTs.antiPLA2rPositive$total)
results_combo_ctrl_v3$diff_freq_PLA2RposMN = abs(results_combo_ctrl_v3$freq_ctrlExclT1D - results_combo_ctrl_v3$freq_PLA2RposMN)
head(results_combo_ctrl_v3[order(-results_combo_ctrl_v3$diff_freq_PLA2RposMN),])

results_combo_ctrl_v3.diff_freq_0.3 <- results_combo_ctrl_v3[results_combo_ctrl_v3$diff_freq_PLA2RposMN>=0.3,]
write.csv(results_combo_ctrl_v3.diff_freq_0.3, "AF_combos_5lead_plus_PLA2RLead_and_AllPredmissense_16total_allGroups_incnonT1D_Hom_Het_v3_0.3.csv", quote=F, 
	row.names=F, col.names=T)


# check combos that contain PLA2R vars
results_combo_ctrl_v3.diff_freq_0.3[grep("rs35771982",results_combo_ctrl_v3.diff_freq_0.3$SNPs),]
# only 1 missense PLA2R or both PLA2R
results_combo_ctrl_v3.diff_freq_0.3[grep("rs3749117",results_combo_ctrl_v3.diff_freq_0.3$SNPs),]
# rs66484345: pla2r1 lead snp
results_combo_ctrl_v3.diff_freq_0.3[grep("rs66484345",results_combo_ctrl_v3.diff_freq_0.3$SNPs),]
# 0
# hla-dqa1 utr3 lead snp:
results_combo_ctrl_v3.diff_freq_0.3[grep("rs1064994",results_combo_ctrl_v3.diff_freq_0.3$SNPs),]
# 0

table(results_combo_ctrl_v3[results_combo_ctrl_v3$diff_freq_PLA2RposMN>=0.3,3])

nrow(results_combo_ctrl_v2)
# 1149016
names(results_combo_ctrl_v2) <- c("no_SNPs", "SNPs", "overlap_ctrls_exclT1D", "no_overlap_ctrls_exclT1D", "freq_ctrls_exclT1D")

mn_pla2r_vs_ctrl.nonMN.MN.nonT1DN.ctrl_NonT1D <- merge(mn_pla2r_vs_ctrl.nonMN.MN.nonT1DN, results_combo_ctrl, by = c("no_SNPs","SNPs"),all.x=T)

write.csv(results_combo_ctrl_v2, "AF_combos_5lead_plus_PLA2RLead_and_AllPredmissense_16total_allGroups_incnonT1D_Hom_Het.csv", quote=F, row.names=F)
write.csv(mn_pla2r_vs_ctrl.nonMN.MN.nonT1DN.ctrl_NonT1D, "AF_combos_5lead_plus_PLA2RLead_and_AllPredmissense_16total_allGroups_incnonT1D_Hom_Het.csv", quote=F, row.names=F)


