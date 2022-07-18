#!/usr/bin/env Rscript
###parse options

Script.Start = Sys.time()
qual1 = 25
qual2 = 20
SNP_Num = 7
out = "/Users/aggars05/Desktop/Analysis_Tools/output/SDA281/Spn_Haplotypes.csv"
ref <- strsplit('TATATCAAAGAGGGAAAAACGACTTCTGATTTTGAAGTAAGTAATCAAGAAAAATCAGCAGCCACTCCTACAAAACAACTGAAACAATGAAACAAAAACTAATAACGTAACGTGACTGGCAAGAGATATTTTTAAAACAATGAATAGGTTTACACTTACTTTAGTTTTATGGAAATGAAAGATCATATCATATATAATCTANAATAAAATTAACTAAAATAATTATTATCTAGATAAAAAATTTAGAAGCCAATGAAATCTATAAATAAACTAAATTAAGTTTATTTAATTAACAACTATGGATATAAAATAGGTACTAATCAAAATAGT',split='')[[1]]

####build a dictionary for the base quality score
IllumQ.Qscore <- 0:40
names(IllumQ.Qscore) <- c('!','"','#','$','%','&',"'",'(',')','*','+',',','-','.','/',0:9,':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I')

###list of the variable position, numbered from the beginning of the amplicon, with an added offset
pos = c(82, 83, 84, 90, 91, 92, 96) + 54

library(data.table); library(tidyr); library(seqinr); library(dplyr)

###load the sam file and process it to be easier to use
data <- read.csv(file = "/Users/aggars05/Desktop/Analysis_Tools/output/SDA281/joined_iga.sam", sep = '\t', header = FALSE, skip=2)
data <- data[,c(1,2,4:6,8:11)] #remove unnecessary columns
colnames(data)<- c('QNAME','FLAG','POS','MAPQ','CIGAR','MPOS','ISIZE','SEQ','QUAL')

#remove unmatched reads (ie low quality of mapping)
cat('number of unaligned reads:', as.numeric(table(data$MAPQ<50)[2]),'\n')

strsplit(as.character(data[1,"SEQ"]),split = '')[[1]]

data = data %>% 
  mutate(POS = as.numeric(POS),
         LEN = SEQ %>% as.character %>% nchar) %>%
  filter(MAPQ > 50, # must be mapped to reference
         POS < 50, POS + LEN >= 64) # must cover full span of SNPs

data <- tidyr::separate(data, col='CIGAR', into=c('C1','C2'),sep="(?<=[A-Za-z])(?=[0-9])")
data$POS %>% table

#extract the alignment length information in a new column

data$M <- data$C1 #matching positions
data$M[!grepl('M',data$M)] <- data$C2[!grepl('M',data$M)]
data$M <- strtoi(gsub('\\D','',data$M))
data$S <- rep(0,nrow(data)) #soft clipping in the beginning of the read
data$S[grepl('S',data$C1)] <- strtoi(gsub('\\D','',data$C1)[grepl('S',data$C1)])

#add the start position of alignment on the read (ie taking into account soft clipping)
data$RPOS <- data$POS - 1 + data$S

#data = data[1:min(nrow(data), 5000),]

# dum = strsplit(as.character(data[2,'SEQ']),split='')[[1]] #sequence
# Hapl
# matches = rep(0,50)
# Hapl
# dum
# pos
# dum
# for (offset in 1:33) {
# 
#   matches[offset] = ((Hapl[2,] %>% as.character) == (dum[pos - offset])) %>% sum
# }
# matches
# which(matches == 12)
# data$S
### Reference sequence
# with the WT plasmid containing the Xho1 restriction site:
#ref <- strsplit('ATCAATTTGCCCTTGGACAGGGAACAACACTAAACAACAGGCATTCAAATGACACAGTACATGATAGGACCCCTTATCGAACCCTATTGATGCTCGAGTTGGGTGTTCCATTTCATTTGGGAACCAAGCAAGTGTGTATAGCATGGTCCAGCTCAAGTTGTCACGATGGAAAAGCATGGCTGCATGTTTGTGTAACTGGGCATGATGAAAATGCAACTGCTAGCTTCATTTACGATGGGAGACTTGTAGATAGTATTGGTTCATGGTCCAAAAAAATCCTCAGGACCCAGG',split='')[[1]]

Hapl <- t(sapply(unique(data$QNAME), function(ampl){
  
  reads <- which(data$QNAME==ampl) #identify all the mapped reads of the pair
  haplM <- matrix(ncol=SNP_Num,nrow=length(reads)); colnames(haplM)<- pos; rownames(haplM)<-reads #create a matrix listing the variants at each variable position, for each read of the pair
  qsc <- matrix(ncol=SNP_Num,nrow=length(reads)); colnames(qsc)<- pos; rownames(qsc)<-reads #create a matrix of the quality score of each variable position for each mapped read

  for (r in reads){
    
    Seq <- strsplit(as.character(data[r,'SEQ']),split='')[[1]] #sequence
    Qual <- IllumQ.Qscore[strsplit(as.character(data[r,'QUAL']),split='')[[1]]] #translated quality score of the sequence (ie numbers 0-40)
    
    #check that read aligns with ref
    SNPposR <- pos-data[r,'RPOS'] ; SNPposR <- SNPposR[SNPposR>0 & SNPposR<=data$S[r]+1+data$M[r]] #variant positions numbered from the beginning of the read, and present on the read
    R <- ref[data$POS[r] : (data$POS[r]+data$M[r]-1) ] [-SNPposR] [Qual[ (data$S[r]+1) : (data$S[r]+data$M[r]) ] [-SNPposR] > qual2] #ref sequence outside the variant positions where the amplicon sequence has a quality > quality required for two reads variant quality
    Q <- Seq[ (data$S[r]+1) : (data$S[r]+data$M[r]) ] [-SNPposR] [Qual[ (data$S[r]+1) : (data$S[r]+data$M[r]) ] [-SNPposR] > qual2] #amplicon sequence outside the variant positions where the amplicon sequence has a quality > quality required for two reads variant quality

    #identify the variants on the variable positions
    for(p in pos){
      if(data[r,'M']+data[r,'RPOS']>p & data[r,'RPOS']<p){ #if the read contains the position of interest
        haplM[as.character(r),as.character(p)] <- Seq[p-data[r,'RPOS']] #Identify the variant
        qsc[as.character(r),as.character(p)]<- Qual[p-data[r,'RPOS']] #Identify its quality
      }
    }
  }
  
  haplM[qsc<qual2]<-NA ; qsc[qsc<qual2]<-NA #remove the variantys identified if the read quality if < quality required for two reads variant quality
  hapl <- rep(NA,SNP_Num) #the final aggregated variantys value
  names(hapl) <- pos

  #Control the quality of the identified variants and aggregate them for each pair
  for(p in pos){
    if(length(unique(haplM[,as.character(p)][!is.na(haplM[,as.character(p)])]))==1 ) { #if the variants identified for the two reads of the pair are matching
      if( length(haplM[,as.character(p)][!is.na(haplM[,as.character(p)])]) == 2 ) {hapl[as.character(p)] <- unique(haplM[,as.character(p)][!is.na(haplM[,as.character(p)])]) #if the two paired reads map the region
      } else if( length(haplM[,as.character(p)][!is.na(haplM[,as.character(p)])]) == 1 & qsc[,as.character(p)][!is.na(qsc[,as.character(p)])]>qual1 )  hapl[as.character(p)] <- unique(haplM[,as.character(p)][!is.na(haplM[,as.character(p)])]) #if only one of the two paired reads map the region, filter quality>32
    }
  }

  hapl
  }))

Hapl2 = (Hapl)[!sign(is.na((Hapl)) %>% rowSums),]

for (i in c(3,6)) {

  Hapl2 = Hapl2[which(Hapl2[,i] %in% c("C","A")),]

}

cat('total number of uncomplete haplotypes:',nrow(Hapl) - nrow((Hapl2)),'\n')

Hapl3 = apply(Hapl2, MARGIN = 1,function(i) paste(i,collapse=''))

Hapl4 = cbind(table(Hapl3) %>% names,
      table(Hapl3) %>% as.numeric)

colnames(Hapl4) = c("Haplotype","Count")

write.csv(Hapl4, file = out, row.names = FALSE)

cat('total number of unique haplotypes:',nrow(Hapl4),'\n')

Script.End = Sys.time()
cat('R .SAM Analysis Time:', Script.End - Script.Start)
Script.End - Script.Start
Hapl4[,2] %>% as.numeric %>% hist