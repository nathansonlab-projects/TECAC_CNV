setwd("~/Documents/nathansonlab/CNV/")
library(data.table)
library(IRanges)
library(stringr)
library(ggplot2)

# ============================== constants ============================ #
# reference data
grch37.ref.dat = data.frame( chromosome = c(seq(1:22), "X", "Y"),
                             centromere.start = c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166,
                                                  58054331, 43838887, 47367679, 39254935, 51644205, 34856694,
                                                  16000000, 16000000, 17000000, 35335801, 22263006, 15460898,
                                                  24681782, 26369569, 11288129, 13000000, 58632012, 10104553),
                             centromere.end = c(124535434, 95326171, 93504854, 52660117, 49405641, 61830166,
                                                61054331, 46838887, 50367679, 42254935, 54644205, 37856694,
                                                19000000, 19000000, 20000000, 38335801, 25263006, 18460898,
                                                27681782, 29369569, 14288129, 16000000, 61632012, 13104553),
                             p.telomere.end = rep(10000, 24),
                             q.telomere.start = c(249240621, 243189373, 198012430, 191144276, 180905260,
                                                  171105067, 159128663, 146354022, 141203431, 135524747,
                                                  134996516, 133841895, 115159878, 107339540, 102521392,
                                                  90344753, 81185210, 78067248, 59118983, 63015520,
                                                  48119895, 51294566, 155260560, 59363566),
                             chr.size = c(249250621, 243199373,  198022430, 191154276, 180915260,
                                          171115067, 159138663, 146364022, 141213431, 135534747,
                                          135006516, 133851895, 115169878, 107349540, 102531392,
                                          90354753, 81195210, 78077248, 59128983, 63025520,
                                          48129895, 51304566, 155270560, 59373566)

)
# ================================================================== #


# ======================= functions ================================= #
# ------------------------- in.telomere ---------------------------- #
in.telomere <- function( start.pos, end.pos, q.tel.start, OFFSET )
# input: start.pos (integer), the starting position of the segment
# end.pos (integer), the ending position of the segment
# q.tel.start (integer), the starting position of the q-arm telomere (from ref data)
# OFFSET (integer), the offset in bp. exclude all snps with in OFFSET of the telomere.
#
# output: boolean value
{
  return(start.pos <= 10000 + OFFSET | end.pos >= q.tel.start - OFFSET)
}
# ------------------------------------------------------------------- #

# --------------------------- in.centromere ------------------------- #
in.centromere <- function (start.pos, end.pos, ref.cent.start, ref.cent.end)
{
  return(start.pos < ref.cent.start & end.pos > ref.cent.end |
           start.pos > ref.cent.start & start.pos < ref.cent.end |
           end.pos > ref.cent.start & end.pos < ref.cent.end)
}
# ------------------------------------------------------------------- #



# ------------------------- inRegion -------------------------------- #
# find regions of B overlapping A
inRegion <- function( region.start, region.end, rangeA )
{
  rangeB <- IRanges(region.start, region.end)
  return( findOverlaps( rangeA, rangeB, type = "within"))
}
# ------------------------------------------------------------------- #



# ------------------------------------------------------------------- #
fixName <- function( id )
{
  return( paste(strsplit(id, "_")[[1]][1:2], collapse = "_"))
}
# ------------------------------------------------------------------- #



# ------------------------------------------------------------------- #
CNV.assocTest <- function( x, dat, min.obs = 50 )
  # run association testing
  # input: x, (integer), the vector of CNV values
{
  if( class(x) != "matrix")
  {
    x <- as.matrix(x)
  }


  if( length(unique(x)) < 2)
  {
    return(NA)
  } else
    if( sum( x > 0 ) < min.obs  )
    {
      return(NA)
    }

  fit <- glm(data = dat, phenotype ~ x + EV1 + EV2 + EV3 + as.factor(Site), family = "binomial")
  fit <- glm(data = dat, phenotype ~ x + as.factor(Site), family = "binomial")
  return(summary(fit)$coefficients[2,4])
}
# ------------------------------------------------------------------- #


# =================================================================== #





# ================== MAIN ========================================== #
# for now, Y chromosome wont be considered, drop it
grch37.ref.dat <- grch37.ref.dat[1:22,]

qc.dat <- read.table("exomehmm.qcsum", header = T)
qc.dat$ID <- unlist(lapply(strsplit(as.character(qc.dat$File), "\\."), function(x) x[2]))

# this contains the BID - SSID mapping
pheno <- read.table("TECAC_CNV_PHENO.csv", header = TRUE, sep = ",")
qc.dat$pheno <- pheno$PHENO[match(qc.dat$ID, pheno$BID)]

# read in the covariate data and match to phenotype data
cov <- read.table("case_ctl_cauc_3ev.sample", header = TRUE, as.is = TRUE)
cov <- cov[-1,]
ind <- unlist(substr(cov$ID_2, 1, 2) == "NO" & str_count(cov$ID_2, "_") == 3)
cov$ID <- "0"
cov$ID[!ind] <- unlist(lapply(strsplit(cov$ID_2[!ind], "_"), function(x) x[1]))
cov <- cov[ substr(cov$ID, 1, 3) != "MDA", ]

#cov$ID[ind] <- unlist(lapply(cov$IID[ind], fixName))
cov <- cov[ cov$ID_2 %in% pheno$SSID,]
cov <- cov[match(pheno$SSID, cov$ID_2),]
#cov <- cov[!duplicated(cov$ID_2),]
cov$BID <- pheno$BID[match(pheno$SSID, cov$ID_2)]
cov$phenotype <- as.integer(cov$phenotype)

p1  <- ggplot(data = qc.dat, aes(x = LRR_SD,  y = NumCNV)) + 
  geom_point() +
  theme_minimal() +
  ggtitle("LRR_SD vs NumCNV")
p1

p2  <- ggplot(data = qc.dat, aes(x = WF,  y = NumCNV)) + 
  geom_point() +
  theme_minimal() +
  ggtitle("Wave Factor vs NumCNV")
p2

p3  <- ggplot(data = qc.dat, aes(x = BAF_drift, y = NumCNV)) +
  geom_point() +
  theme_minimal() +
  ggtitle("BAF_drift vs NumCNV")
p3

# remove subjects w/ excess LRR_SD
qc.dat <- qc.dat[ qc.dat$LRR_SD < 0.275, ]

# alternatively, could use 95th percentile
# lrrsd <- qc.dat[,4]
# qc.dat$LRR_SD <- as.numeric(qc.dat$LRR_SD)
# LRR.thresh <- quantile(lrrsd, probs = 0.95)
# qc.dat <- qc.dat[ qc.dat$LRR_SD < LRR.thresh, ]

# this seems too restrictive
# get the 95th percentile value for BAF drift
# exclude samples with BAF in the 95th+ percentile
# qc.dat$BAF_mean <- as.numeric(qc.dat$BAF_mean)
# BAF.thresh <- quantile(qc.dat$BAF_drift, probs = 0.95)
# qc.dat <- qc.dat[ qc.dat$BAF_drift < BAF.thresh, ]
# try using simple outliers
qc.dat  <- qc.dat[ qc.dat$BAF_drift < 0.002,]


# ---
# parse all the extra text from penncnv; this sets up the snp/cnv data
# for further processing
dat <- as.data.frame(fread("clean_bp_exomehmm.rawcnv", header = F))
dat$V2 <- as.integer(gsub("numsnp=", "", dat$V2))
 
dat$V3 <- as.integer( gsub(",", "", gsub("length=", "", dat$V3)))
 
dat$V4 <- as.character(dat$V4)
dat$state <- unlist(lapply(strsplit(dat$V4, ","), function(x) x[1]))
dat$CN <- as.integer(gsub("cn=", "", unlist(lapply(strsplit(dat$V4, ","), function(x) x[2]))))

dat$V5 <- as.character(dat$V5)


dat$start <- as.character(dat$V6)
dat$start <- gsub("startsnp=", "", dat$start)
dat$chr <- unlist(lapply(strsplit(dat$start, ":"), function(x) x[1]))
dat$start <- as.integer(unlist(lapply(strsplit(dat$start, ":"), function(x) x[2])))

dat$end <- as.character(dat$V7)
dat$end <- gsub("endsnp=", "", dat$end)
dat$end <- as.integer(unlist(lapply(strsplit(dat$end, ":"), function(x) x[2])))

colnames(dat) <- c("V1", "nsnp", "seglen", "V4", "ID", "V6", "V7", "conf", "state", "CN", "start", "chr", "end")
dat <- dat[, -grep("V", colnames(dat))]

dat$conf <- as.character(dat$conf)
dat$conf <- gsub("conf=", "", dat$conf)
dat$conf <- as.numeric(dat$conf)

dat$IID <- gsub("split_files/sample.", "", dat$ID)

# ---

n.cnv <- dim(dat)[1]

# ----- qc ----

# truncate based on subject qc
dat <- dat[ dat$ID %in% qc.dat$File, ]

# retain CNVs with >= N snps
dat <- dat[ dat$nsnp >= 5,]

# remove CNVs within OFFSET of centromere/telomere
# offset is 1mb 
dat$in.centromere <- FALSE
dat$in.telomere <- FALSE
OFFSET <- 1000000 # 1mbp by default

# no x data?
# for( i in c(1:22))
# {
#   chr=i
#   
#   ref.dat <- grch37.ref.dat[grch37.ref.dat$chromosome == chr,]
#   
#   dat$in.centromere[dat$chr == chr] <- mapply(in.centromere, dat$start[dat$chr == chr], dat$end[dat$chr == chr],
#                                               ref.dat$centromere.start - OFFSET, ref.dat$centromere.end + OFFSET)
#   dat$in.telomere[dat$chr == chr] <- mapply(in.telomere, dat$start[dat$chr == chr], dat$end[dat$chr == chr],
#                                             ref.dat$q.telomere.start, OFFSET)
# }

# dat <- dat[ !dat$in.centromere & !dat$in.telomere, ]

# remove cnvs in HLA
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
# chr6, 28477797-33448354
dat$inHLA <- dat$start >= 28477797 & dat$end <= 33448354 & dat$chr == 6
dat <- dat[ !dat$inHLA, ]

print("QC complete")
print(paste("Started with: ", n.cnv, " snps", sep = ""))
print(paste(n.cnv - dim(dat)[1], " snps removed."))
print(paste(dim(dat)[1], " snps remaining", sep = ""))

# write.csv(dat,"C:/Users/kxrua/Documents/Lab/QC_CNV/pdat2_clean_bp_exomehmm.csv", row.names = FALSE)

# plot cnvs per chromosome after qc
dat$chr <- as.integer(dat$chr)
p1 <- ggplot(data = dat, aes(sort(chr))) + 
   geom_bar() +
   theme_minimal() +
   scale_x_discrete("chr",breaks =  seq(1:22), labels = as.character(seq(1:22)), limits = seq(1:22))
p1
 
 


cov <- cov[cov$BID %in% qc.dat$ID,]
#cov <- merge(cov, pheno, by.x = "ID", by.y = "SSID", all.x = T)
dat <- dat[ dat$IID %in% qc.dat$ID,]
for (chromosome in 1:22)
{
 
   # convert CNVs to Copy Number Polymorphic Regions (CNPR)
   # need to standardize regions across subjects
   dat1 <- dat[ dat$chr == chromosome,]
 
   # define start and end points for this chromosome
   # start points are the unique set of starting points from all subjects
   # end points are 1bp before the next start point
   start.pts <- sort(unique(dat1$start))
   end.pts <- start.pts[2:(length(start.pts))] - 1
   end.pts <- c(end.pts, grch37.ref.dat$chr.size[grch37.ref.dat$chromosome == i])
   rangeA <- IRanges(start.pts, end.pts)
 
   # subject x CNPR matrix, cells are copy number value
   m <- matrix(0, nrow = length(unique(dat1$ID)), ncol = length(start.pts))
 
   # can be multiple regions per subject, so have to iterate by subject
   # can this be sped with apply?
   for( i in 1:length(unique(dat1$ID)) )
   {
     tmp <- dat1[ dat1$ID == unique(dat1$ID)[i],]
     ind <- c()
 
     # for subject i, make a vector of all regions with non-zero copy number status
     for(j in 1:dim(tmp)[1])
     {
       ind <- c(ind, inRegion(tmp$start[j], tmp$end[j], rangeA )@from)
     }
 
     m[i, ind] <- dat1$CN[ind]
   }
 
   rownames(m) <- unique(dat1$IID)
 
 
   df <- data.frame( start= rangeA@start,
                     end = rangeA@start + rangeA@width,
                     p.value = apply(m, 2, CNV.assocTest, cov[cov$BID %in% row.names(m),], 40))
   df$p.value[is.na(df$p.value)] <- 1
 
 
   # manhattan-esq plot
   p1 <- ggplot() + geom_point(data=df, aes(x = start, y = -log10(p.value), colour = end - start)) +
     theme_minimal()
   #png(paste('clean_exomehmm_p',chromosome,'.png',sep=''))
   print(p1)
  # dev.off()
 
}
