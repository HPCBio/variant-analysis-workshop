# This script generates the example dataset for this workshop.

library(VariantAnnotation)

# VCFs from http://cbsusrv04.tc.cornell.edu/users/panzea/download.aspx?filegroupid=34
vcfpath <- "/mnt/lvclark1/maize_snp/" # point to where you have the VCF files
vcffiles <- paste0("hmp321_agpv4_chr", 1:10, ".vcf.gz")

myfile <- VcfFile(paste0(vcfpath, vcffiles[1]))

mysam <- samples(scanVcfHeader(myfile))

# Just import 200 samples
set.seed(1110)
subsam <- mysam[sort(sample(length(mysam), 200))]

# Get one portion of one chromosome
myrange <- GRanges("1", IRanges(start = 21e6, end = 24e6))

# read this in
svp <- ScanVcfParam(samples = subsam, which = myrange)
myvcf <- readVcf(myfile, genome = "1", param = svp)

myvcf # 120,304 markers
rowRanges(myvcf)
info(myvcf) # LLD indicates high quality markers according to Panzea.org

writeVcf(myvcf, filename = "results/hmp321_agpv4_chr1_subset.vcf",
         index = TRUE)

test <- readVcf("results/hmp321_agpv4_chr1_subset.vcf.bgz")
