#Run Bambu and record execution time _ dataset= LSK114_sequins_sub150k

test.bam <- "LSK114_chrIS_mixA_cDNA_sub150k.bam"
fa.file <- "chrIS.fa"
gtf.file <- "D:/Anaquin-main/data/transcriptome/rnasequin_annotation_2.4.gtf"

start.time1 <- Sys.time()
bambuAnnotations <- prepareAnnotations(gtf.file)
end.time1 <- Sys.time()

start.time2 <- Sys.time()
se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file)
end.time2 <- Sys.time()

start.time3 <- Sys.time()
se_noRef <- bambu(reads = test.bam, genome = fa.file, NDR=1)
end.time3 <- Sys.time()

start.time4 <- Sys.time()
writeBambuOutput(se, path = "test/")
end.time4 <- Sys.time()

start.time5 <- Sys.time()
writeBambuOutput(se_noRef, path = "test_noRef/")
end.time5 <- Sys.time()

time.taken1 <- round(end.time1 - start.time1,2)
time.taken2 <- round(end.time2 - start.time2,2)
time.taken3 <- round(end.time3 - start.time3,2)
time.taken4 <- round(end.time4 - start.time4,2)
time.taken5 <- round(end.time5 - start.time5,2)

time.taken1 #Time to prepare annotations (pre-processing)
time.taken2 #Time to run Bambu
time.taken3 #Time to run Bambu-noRef
time.taken4 #Time to write Bambu output (post-processing)
time.taken5 #Time to write Bambu-noRef output (post-processing)
