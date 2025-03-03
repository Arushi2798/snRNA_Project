# Assuming 'missing_barcodes' is a vector of the missing barcodes
# and 'new_sample_ids' is a data frame containing the new SampleIDs

# Create a named vector or a data frame that maps missing barcodes to new SampleID
new_sample_ids <- data.frame(Barcode = c("AACAAAGGTGTTGAGG-2","AAGGAATCAATTGGTC-2"  ,"ACCTGAAGTTGGTGTT-2" , "AGGGTGAGTATGTCAC-2"  ,"ATCCTATTCGACATAC-2" ,
                                         "ATGTCTTGTCACAGAG-2" , "CAGCCAGTCCTCAGAA-2",  "CTGAATGTCCCTCTCC-2",  "GACCGTGTCAAATAGG-2" , "GTGACGCGTGGTTCTA-2", 
                                         "GTTCATTCATCCAATG-2"  ,"GTTGAACCATGACGAG-2",  "TATCCTAAGCCTTTGA-2",  "AACCTTTGTGCCTAAT-3",  "ACTCTCGGTATCGTGT-3" ,
                                         "AGAAGCGGTGCTTCAA-3"  ,"AGCTCAATCGTAACTG-3",  "CATAAGCCATACAGCT-3",  "CTCAGTCGTCAGTTTG-3" , "CTCGAGGTCGGAAGGT-3" ,
                                         "GACCGTGGTCCCGTGA-3" ,"GAGTCATTCCCAGGCA-3" , "GATTTCTCACAACGTT-3",  "GTCAAGTTCTTTGGAG-3" , "TCGATTTAGACGGTCA-3" ,
                                         "TCTGTCGAGGTTCTAC-3"  ,"TCTTAGTAGACAGCGT-3",  "TGCAGATGTCAGATTC-3",  "TGCGATAGTTTCCAAG-3" , "TGGGAGATCTTTCTAG-3" ,
                                         "TTCATTGAGAATGTTG-3" , "AAACGCTAGATGAATC-6",  "AATGGCTGTCGCAACC-6" , "AGTGCCGCACGGCACT-6",  "ATCACGACATAGATCC-6" ,
                                         "ATCTTCATCCGTAGTA-6" , "CCAATTTGTATATGGA-6",  "CCTGTTGTCTTCACGC-6" , "CGAGTTAAGTGTTCCA-6" , "CGGACACTCCTCTGCA-6" ,
                                         "GACAGCCTCCGCACGA-6" , "GACATCAGTTTGTGGT-6",  "GCCCAGAGTTCGGTTA-6" , "GTTAGACCACACCTTC-6"  ,"TAGACTGGTGCATACT-6" ,
                                         "TATCCTAGTTTAAGGA-6" , "TCAGGTACAGAGCCCT-6",  "TGACAGTTCGTGGAAG-6" , "TGCATGAGTCGTATTG-6",  "CTGCTCATCAAGGTGG-8" ,
                                         "GAGTCATCAGGCAATG-8" , "GCCATGGAGGTCGAGT-8",  "GCCGTGATCTCTGACC-8" , "TAGGAGGAGCAGCCTC-8" , "TCTTAGTGTCCACTCT-8" ,
                                         "TGACAGTGTTGTGTAC-8" , "TGGGCGTTCTGGTGCG-8",  "TTCTTGACAGGAATCG-8" , "AGAGCAGCACTCCGAG-9"  ,"AGCGATTTCATTACCT-9" ,
                                         "ATGAGGGCACGTGAGA-9" , "CATCCACAGGGTTTCT-9" , "CTGAGCGGTAGTCACT-9" , "CTTTCGGAGCGTGAGT-9",  "GACGCTGAGTTCACTG-9" ,
                                         "GAGTCATAGGTAGCCA-9" , "GGGACTCCACGATTCA-9",  "GTCACTCTCAAAGGTA-9"  ,"GTCATTTCAGATACTC-9" , "TACTGCCAGATGAAGG-9" ,
                                         "TTGTGGAGTGAGATAT-9" , "AACGTCAGTGAGTTGG-10", "AGCGCCAGTTGGGATG-10", "AGGCATTTCGCTTAAG-10" ,"AGGGTGAGTGGCGTAA-10",
                                         "ATAGACCTCCACGGAC-10", "ATATCCTTCGACACCG-10" ,"ATGATCGAGCTATCTG-10" ,"CAGCAGCAGTTGGAAT-10", "CGAAGGAAGGTACTGG-10",
                                         "CTATAGGTCGACACTA-10", "CTGATCCCACGTTCGG-10" ,"GAATCACTCCTGTAGA-10" ,"GCCAGCAAGGAAGTGA-10" ,"GGGACAAAGGAAAGGT-10",
                                         "TACTTACCAAGGCGTA-10" ,"TCGTCCAAGTCCGTCG-10", "TGAACGTTCGACCAAT-10" ,"TGACAGTCATAGATGA-10" ,"TTACGTTAGGTCACAG-10",
                                         "TTCGCTGGTCGTGGAA-10" ,"TTTAGTCGTATGACAA-10", "TTTGGTTTCCAAATGC-10", "ACTATGGAGGCCTGCT-11", "AGATCCAAGTCCCGGT-11",
                                         "AGCCAATCAAGAAACT-11" ,"CAGATCAGTATGAGGC-11" ,"GAGAAATGTCGTCATA-11", "GGCTGTGAGTCGCCAC-11" ,"TACACCCGTGGAAATT-11",
                                          "TCCATGCAGGCAGTCA-11", "TTACAGGTCGCAGTGC-11", "TTCAATCTCAGACCTA-11", "ACGTCCTGTTGCTCGG-12", "AGGCCACTCGGACGTC-12",
                                         "CACTGGGTCGGACTGC-12" ,"CAGCGTGGTCCCTAAA-12", "CGATGGCAGACCAAAT-12" ,"CGGCAGTAGCTTCATG-12", "CGGCAGTGTTCTGACA-12",
                                         "CTCAACCCATCGAGCC-12" ,"GAGTTACGTCTCTCCA-12", "GATCACAGTGCCCGTA-12", "GCCAGCAGTGAATGTA-12" ,"GGCTTTCAGCGTCAAG-12",
                                         "GTCTACCAGAAGAACG-12" ,"GTTATGGCATGGCTAT-12", "TATTTCGAGCTCGTGC-12", "TGCCGAGGTGGAAATT-12", "TTCATTGTCGTTCTCG-12",
                                         "AGTCACACAAGCTGCC-13" ,"CAAGCTACACTCAAGT-13", "CGAGAAGAGACAACTA-13", "CGTCCATAGCGACTGA-13", "CTCATTAAGTACTGTC-13",
                                         "GACTCAAGTTCTTAGG-13" ,"GCAGTTAGTGCCCGTA-13", "GCAGTTATCTCTGCCA-13", "GGGAGTAGTGAGCTCC-13", "GGGCCATTCGAACTCA-13",
                                         "GTTTACTAGCCAGAGT-13" ,"TATGTTCCAAAGGATT-13", "TCATGTTTCGGATTAC-13", "TCGGATAAGCCATTGT-13" ,"TGTGTGACAGTAACCT-13",
                                         "TTCTCTCAGATAACAC-13" ,"AAGACAACAGCGCGTT-14", "AGAGAGCGTGGACTGA-14", "AGAGCAGTCACCGCTT-14", "AGTACTGTCATAAGGA-14",
                                         "ATCACAGCACCTAAAC-14" ,"ATCGTGAAGGGCTAAC-14", "ATTGTTCCAATTTCCT-14", "CATGAGTGTGTTATCG-14", "CATTCATTCAGGACAG-14",
                                         "CGCCATTTCCATTTAC-14" ,"CGCGTGACATTGACAC-14", "GACGCTGTCATGAGGG-14", "GAGTGAGTCAGTAGGG-14" ,"GTAATGCTCGAGTCCG-14",
                                         "GTTAGTGTCAAGCCTA-14" ,"TACGTCCGTCTCGACG-14", "TCCGTGTCACTCTCGT-14", "TCCTTTCGTCGCCTAG-14" ,"TGGATGTGTGGAATGC-14",
                                         "TGGGAGACAACAGAGC-14" ,"TGGGCGTCACAAGTGG-14", "TTCATGTCATAGATGA-14", "TTGACCCCAGGAATAT-14" ,"TTGCCTGCAAGTCGTT-14",
                                         "TTTACGTTCTATTCGT-14", "AACCATGCACGCTGAC-15", "AAGCCATTCGAATGCT-15", "AGAACAATCGTTAGTG-15", "CAGATCATCACGGGCT-15",
                                         "CATCGGGTCTAGACCA-15", "CATGCTCAGTTACGGG-15", "CCATCACGTCAGATTC-15" ,"CCTCACAGTAGGTCAG-15", "CGCATGGAGGTAGTCG-15",
                                         "CGGGCATCAGCGGATA-15" ,"CTAACCCCAACCGCCA-15", "GAAATGACAGCAGTAG-15", "GTAGTACGTTTCTATC-15", "TCACTATGTCAGTCTA-15",
                                         "TCGCTCAGTACACGTT-15" ,"TCGTGGGCACTGGAAG-15", "TGATGCAGTCTTGGTA-15", "TGCTTGCAGGACAAGA-15", "TGGAGAGTCATCGCTC-15",
                                         "TGTTGGAGTCAGCTTA-15", "ACATTTCAGTGAGTTA-16", "AGAGAATAGCTCCACG-16", "AGTCTCCTCTCGGTCT-16", "CAAAGAACACAGAGCA-16",
                                         "CAAGAGGGTACTCCGG-16", "CATCGTCCATCCCACT-16", "CCATAAGCAGAATCGG-16", "CTCACTGCAGAGGCTA-16", "CTGGACGGTACACGCC-16",
                                         "GAACTGTCAACGTATC-16", "GAACTGTGTACCTAAC-16", "GAGGGATCAAATTGGA-16", "GGTGTCGCATTCCTAT-16", "TACTTGTGTACCTAGT-16",
                                         "TATCGCCCACCTCTGT-16", "TCCATCGGTTGTACGT-16", "TGTGCGGGTGAGCAGT-16", "TTAGGGTAGCAGTCTT-16" ,"TTCATTGAGATGCCGA-16",
                                         "TTGAGTGCAGTCTTCC-16", "TTTACGTCAACTACGT-16", "TTTCCTCCAATACGAA-16", "ATGAGGGAGGCGATAC-17", "GAATAGAGTACTCCGG-17",
                                         "TGCCGAGGTTGACTGT-17", "ACTGCAACACTGAGGA-18" ,"AGAGAGCAGGATTTCC-18", "AGGGAGTCATGACACT-18" ,"ATACCGATCCGTGTGG-18",
                                         "CAGATTGCACACCGAC-18", "GACTCAAAGGTTAGTA-18", "TCATGAGAGCATGATA-18", "TCTATCATCGACATAC-18" ,"TGTAAGCCAGAGTTGG-18",
                                         "TTTGACTTCCTAGCGG-18"), 
                             SampleID = c(rep("Sample-43", 13), rep("Sample-50", 18), rep("Sample-37", 18), rep("Sample-52", 9), rep("Sample-96", 13), 
                                        rep("Sample-33", 22),rep("Sample-46", 10), rep("Sample-90", 17), rep("Sample-58", 16), rep("Sample-82", 25), 
                                        rep("Sample-45", 20), rep("Sample-100", 22),rep("Sample-66", 3), rep("Sample-47", 10)),
                             
                             Diagnosis = c(rep("AD", 13), rep("AD", 18), rep("AD", 18), rep("Control", 9), rep("Control", 13), rep("AD", 22), 
                                        rep("AD", 10), rep("Control", 17), rep("Control", 16), rep("Control", 25), rep("AD", 20), rep("Control", 22),
                                        rep("Control", 3), rep("AD", 10)),
                             
                             Batch = c(rep("1", 13), rep("2", 18), rep("3", 18), rep("3", 9), rep("1", 13), rep("3", 22), 
                                            rep("2", 10), rep("3", 17), rep("3", 16), rep("2", 25), rep("1", 20), rep("1", 22),
                                            rep("2", 3), rep("3", 10)),
                             
                             Age = c(rep("90", 13), rep("89", 18), rep("87", 18), rep("83", 9), rep("79", 13),
                                      rep("80", 22), rep("90", 10), rep("79", 17), rep("90", 16), rep("79", 25),
                                      rep("89", 20), rep("79", 22),rep("90", 3), rep("90",10)),
                             
                             Sex = c(rep("F", 13), rep("F", 18), rep("F", 18), rep("M", 9), rep("M", 13), rep("M", 22), 
                                            rep("M", 10), rep("F", 17), rep("M", 16), rep("M", 25), rep("F", 20), 
                                      rep("M", 22), rep("F", 3), rep("M", 10)),
                             
                             PMI = c(rep("4.17", 13), rep("4.08", 18), rep("4.25", 18), rep("1.8", 9), rep("7", 13),
                                      rep("3.75", 22), 
                                            rep("4.23", 10), rep("4.25", 17), rep("4", 16), rep("5", 25), rep("3.88", 20),
                                      rep("NA", 22),rep("2.92", 3), rep("5", 10)),
                             
                             Tangle.Stage = c(rep("Stage 6", 13), rep("Stage 6", 18), rep("Stage 6", 18),
                                               rep("Stage 2", 9), rep("Stage 1", 13), rep("Stage 5", 22), 
                                            rep("Stage 5", 10), rep("Stage 1", 17), rep("Stage 2", 16), rep("Stage 1", 25), 
                                            rep("Stage 6", 20), rep("Stage 2", 22),rep("Stage 2", 3), rep("Stage 5", 10)),
                             
                             Plaque.Stage = c(rep("Stage B", 13), rep("Stage C", 18), rep("Stage C", 18), 
                                               rep("None", 9), rep("Stage A", 13), rep("Stage C", 22), 
                                            rep("Stage C", 10), rep("Stage A", 17), rep("Stage B", 16), rep("Stage A", 25), 
                                            rep("Stage B", 20), rep("Stage B", 22),rep("None", 3), rep("Stage B", 10)),
                             
                             RIN = c(rep("8.9", 13), rep("8.3", 18), rep("7", 18), rep("7.2", 9), rep("7.1", 13), 
                                      rep("7.3", 22),rep("7.9", 10), rep("7.3", 17), rep("8.1", 16), rep("9", 25), 
                                      rep("7.8", 20), rep("8.7", 22),rep("10", 3), rep("8", 10)))


# Update the SampleID and other column info in the metadata for the missing barcodes
for (i in 1:nrow(new_sample_ids)) {
  barcode <- new_sample_ids$Barcode[i]
  sample <- new_sample_ids$SampleID[i]
  diagnose <- new_sample_ids$Diagnosis[i]
  batch <- new_sample_ids$Batch[i]
  age <- new_sample_ids$Age[i]
  sex <- new_sample_ids$Sex[i]
  pmi <- new_sample_ids$PMI[i]
  tangle <- new_sample_ids$Tangle.Stage[i]
  plaque <- new_sample_ids$Plaque.Stage[i]
  rin <- new_sample_ids$RIN[i]

  # Check if the barcode exists in the Seurat object
  if (barcode %in% rownames(seurat_hdf5@meta.data)) {
    seurat_hdf5@meta.data[barcode, "SampleID"] <- sample
    seurat_hdf5@meta.data[barcode, "Diagnosis"] <- diagnose
    seurat_hdf5@meta.data[barcode, "Batch"] <- batch
    seurat_hdf5@meta.data[barcode, "Age"] <- age
    seurat_hdf5@meta.data[barcode, "Sex"] <- sex
    seurat_hdf5@meta.data[barcode, "PMI"] <- pmi
    seurat_hdf5@meta.data[barcode, "Tangle.Stage"] <- tangle
    seurat_hdf5@meta.data[barcode, "Plaque.Stage"] <- plaque
    seurat_hdf5@meta.data[barcode, "RIN"] <- rin
  }
}


# Verify that the SampleID has been updated
updated_metadata <- seurat_hdf5@meta.data[missing_barcodes, ]
head(updated_metadata)

dim(updated_metadata)

rm(batch,sex, age,rin,pmi,Age,Batch,diagnose,Diagnosis,Sex,RIN,PMI,Tangle.Stage,tangle,pla4
   ,Plaque.Stage)
rm(plaque,barcode,i,missing_barcodes,sample)