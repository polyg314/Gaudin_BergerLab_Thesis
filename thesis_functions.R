##merge snp-pileup output files into single table
merge_count_files <- function(sample_paths){
  sample_counts_merged <- read.table(sample_paths[1], stringsAsFactors = F, header = TRUE, sep=",")
  sample_counts_merged$Position <- as.numeric(sample_counts_merged$Position)
  sample_counts_merged <- subset(sample_counts_merged, select = -c(3, 4))
  for(i in 2:length(sample_paths)){
    print(i)
    new_counts_merged <- read.table(sample_paths[i], stringsAsFactors = F, header = TRUE, sep=",")
    colnames(new_counts_merged)[5] <- paste("File",i,"R",sep="")
    colnames(new_counts_merged)[6] <- paste("File",i,"A",sep="")
    colnames(new_counts_merged)[7] <- paste("File",i,"E",sep="")
    colnames(new_counts_merged)[8] <- paste("File",i,"D",sep="")
    new_counts_merged <- subset(new_counts_merged, select = -c(3, 4))
    #new_counts_merged$Position <- as.numeric(new_counts_merged$Position)
    sample_counts_merged <- merge(sample_counts_merged, new_counts_merged, by=c("Chromosome", "Position"), all = TRUE)
  }
  sample_counts_merged[is.na(sample_counts_merged)] <- 0 
  return(sample_counts_merged)
}

#format table for further analysis
format_het_table <- function(input_countsmerged,min_alt_freq,max_alt_freq){
  counts_table <- input_countsmerged
  counts_table <- counts_table[grep("\\Ref|Alt|E|D", names(counts_table), invert=TRUE)]
  ref_table <- counts_table[,3:ncol(counts_table)]
  ref_table <- ref_table[grep("\\A", names(ref_table), invert=TRUE)]
  alt_table <- counts_table[,3:ncol(counts_table)]
  alt_table <- alt_table[grep("\\R", names(alt_table), invert=TRUE)]
  het_table <- counts_table[,1:ncol(counts_table)]
  het_table <- het_table[grep("\\A", names(het_table), invert=TRUE)]
  het_table[,3:ncol(het_table)] <- 100*(alt_table/(ref_table + alt_table))
  het_table[is.na(het_table)] <- 0 ## NA results from 0/0
  coverage_table <- counts_table[,1:ncol(counts_table)]
  coverage_table <- coverage_table[grep("\\A", names(coverage_table), invert=TRUE)]
  coverage_table[,3:ncol(coverage_table)] <- ref_table + alt_table
  data_frame_object <- list(het_table, coverage_table)
  return(data_frame_object)
}

##Filter het table - only include het sites above min_coverage
filter_by_min_coverage <- function(counts_tables,min_coverage){
  het_table <- counts_tables[[1]]
  coverage_table <- counts_tables[[2]]
  het_table[,3:ncol(het_table)][coverage_table[,3:ncol(coverage_table)] < min_coverage] <- 0 ##change VAF frequency not meeting min coverage to 0
  return(het_table)
}

##Condense Tumor File to Just Het Sites
condense.tumor <- function(het.freqs, min.cutoff, max.cutoff){
  het.freqs <- het.freqs %>% filter(Ref < max.cutoff)
  het.freqs <- het.freqs %>% filter(Ref > min.cutoff)
  colnames(het.freqs)[3] <- "alt_allele_freq"
  return(het.freqs)
}

##SNP table with normal means and SDs

het.analysis.table <- function(tumor.table, normal.het.table, min_number_normal_hets, min_alt_freq, max_alt_freq){
  full.table <- tumor.table
  means <- c()
  sds <- c()
  ns <- c()
  devs <- c()
  for(i in 1:nrow(full.table)){
    tmp <- normal.het.table[which(normal.het.table[,1] == full.table[i,1]),]
    row <- tmp[which(tmp[,2] == full.table[i,2] ),]
    #if(nrow(row) > 0){
    #print(dim(row[,3:ncol(row)]))
    hets <- row[,3:ncol(row)]
    het_list <- c()
    for(j in 1:ncol(hets)){
      het_list <- c(het_list, hets[1,j])
    }
    het_list <- het_list[het_list > min_alt_freq]
    het_list <- het_list[het_list < max_alt_freq]
    mn <- mean(het_list)
    means <- c(means, mn)
    sd <- sd(het_list)
    sds <- c(sds, sd)
    ns <- c(ns, length(het_list))
    dev <- (full.table[i,3] - mn) / sd
    devs <- c(devs, dev)
    #}
  }
  colnames(full.table)[3] <- c("AA_freq")
  full.table$normal_mean <- means
  full.table$normal_sd <- sds
  full.table$normal_ns <- ns
  full.table$tdev_from_n <- devs
  full.table <- full.table %>% filter(!is.na(normal_sd))
  full.table <- full.table %>% filter(normal_ns > min_number_normal_hets)
  return(full.table)
}

##Compute mean avg absolute Z-scores for each normal and tumor sample
VAF_Z_Score_analysis <- function(normal_table, tumor_table, minimum_coverage, min_alt_freq_n, max_alt_freq_n, min_alt_freq_t, max_alt_freq_t,min_normal_hets, normal_and_tumor_names){
  coverage_filtered_normal_table <- filter_by_min_coverage(normal_table,minimum_coverage)
  coverage_filtered_tumor_table <- filter_by_min_coverage(tumor_table,minimum_coverage)
  number_tumor_samples <- ncol(coverage_filtered_tumor_table) - 2
  number_normal_samples <- ncol(coverage_filtered_normal_table) - 2
  standard_devs <- c()
  index <- c()
  het_counts <- c()
  for(i in 3:ncol(coverage_filtered_normal_table)){
    index <- c(index,i)
    current_test <- coverage_filtered_normal_table[,1:2]
    current_test$Ref <- coverage_filtered_normal_table[,i]
    current_pool <- coverage_filtered_normal_table[,-c(i)] ##remove current normal sample
    current_test <- condense.tumor(current_test, min_alt_freq_t, max_alt_freq_t)
    sample_v_normal_pool_table <- het.analysis.table(current_test, current_pool, min_normal_hets, min_alt_freq_n, max_alt_freq_n)
    standard_devs <- c(standard_devs,sum(abs(sample_v_normal_pool_table$tdev_from_n))/nrow(sample_v_normal_pool_table))
    het_counts <- c(het_counts, nrow(sample_v_normal_pool_table))
  }
  for(i in 3:ncol(coverage_filtered_tumor_table)){
    index <- c(index,i)
    current_test <- coverage_filtered_tumor_table[,1:2]
    current_test$Ref <- coverage_filtered_tumor_table[,i]
    current_test <- condense.tumor(current_test, min_alt_freq_t, max_alt_freq_t)
    sample_v_normal_pool_table <- het.analysis.table(current_test, coverage_filtered_normal_table, min_normal_hets, min_alt_freq_n, max_alt_freq_n)
    standard_devs <- c(standard_devs,sum(abs(sample_v_normal_pool_table$tdev_from_n))/nrow(sample_v_normal_pool_table))
    het_counts <- c(het_counts, nrow(sample_v_normal_pool_table))
  }
  normal_tumor_seq <- c(rep("normal", number_normal_samples), rep("tumor", number_tumor_samples))
  Z_score_df <- data.frame(pool = normal_tumor_seq, index = index, sample_names = normal_and_tumor_names, VAF_Z_Score = standard_devs, number_of_hets_avgd = het_counts)
  return(Z_score_df)
}




