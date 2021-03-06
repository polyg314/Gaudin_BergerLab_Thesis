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
  standard_devs_NA <- c() #not absolute
  qmin <- c()
  q25 <- c()
  q50 <- c()
  q75 <- c()
  qmax <- c()
  chrom_single_max <- c()
  chrom_double_max <- c()
  chrom_triple_max <- c()
  for(i in 3:ncol(coverage_filtered_normal_table)){
    index <- c(index,i)
    current_test <- coverage_filtered_normal_table[,1:2]
    current_test$Ref <- coverage_filtered_normal_table[,i]
    current_pool <- coverage_filtered_normal_table[,-c(i)] ##remove current normal sample
    current_test <- condense.tumor(current_test, min_alt_freq_t, max_alt_freq_t)
    sample_v_normal_pool_table <- het.analysis.table(current_test, current_pool, min_normal_hets, min_alt_freq_n, max_alt_freq_n)
    #print(sample_v_normal_pool_table)
    sample_v_normal_pool_table <- sample_v_normal_pool_table[which(sample_v_normal_pool_table$normal_mean > 25 & sample_v_normal_pool_table$normal_mean < 75),]
    #print(sample_v_normal_pool_table)
    standard_devs <- c(standard_devs,sum(abs(sample_v_normal_pool_table$tdev_from_n))/nrow(sample_v_normal_pool_table))
    standard_devs_NA <- c(standard_devs_NA,sum(sample_v_normal_pool_table$tdev_from_n)/nrow(sample_v_normal_pool_table)) 
    het_counts <- c(het_counts, nrow(sample_v_normal_pool_table))
    x <- quantile(abs(sample_v_normal_pool_table$tdev_from_n))
    qmin <- c(qmin, x[1])
    q25 <- c(q25, x[2])
    q50 <- c(q50, x[3])
    q75 <- c(q75, x[4])
    qmax <- c(qmax, x[5])
    chrom_zs <- c()
    for(i in 1:23){
      current_s_v_n <- sample_v_normal_pool_table[which(sample_v_normal_pool_table$Chromosome == i),]
      if(nrow(current_s_v_n) > 5){
        chrom_zs <- c(chrom_zs, sum(abs(current_s_v_n$tdev_from_n))/nrow(current_s_v_n))
      }else{
        chrom_zs <- c(chrom_zs, 1)
      }
    }
    n <- length(chrom_zs)
    chrom_zs <- sort(chrom_zs,partial=n-1)
    chrom_single_max <- c(chrom_single_max, chrom_zs[n])
    chrom_double_max <- c(chrom_double_max, mean(c(chrom_zs[n], chrom_zs[n-1])))
    chrom_triple_max <- c(chrom_triple_max, mean(c(chrom_zs[n], chrom_zs[n-1], chrom_zs[n-2])))
  }
  for(i in 3:ncol(coverage_filtered_tumor_table)){
    index <- c(index,i)
    current_test <- coverage_filtered_tumor_table[,1:2]
    current_test$Ref <- coverage_filtered_tumor_table[,i]
    current_test <- condense.tumor(current_test, min_alt_freq_t, max_alt_freq_t)
    sample_v_normal_pool_table <- het.analysis.table(current_test, coverage_filtered_normal_table, min_normal_hets, min_alt_freq_n, max_alt_freq_n)
    sample_v_normal_pool_table <- sample_v_normal_pool_table[which(sample_v_normal_pool_table$normal_mean > 25 & sample_v_normal_pool_table$normal_mean < 75),]
    standard_devs <- c(standard_devs,sum(abs(sample_v_normal_pool_table$tdev_from_n))/nrow(sample_v_normal_pool_table))
    standard_devs_NA <- c(standard_devs_NA,sum(sample_v_normal_pool_table$tdev_from_n)/nrow(sample_v_normal_pool_table)) 
    het_counts <- c(het_counts, nrow(sample_v_normal_pool_table))
    x <- quantile(abs(sample_v_normal_pool_table$tdev_from_n))
    qmin <- c(qmin, x[1])
    q25 <- c(q25, x[2])
    q50 <- c(q50, x[3])
    q75 <- c(q75, x[4])
    qmax <- c(qmax, x[5])
    chrom_zs <- c()
    for(i in 1:23){
      current_s_v_n <- sample_v_normal_pool_table[which(sample_v_normal_pool_table$Chromosome == i),]
      if(nrow(current_s_v_n) > 5){
        chrom_zs <- c(chrom_zs, sum(abs(current_s_v_n$tdev_from_n))/nrow(current_s_v_n))
      }else{
        chrom_zs <- c(chrom_zs, 1)
      }
    }
    n <- length(chrom_zs)
    chrom_zs <- sort(chrom_zs,partial=n-1)
    chrom_single_max <- c(chrom_single_max, chrom_zs[n])
    chrom_double_max <- c(chrom_double_max, mean(c(chrom_zs[n], chrom_zs[n-1])))
    chrom_triple_max <- c(chrom_triple_max, mean(c(chrom_zs[n], chrom_zs[n-1], chrom_zs[n-2])))
  }
  normal_tumor_seq <- c(rep("normal", number_normal_samples), rep("tumor", number_tumor_samples))
  Z_score_df <- data.frame(pool = normal_tumor_seq, index = index, sample_names = normal_and_tumor_names, VAF_Z_Score = standard_devs, number_of_hets_avgd = het_counts, VAF_Z_NA = standard_devs_NA, chrom_single_max = chrom_single_max, chrom_double_max = chrom_double_max, chrom_triple_max = chrom_triple_max, qmin = qmin, q25 = q25, q50 = q50, q75 = q75, qmax = qmax)
  return(Z_score_df)
}

## gene specific table function 

gene_spec_table <- function(counts.table, chromosome, start_position, end_position){
  chrom.spec.table <- counts.table %>% filter(Chromosome == chromosome)
  position.spec.table <- chrom.spec.table %>% filter(Position > start_position)
  position.spec.table <- position.spec.table %>% filter(Position < end_position)
  return(position.spec.table)
}


## get het counts (only) for each sample within a gene 
gene_het_count_table <- function(gene_table,alt_cov_table, min_alt_freq, max_alt_freq, minimum_coverage){
  coverage_filtered_table <- filter_by_min_coverage(alt_cov_table,minimum_coverage)
  for(j in 1:nrow(gene_table)){
    current_gene_table <- gene_spec_table(coverage_filtered_table,gene_table$Chromosome[j],gene_table$Start[j],gene_table$Stop[j])
    het_counts <- c()
    for(i in 3:ncol(current_gene_table)){
      het_count <- length(which(current_gene_table[,i] > min_alt_freq & current_gene_table[,i] < max_alt_freq))
      het_counts <- c(het_counts,het_count)
    }
    if(j == 1){
      gene_counts_mat <- matrix(data = het_counts, ncol = 1)
    }else{
      gene_counts_mat <- cbind(gene_counts_mat,het_counts)
    }
  }
  gene_counts_df <- data.frame(gene_counts_mat)
  colnames(gene_counts_df) <- gene_table$GENE
  gene_counts_table <- gene_counts_df
  return(gene_counts_table)
}

## Gene Specific VAF Z-scores
gene_specific_tables <- function(gene_table, tumor_tables, normal_tables, tumor_names, min_alt_freq_n, max_alt_freq_n, min_alt_freq_t, max_alt_freq_t, min_normal_hets, minimum_coverage){
  coverage_filtered_normal_table <- filter_by_min_coverage(normal_tables,minimum_coverage)
  coverage_filtered_tumor_table <- filter_by_min_coverage(tumor_tables,minimum_coverage)
  for(j in 1:nrow(gene_table)){
    current_gene_tumor_table <- gene_spec_table(coverage_filtered_tumor_table,gene_table$Chromosome[j],gene_table$Start[j],gene_table$Stop[j])
    current_gene_normal_table <- gene_spec_table(coverage_filtered_normal_table,gene_table$Chromosome[j],gene_table$Start[j],gene_table$Stop[j])
    standard_devs <- c()
    het_counts <- c()
    avg_percent_dif <- c()
    for(i in 3:ncol(current_gene_tumor_table)){
      current_test <- current_gene_tumor_table[,1:2]
      current_test$Ref <- current_gene_tumor_table[,i]
      current_test <- condense.tumor(current_test, min_alt_freq_t, max_alt_freq_t)
      if(nrow(current_test) < 1){
        standard_devs <- c(standard_devs,0)
        het_counts <- c(het_counts,0)
        avg_percent_dif <- 0
      }
      else{
        sample_v_normal_pool_table <- het.analysis.table(current_test, current_gene_normal_table, min_normal_hets, min_alt_freq_n, max_alt_freq_n)
        sample_v_normal_pool_table <- sample_v_normal_pool_table[which(sample_v_normal_pool_table$normal_mean > 25 & sample_v_normal_pool_table$normal_mean < 75),]
        if(nrow(sample_v_normal_pool_table) > 0){
          avg_abs_z <- sum(abs(sample_v_normal_pool_table$tdev_from_n))/nrow(sample_v_normal_pool_table)
          number_hets <- nrow(sample_v_normal_pool_table)
          standard_devs <- c(standard_devs,avg_abs_z)
          het_counts <- c(het_counts,number_hets)            
        }
        else{
          standard_devs <- c(standard_devs,0)
          het_counts <- c(het_counts,0)
          avg_percent_dif <- 0          
        }
      }
    }
    if(j == 1){
      gene_z_mat <- matrix(data = standard_devs, ncol = 1)
      gene_counts_mat <- matrix(data = het_counts, ncol = 1)
    }else{
      gene_z_mat <- cbind(gene_z_mat,standard_devs)
      gene_counts_mat <- cbind(gene_counts_mat,het_counts)
    }
  }
  gene_z_df <- data.frame(gene_z_mat)
  gene_counts_df <- data.frame(gene_counts_mat)
  colnames(gene_z_df) <- gene_table$GENE
  colnames(gene_counts_df) <- gene_table$GENE
  gene_z_df$sample_name <- tumor_names
  gene_counts_df$sample_name <- tumor_names
  gene_counts_table <- list(gene_z_df, gene_counts_df)
  return(gene_counts_table)
}

##bin fragments into bins based on bin size 
bin_fragments = function(frag_lengths,bin_size,lower,upper,sample_name){
  sample_name <- as.character(sample_name)
  # make absolute fragment length for reversed segments
  frag_lengths = abs(frag_lengths)
  # filter for upper and lower range fragment sizes
  frag_lengths_vector <- frag_lengths[frag_lengths <= as.numeric(upper+1) & frag_lengths >= as.numeric(lower)]
  # tallying by each size
  frag_lengths_vector_tally = table(frag_lengths_vector)
  # making bin size table explicit
  if(upper %% bin_size == 0){
    bins_start = seq(lower,upper,bin_size)
    bins_end = bins_start+bin_size-1
  }else{
    bins_start = seq(lower,upper,bin_size)
    bins_end = ifelse(bins_start+bin_size-1 > upper, upper,bins_start+bin_size-1)
  }
  bins_dt = data.table(start = bins_start,end = bins_end)
  # counting the sum of fragments within each bin
  bins_dt[,eval(sample_name) := sum(frag_lengths_vector_tally[start:end]),.(start,end)]
}


normalize_frag_table <- function(frag_table){
  end <- ncol(frag_table)-1 #don't include sample name column 
  for(i in 1:nrow(frag_table)){
    frag_table[i,1:end] <- as.numeric(frag_table[i,1:end])/sum(as.numeric(frag_table[i,1:end]))
  }
  return(frag_table)
}


calculate_frag_z_scores <- function(tumor_frag_table, normal_frag_table){
  end <- ncol(normal_frag_table) - 1
  normal_z_table <- normal_frag_table[,-c(ncol(normal_frag_table))]
  for(i in 1:end){
    for(j in 1:nrow(normal_frag_table)){
      temp_normal_frag_table <- normal_frag_table[-(j),]
      normal_mean <- mean(temp_normal_frag_table[,i])
      normal_sd <- sd(temp_normal_frag_table[,i])
      normal_z_table[j,i] <- abs((normal_frag_table[j,i] - normal_mean)/normal_sd)
    }
  }
  normal_z_avg_table <- data.frame(sample_names = normal_frag_table$sample_names, abs_frag_z_score = rowMeans(normal_z_table))
  end <- ncol(tumor_frag_table)-1
  tumor_z_table <- tumor_frag_table[,-c(ncol(tumor_frag_table))]
  for(i in 1:end){
    normal_mean <- mean(normal_frag_table[,i])
    normal_sd <- sd(normal_frag_table[,i])
    for(j in 1:nrow(tumor_frag_table)){
      tumor_z_table[j,i] <- abs((tumor_frag_table[j,i] - normal_mean)/normal_sd)
    }
  }
  tumor_z_avg_table <- data.frame(sample_names = tumor_frag_table$sample_names, abs_frag_z_score = rowMeans(tumor_z_table))
  return(rbind(normal_z_avg_table,tumor_z_avg_table))
}


ROC_curve_function <- function(sd_df, Z_score){
  normal_df <- sd_df[which(sd_df$pool =="normal"),]
  tumor_df <- sd_df[which(sd_df$pool =="tumor"),]
  max_normal_vaf_z <- max(normal_df$VAF_Z_Score)
  max_normal_frag_z <- max(normal_df$abs_frag_z_score)
  triple_max_normal_vaf_z <- max(normal_df$chrom_triple_max)
  if(Z_score == "VAF"){
    sd_df$sd_from_normals <- sd_df$VAF_Z_Score
    print("ctDNA flagged positive at 100% specificity")
    sum(tumor_df$VAF_Z_Score > max_normal_vaf_z)
    print("Number of tumor samples")
    print(length(tumor_df$VAF_Z_Score))
    print("Number of tumor samples above max normal vaf z")
    print(sum(tumor_df$VAF_Z_Score > max_normal_vaf_z))
    print("Sensitivity")
    print(sum(tumor_df$VAF_Z_Score > max_normal_vaf_z)/length(tumor_df$VAF_Z_Score))
  }
  else if(Z_score == "top_3"){
    triple_max_normal_vaf_z <- max(normal_df$chrom_triple_max)
    sd_df$sd_from_normals <- sd_df$chrom_triple_max
    print("ctDNA flagged positive at 100% specificity")
    sum(tumor_df$chrom_triple_max > triple_max_normal_vaf_z)
    print("Number of tumor samples")
    print(length(tumor_df$chrom_triple_max))
    print("Number of tumor samples above max normal vaf z")
    print(sum(tumor_df$chrom_triple_max > triple_max_normal_vaf_z))
    print("Sensitivity")
    print(sum(tumor_df$chrom_triple_max > triple_max_normal_vaf_z)/length(tumor_df$chrom_triple_max))
  }
  else if(Z_score == "fragment_length"){
    sd_df$sd_from_normals <- sd_df$abs_frag_z_score
    print("ctDNA flagged positive at 100% specificity")
    sum(tumor_df$VAF_Z_Score > max_normal_vaf_z )
    print("Number of tumor samples")
    print(length(tumor_df$abs_frag_z_score))
    print("Number of tumor samples above max frag z")
    print(sum(tumor_df$abs_frag_z_score > max_normal_frag_z))
    print("Sensitivity")
    print(sum(tumor_df$abs_frag_z_score > max_normal_frag_z)/length(tumor_df$abs_frag_z_score))
  }
  else if(Z_score == "frag_with_VAF_cutoff"){
    sd_df$sd_from_normals <- sd_df$abs_frag_z_score
    print("ctDNA flagged positive at 100% specificity")
    sum(tumor_df$abs_frag_z_score > max_normal_frag_z | tumor_df$VAF_Z_Score > max_normal_vaf_z | tumor_df$chrom_triple_max > triple_max_normal_vaf_z)
    print("Number of tumor samples")
    print(length(tumor_df$abs_frag_z_score))
    print("Number of tumor samples above max frag z or max VAF frag z")
    print(sum(tumor_df$abs_frag_z_score > max_normal_frag_z | tumor_df$VAF_Z_Score > max_normal_vaf_z | tumor_df$chrom_triple_max > triple_max_normal_vaf_z))
    print("Sensitivity")
    print(sum(tumor_df$abs_frag_z_score > max_normal_frag_z | tumor_df$VAF_Z_Score > max_normal_vaf_z | tumor_df$chrom_triple_max > triple_max_normal_vaf_z)/length(tumor_df$VAF_Z_Score))
  }
  labels_vector <- as.character(sd_df$pool)
  labels_vector[which(labels_vector == "tumor")] <- "1"
  labels_vector[which(labels_vector == "normal")] <- "0"
  labels_vector <- as.numeric(labels_vector)
  number_normal_samples <- 47
  sd_cutoff <- seq(0,15,.001)
  output_performance <- matrix(rep(0,nrow(sd_df)*length(sd_cutoff)),nrow(sd_df),length(sd_cutoff))
  label_matrix <- matrix(rep(0,nrow(sd_df)*length(sd_cutoff)),nrow(sd_df),length(sd_cutoff))
  sd_matrix <- matrix(rep(0,nrow(sd_df)*length(sd_cutoff)),nrow(sd_df),length(sd_cutoff))
  for(i in 1:length(sd_cutoff)){
    current_cuttoff <- sd_cutoff[i]
    predictions_vector <- c()
    for(j in 1:nrow(sd_df)){
      if(Z_score == "frag_with_VAF_cutoff"){
        if(sd_df$VAF_Z_Score[j] > max_normal_vaf_z | sd_df$chrom_triple_max[j] > triple_max_normal_vaf_z){
          predictions_vector <- c(predictions_vector, 1)
        }
       else if(sd_df$sd_from_normals[j] < current_cuttoff){
          predictions_vector <- c(predictions_vector, 0)
       }
       else{
         predictions_vector <- c(predictions_vector, 1)
        }
      }
      else{
        if(sd_df$sd_from_normals[j] < current_cuttoff){
          predictions_vector <- c(predictions_vector, 0)
        }
        else{
          predictions_vector <- c(predictions_vector, 1)
        }
      }
    }
    output_performance[,i] <- predictions_vector
    label_matrix[,i] <- labels_vector 
    sd_vector <- c(labels_vector, c(rep(sd_cutoff[i], number_normal_samples), rep(sd_cutoff[i], nrow(sd_df) - number_normal_samples)))
    sd_matrix[,i] <- labels_vector 
  }
  tpv <- c()
  fpv <- c()
  accuracy_v <- c()
  for(i in 1:ncol(output_performance)){
    true_pos <- 0
    false_pos <- 0
    true_neg <- 0
    false_neg <- 0
    for(j in 1:nrow(output_performance)){
      if(output_performance[j,i] == 1 & label_matrix[j,i] == 1){
        true_pos <- true_pos + 1
      } 
      else if(output_performance[j,i] == 0 & label_matrix[j,i] == 1){
        false_neg <- false_neg + 1
      }
      else if(output_performance[j,i] == 0 & label_matrix[j,i] == 0){
        true_neg <- true_neg + 1
      }
      else{
        false_pos <- false_pos + 1
      }
    }
    tpr <- true_pos/(true_pos + false_neg)
    fpr <- false_pos/(false_pos + true_neg)
    accuracy <- (true_pos + true_neg)/nrow(output_performance)
    tpv <- c(tpv, tpr)
    fpv <- c(fpv, fpr)
    accuracy_v <- c(accuracy_v, accuracy)
  }
  roc_df <- data.frame(True_positive_rate = tpv, False_positive_rate = fpv, accuracy = accuracy_v, abs_z_score_cutoffs = sd_cutoff)
  print_ROC_AUC(roc_df)
  return(roc_df)
}

print_ROC_AUC <- function(roc_df){
  fp_rates <- roc_df$False_positive_rate
  fp_rates_2 <- c(1,fp_rates)
  horizontal_distance <- abs(fp_rates - fp_rates_2[1:length(fp_rates)])
  AUC <- sum(roc_df$True_positive_rate * horizontal_distance)
  print("AUC")
  print(AUC)
}
