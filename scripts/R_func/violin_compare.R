library(ggplot2)
library(dplyr)
library(ggsignif) 

plot_comparison <- function(overlap_data1, overlap_data2, 
                          plot_title = "Comparison", data1_lab='data1',data2_lab='data2',
                          y_axis = "log2_value",transform_type='log2',
                          feature_col = 8,
                          args_bw = 0.3,
                          fill_colors = c("indianred", "lightblue"),
                          y_limits_fold = c(0.5, 1.5)) {
  
    # Prepare input data
    in_file_1 <- overlap_data1[, c(1, 2, 3, feature_col)]
    colnames(in_file_1) <- c('V1', 'V2', 'V3', 'V4')
    in_file_1$subgroup <- data1_lab
    
    in_file_2 <- overlap_data2[, c(1, 2, 3, feature_col)]
    colnames(in_file_2) <- c('V1', 'V2', 'V3', 'V4')
    in_file_2$subgroup <- data2_lab

    # Combine data
    all_data <- rbind(in_file_1, in_file_2)

    # Transform feature column based on transform_type
    if (transform_type == 'log2') {
        all_data$feature_val <- ifelse(all_data$V4 > 0, log2(all_data$V4), NA)
    } else if (transform_type == 'log10') {
        all_data$feature_val <- ifelse(all_data$V4 > 0, log10(all_data$V4), NA)
    } else if (transform_type == 'NA') {
        all_data$feature_val <- all_data$V4
    }
  
    # Statistical test
    wilcox_result <- wilcox.test(in_file_2$V4, in_file_1$V4, alternative = "two.sided")
    y_limits = c(0, 0)
    y_limit_val_hi <- quantile(all_data$feature_val, 0.999, na.rm = TRUE)
    y_limits[2] <- ifelse(y_limit_val_hi < 0, y_limit_val_hi * y_limits_fold[1], y_limit_val_hi * y_limits_fold[2])
    y_limit_val_low <- quantile(all_data$feature_val, 0.001, na.rm = TRUE)
    y_limits[1] <- ifelse(y_limit_val_low < 0, y_limit_val_low * y_limits_fold[2], y_limit_val_low * y_limits_fold[1])
    y_pos <- ifelse(y_limit_val_hi < 0, y_limit_val_hi * 0.8, y_limit_val_hi * 1.1)

    # Summary stats
    summ <- all_data %>%
        group_by(subgroup) %>%
        dplyr::summarize(
        n = n(), 
        mean = round(mean(feature_val, na.rm = TRUE), 2),
        max_val = round(max(feature_val, na.rm = TRUE), 2),
        sd = sd(feature_val, na.rm = TRUE),
        median = round(median(feature_val, na.rm = TRUE), 2)
        )

    # Create plot
    p <- ggplot(all_data, aes(x = subgroup, y = feature_val, fill = subgroup)) +
        geom_violin(trim = FALSE, bw = args_bw, na.rm = TRUE) +
        geom_boxplot(width = 0.1, outlier.shape = NA, na.rm = TRUE) +
        geom_text(aes(label = paste0('N=', n), y = median), 
                data = summ, size = 4, vjust = 2, hjust = 1.5) +
        geom_text(aes(label = paste0('mean= ', mean), y = median), 
                data = summ, size = 4, vjust = 0, hjust = 1.5) +
        scale_fill_manual(values = fill_colors) + theme_classic() +
        geom_signif(comparisons = list(levels(factor(all_data$subgroup))), 
                    test = "wilcox.test", map_signif_level = TRUE,
                    y_position = y_pos,size=1) +
        labs(title = plot_title,
            subtitle = paste0("Wilcoxon p-value = ", format.pval(wilcox_result$p.value)),
            y = y_axis) +
        theme(axis.text = element_text(face = 'bold'),
            axis.title = element_text(face = "bold"),
            plot.title = element_text(face = "bold"),
            legend.position = "top") +
        ylim(y_limits[1], y_limits[2])
#   pdf(paste0(plot_title, '.pdf'))
#   print(p)
#   dev.off()
  return(p)
}