library(eulerr)
library(ggsci)
library(UpSetR)
library(ComplexHeatmap)

plot_overlaps <- function(multi_overlap, set_names = NULL) {
    multi_overlap_mat<-multi_overlap[,6:length(colnames(multi_overlap))]
    colnames(multi_overlap_mat) <- set_names
    # Get npg colors
    npg_colors <- pal_npg()(10)
    
    # Create Euler plot
    combinations <- apply(multi_overlap_mat, 1, function(row) {
        paste(colnames(multi_overlap_mat)[row == 1], collapse = "&")
    })
    counts <- table(combinations)
    set_counts <- as.numeric(counts)
    names(set_counts) <- names(counts)
    fit2 <- euler(set_counts)
    euler_plot <- plot(fit2, quantities = TRUE, fills = npg_colors[c(1:length(set_names))])
    # Create UpSet plot
    m=make_comb_mat(multi_overlap_mat)
    upset_plot <- UpSet(m, comb_order = order(comb_size(m),decreasing = TRUE),
        # data=as.data.frame(multi_overlap_mat),
        # sets = set_names,
        # order.by = "freq",
        # mainbar.y.label = "Intersection Size",
        # sets.x.label = "Set Size"
    )
    
    # Return both plots in a list
    return(list(
        euler = euler_plot,
        upset = upset_plot
    ))
}
?order

# Example usage:
# multi_overlap<-bedtoolsr::bt.multiinter(list(paused_tss_bed,noPaused_tss_bed,MICC_1d),cluster=TRUE)
# head(multi_overlap)
# multi_overlap_mat<-multi_overlap[multi_overlap$V5!='1,2,3',6:length(colnames(multi_overlap))]
# cols<-c('paused_tss_bed','noPaused_tss_bed','MICC_1d')
# plots <- plot_overlaps('multi_overlap.xls',set_names=cols)
# plots$euler  # View Euler plot
# plots$upset  # View UpSet plot