library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)


calculate_motif_density <- function(bed_df,motif, genome="hg19") {
    # Read BED file
    regions <- bed_df
    colnames(regions)[1:3] <- c("chr", "start", "end")
    
    # Create GRanges object
    gr <- makeGRangesFromDataFrame(regions, keep.extra.columns=TRUE)
    
    # Get genome sequence (change according to your genome)
    genome_seq <- getBSgenome(paste0("BSgenome.Hsapiens.UCSC.", genome))
    
    # Get sequences for each region
    sequences <- getSeq(genome_seq, gr)
    
    # Count motif motifs (motif dinucleotides)
    motif_counts <- vcountPattern(motif, sequences)
    
    # Calculate density (motif per base)
    region_lengths <- width(gr)
    motif_density <- motif_counts / region_lengths
    
    # Add results to original data
    regions$motif_count <- motif_counts
    regions$motif_density <- motif_density
    
    return(regions)
}
