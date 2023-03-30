
rm(list = ls())
graphics.off()

if(!require(Biostrings)){
  install.packages("Biostrings")
  library(Biostrings)
}

args = commandArgs(trailingOnly = TRUE)

fastafile_path <- args[1]
# example: "C:\\Users\\JiaYing\\GP\\Rubus_Idaeus.genome.fa.gz"
outputfile_path <- args[2]
#example: "C:\\Users\\JiaYing\\GP\\Rubus_Idaeus_result.txt"

calculate_fasta_stats <- function(file_path) {
  # Read in the fasta file
 
  fasta_seqs <- readDNAStringSet(file_path)
  
  # Calculate the total number of sequences
  num_seqs <- length(fasta_seqs)
  
  # Calculate the length of each sequence
  seq_lengths <- width(fasta_seqs)
  
  # Calculate the total length of all sequences
  total_length <- sum(seq_lengths)
  
  # Calculate the average length of a sequence
  avg_length <- mean(seq_lengths)
  
  # Calculate the index of the longest and shortest sequence
  longest_seq_index <- which.max(seq_lengths)
  shortest_seq_index <- which.min(seq_lengths)
  
  # Extract the longest and shortest sequences
  longest_seq <- width(fasta_seqs[longest_seq_index])
  shortest_seq <- width(fasta_seqs[shortest_seq_index])
  
  # Calculate the N50
  sorted_lengths <- sort(seq_lengths, decreasing = TRUE)
  cumulative_lengths <- cumsum(sorted_lengths)
  n50_index <- which.max(cumulative_lengths >= 0.5 * total_length)
  n50 <- sorted_lengths[n50_index]
  
  # Calculate the number of N positions
  num_n <- sum(vcountPattern("N", fasta_seqs))
  
  # Calculate the GC content for each sequence
  gc_content <- sapply(fasta_seqs, function(x) {
    (countPattern("G", x) + countPattern("C", x)) / nchar(x)
  })
  
  # Calculate the mean GC content
  mean_gc_content <- ((mean(gc_content))*100)
  
  # Return a named list with the results
  return(list(Total_Number_Of_Sequences = num_seqs,
              avg_length = avg_length,
              shortest_seq = shortest_seq,
              longest_seq = longest_seq,
              n50 = n50,
              mean_gc_content = mean_gc_content,
              num_n = num_n))
}

result <- calculate_fasta_stats(fastafile_path)

fileConn<-file(outputfile_path)
writeLines(c(paste("Total number of sequences:", result[[1]], sep = "\t"),
             paste("Average length:", result[[2]], sep = "\t"),
             paste("Shortest sequence length:", result[[3]], sep = "\t"),
             paste("Longest sequence length:", result[[4]], sep = "\t"),
             paste("N50:", result[[5]], sep = "\t"),
             paste("Mean GC content:", result[[6]], sep = "\t"),
             paste("Number of Ns:", result[[7]], sep = "\t")), 
           fileConn)
close(fileConn)

##write.table(result, outputfile_path, sep = "\t", quote = FALSE)
