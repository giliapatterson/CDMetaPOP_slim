# Function to extract loci, alleles, and frequencies
split_genes <- function(gene_file, PatchID){
  genes <- read_csv(gene_file, show_col_types = FALSE)
  loci_and_allele <- genes$`Allele List`
  locus = str_extract(loci_and_allele, "^[^A]+(?=A\\d+$)")
  allele = as.numeric(str_extract(loci_and_allele, "(?<=A)\\d+$"))
  genes <- mutate(genes, locus = locus, allele = allele, PatchID = PatchID)
  return(genes)
}

# Function to merge two gene files
merge_genes <- function(genes1, genes2){
  patchID2 <- genes2$PatchID[1]
  new_df <- full_join(genes1, genes2, by = c("locus", "allele", "Allele List"),
                      suffix = c("", paste0("_", patchID2)))
  return(new_df)
}

patch_genes <- function(out_file, PatchID, all_genes){
  patch_genes <- select(all_genes,
                        "Allele List",
                        locus, allele, position,
                        paste0("Frequency_", PatchID),
                        paste0("PatchID_", PatchID)) |>
    rename(Frequency = paste0("Frequency_", PatchID))
  write_csv(patch_genes, out_file)
}

new_file_name <- function(new_directory, old_names){
  new_names <- rep("", length(old_names))
  for (i in seq_along(old_names)){
    old_name <- old_names[i]
    old_name_parts <- str_split_1(old_name, "/")
    new_names[i] <- paste0(new_directory, "/", old_name_parts[length(old_name_parts)])
  }
  return(new_names)
}



