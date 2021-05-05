library(Biostrings)

library(tidyr)
library(readr)

# Saves amino acid substitution matrices to file

# Writes the long format of the matrix to file
write_table <- function(mat, filename) {
	dataf <- as.data.frame(mat)
	dataf$Row <- rownames(dataf)
	long <- pivot_longer(dataf, cols=which(colnames(dataf) != "Row"))
	write_csv(long, file=filename, col_names=FALSE, quote_escape="none")
}

data(BLOSUM45)
data(BLOSUM50)
data(BLOSUM62)
data(BLOSUM80)
data(BLOSUM100)
data(PAM30)
data(PAM40)
data(PAM70)
data(PAM120)
data(PAM250)

matrix_names <- c(
"BLOSUM45",
"BLOSUM50",
"BLOSUM62",
"BLOSUM80",
"BLOSUM100",
"PAM30",
"PAM40",
"PAM70",
"PAM120",
"PAM250")

for (mat_name in matrix_names) {
	mat <- get(mat_name)
	filename <- file.path("Matrices", paste0(mat_name, ".csv"))
	write_table(mat, filename)
}

nt_subs <- nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")
write_table(nt_subs, file.path("Matrices", "Nucleotides.csv"))

