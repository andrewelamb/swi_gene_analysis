library(gridExtra)
library(grid)

matrix <- round(read.table("CO_ME_results/pw_all_matrix.tsv"), 2)

jpeg("images/co_vs_me_tbl.jpg", quality = 100, width = 1200, height = 900)
grid.table(matrix)
dev.off()
