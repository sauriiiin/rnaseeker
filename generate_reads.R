args <- commandArgs(trailingOnly = TRUE)
options(warn=-1)

# Get input arguments
input_dir <- args[1]
output_dir <- args[2]


##### Gather directories with samlon quants and read raw output
dirs<-list.dirs(path = sprintf("%s/salmon_quant",input_dir), full.names = TRUE, recursive = FALSE)
temp<-read.delim(paste(dirs[1],'quant.sf',sep = '/'),stringsAsFactors = FALSE)
raw_counts <- temp[,c(1,5)]

##### Compile raw count data.frame
for (i in seq(2,length(dirs))) {
  temp<-read.delim(paste(dirs[i],'quant.sf',sep = '/'),stringsAsFactors = FALSE)
  raw_counts <- merge(raw_counts, temp[,c(1,5)], by = 'Name', all = T)
}

##### Update column names with file names from foo.bar
sample_names <- readLines(sprintf('%s/foo.bar',input_dir))
colnames(raw_counts) <- c('ORF', sample_names)

# Write counts table to CSV file
write.csv(file = sprintf('%s/raw_counts.csv',output_dir), raw_counts, row.names = FALSE)
