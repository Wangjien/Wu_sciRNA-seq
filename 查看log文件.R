library(dplyr)
library(ggplot2)
library(stringr)

## merge log file
flist = list()
for (file in list.files(path='./', pattern = 'log.txt$')) {
#    print(file)
   tmp_name = str_split_fixed(file, '_log', n= 2)[,1]
   print(tmp_name)
   tmp = read.table(file, sep = '\t', fill = T, header=F)
   tmp$V1 = trimws(tmp$V1)
   tmp$name = tmp_name
   flist[[tmp_name]] = tmp
}