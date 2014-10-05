m = read.csv('cspn/temp/p-values.txt', sep="\t", header=FALSE)
m [,3] = p.adjust(m[,3], method='fdr')
write.table(m, file='cspn/temp/adj.p-values.txt', 
  sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
