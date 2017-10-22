
library(reshape)
library(methods)
library(mygene)

chasmResult <- setClass("chasmResult", slots=
		c(data.by.sample="data.frame", 
		sig.muts="character"))

parse.chasm <- function(events.file, ensembl.file=NULL) {

	library(mygene)

	ensembl.map <- NULL
	if (!is.null(ensembl.file)) {	
		ensembl.map <- read.table(ensembl.file, header=T, sep='\t')
	}
	
	data <- read.table(events.file, header=T, sep='\t')
	sig.data <- data[which(data$PValue < 0.05),]
	
	names <- as.character(apply(sig.data, 1, function(x) {
	   	name <- x[2]
		if (grepl('^N', name)) {
			name <- strsplit(name, '\\.')[[1]][1]
		} else if (grepl('^ENST', name)) {
			name <- strsplit(name, '_')[[1]][1]
			# try to map to ENTREZ Gene from Transcript
			name <- as.character(ensembl.map[which(ensembl.map[,2] == name),1])
			if (length(name)==0) { name <- NA }

		}
		name
	}))

	ids <- queryMany(names, scopes=c("refseq", "accession", "reporter", "ensembl.gene"), fields="entrezgene", species="human")
	
	sig.events <- as.character(ids$entrezgene)
	names(sig.events) <- as.character(substr(sig.data[,1], 1, 15))
	sig.events <- na.omit(sig.events)

	res <- chasmResult(
		data.by.sample=data,
		sig.muts=sig.events)

	res
}

## FIXME: need acessors to get datamatrix, subset based on confidence or not, collapse based on cytoband. 
