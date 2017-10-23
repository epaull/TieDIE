
library(reshape)
library(methods)

rootFolder <- "metadata"
hugo2entrez <- as.matrix(read.table(paste0(rootFolder, "/hgnc.map"), sep='\t', header=F))

map.entrez <- function(entrez.ids) {

	mapped <- map.entrez._internal(entrez.ids)
	if (all(is.na(mapped))) {
		return (map.entrez._internal(names(entrez.ids)))
	}
	mapped
}

map.entrez._internal <- function(entrez.ids) {
        mapped <- c()
        for (name in entrez.ids) {
                idx <- which(as.numeric(hugo2entrez[,2]) == as.numeric(name))
                if (length(idx) < 1) {
                        mapped <- c(mapped, NA)
                } else {
                        mapped <- c(mapped, hugo2entrez[idx, 1])
                }
        } 
        return (mapped)
}
map.hugo <- function(hugo.ids) {
        mapped <- c()
        for (name in hugo.ids) {
                idx <- which(hugo2entrez[,1] == name)
                if (length(idx) != 1) {
			# don't map to Entrez if it fails: just keep the HUGO id
                        mapped <- c(mapped, name)
                } else {
                        mapped <- c(mapped, as.numeric(hugo2entrez[idx, 2]))
                }
        } 
        return (mapped)
}


proportions <- function(data, samples.POS, samples.NEG, FDR=0.05, correct=T) {


	samples.POS <- intersect(colnames(data), samples.POS)
	samples.NEG <- intersect(colnames(data), samples.NEG)

	if (length(samples.POS)==0) { print ("Error, no positive class samples overlapping with gistic data!"); q(); }
	if (length(samples.NEG)==0) { print ("Error, no negative class samples overlapping with gistic data!"); q(); }

	# compute p-values using Fisher's exact test for differences in the proportions of 
	# high-level amp or deletion events. Get genomic locations for events.
	pvals <- sort(apply(data, 1, function(x) {

		eventsCount.A <- length(which(x[samples.POS] > 0))
	 	sampleSize.A <- length(samples.POS)
		eventsCount.B <- length(which(x[samples.NEG] > 0))
	 	sampleSize.B <- length(samples.NEG)
		
		test.res <- fisher.test( matrix(c(eventsCount.A, sampleSize.A, eventsCount.B, sampleSize.B), 
			ncol=2), alternative='greater')
		test.res$p.value
	}))
	sig.pvals <- NULL
	if (correct) {
		sig.pvals <- pvals[which(p.adjust(pvals, method='BY') < FDR)]
	} else {
		sig.pvals <- pvals[which(pvals < 0.05)]
	}

	# return gene names
	return (as.character(map.entrez(names(sig.pvals))))
}

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
