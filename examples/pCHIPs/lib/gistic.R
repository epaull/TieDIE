
library(reshape)
library(methods)

rootFolder <- "/ifs/data/c2b2/ac_lab/shares/mappings"
if (!file.exists(rootFolder)) {
	rootFolder <- '/Volumes/pancancer/data/mappings'
}
hugo2entrez <- as.matrix(read.table(paste0(rootFolder, "/hgnc.map"), sep='\t', header=F))

# add interactions in intB to intA
# !!not symmetrical!!
mergeMRScores <- function(intA, intB) {

	res <- lapply(names(intA), function(index) {

		scoresA <- intA[[index]]
		new <- scoresA
		if (index %in% names(intB)) {
			if (length(intersect(names(intB[[index]]), names(scoresA)))>0) {
				print ('Error: name collision in mergeMRScores()')
				print (intA[[index]])
				print (intB[[index]])
				q();
			}
			new <- c(scoresA, intB[[index]])	
		}
		new
	})
	names(res) <- names(intA)
	res
}


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


gisticResult <- setClass("gisticResult", slots=
		c(data.by.sample="data.frame", 
		data.by.band="matrix", 
		genes.to.band="list", 
		band.to.genes="list", 
		locus="numeric", 
		sig.dels="data.frame",
		sig.amps="data.frame",
		cytoband="character"))

# Register class method
mapScores <- function(gisticObj, diggit.interactions, mapping='cytoband', from.p=FALSE, pos.nes.only=TRUE) {
	UseMethod("mapScores", gisticObj)
}

# convert diggit interactions to summaries based on genomic location
mapScores.gisticResult <- function(gisticObj, diggit.interactions, mapping='cytoband', from.p=FALSE, pos.nes.only=TRUE) {

	# apply over each TF/MR:
	mapped.scores <- lapply(diggit.interactions, function(x) {

		print (paste("Integrating: ", length(x), " interactions "))
		# check for zero interactions case
		if (length(x) == 0) {
			return (0)
		}	

		scores <- NULL
		# corner case of 1 interaction only...
		if (length(x) == 1) {
			score <- 0
			if (from.p) {
				score <- qnorm(min(x), lower.tail=F)
			} else {
				if (pos.nes.only) {
					# Only include positive scores
					score <- x
					if (score < 0) { score <- 0 }
				} else {
					# Only include positive scores
					score <- x
				}
			}
			names(score) <- names(x)
			scores <- score
		} else {

			# gather either cytoband locations or Locus
			# names are the genes in ENTREZ space
			# 
			locations <- NULL
			if (is.null(mapping)) {
				# do no mapping, use the gene names
				locations <- as.character(names(x))
				names(locations) <- as.character(names(x))
			} else if (mapping=='cytoband') {
				locations <- na.omit(gisticObj@cytoband[as.character(names(x))])
			} else if (mapping=='location') {
				locations <- na.omit(gisticObj@locus[as.character(names(x))])
			}  else {
				# exception
				print ('Error: no valid mapping value supplied to map aggregate function!')
			}

			if (length(locations) == 0) { 
				print (paste("Warning: didn't find any genomic locations for", names(x)))
				return (0)
			}
		
			# find the highest MR-gene score in this location and use that as a summary
			scores <- unlist(lapply(unique(locations), function(loc) {
				# for each location, get entrez gene IDs
				genes.this.loc <- as.character(names(locations[which(locations == loc)]))
				# determine if we're using z-scores or p-values here
				scores.thisLoc <- x[genes.this.loc]
				if (length(scores.thisLoc) == 0) { return (NA) }
				score <- NULL
				if (from.p) {
					score <- qnorm(min(scores.thisLoc), lower.tail=F)
				} else {
					if (pos.nes.only) {
						# Only include positive scores
						score <- as.numeric(max(scores.thisLoc))
						if (score < 0) { score <- 0 }
					} else {
						# take the max abs value of the score
						score <- scores.thisLoc[which.max(abs(scores.thisLoc))]
					}
				}
				score
			}))
			names(scores) <- unique(locations)
		}

		# add fusion events if necessary
		scores
	})
	names(mapped.scores) <- names(diggit.interactions)
	mapped.scores
}

#'
#' Filter TF-genomic event interactions to include only the highest scoring 
#' interaction in a given locus or cytoband. 
#'
#'
#' @param diggit.interaction.pvals : list indexed by TF with associated events and p-values
#' vector each
#' @param 
#' 

# Register the method to this class here
cnvEventsGistic <- function(gisticObj, sample, cytoband=TRUE) {
	UseMethod("cnvEventsGistic", gisticObj)
}
# class method definition: 
#'
#' events : a set of events in gene-space to cover. These are entrez IDs so will need to be mapped. 
cnvEventsGistic.gisticResult <- function(gisticObj, sample, cytoband=TRUE, from.p=FALSE, pos.nes.only=TRUE, threshold=0.5) {

	# gather either cytoband locations or Locus
	locations <- NULL
	if (cytoband) {
		locations <- na.omit(gisticObj@cytoband)
	} else {
		locations <- na.omit(gisticObj@locus)
	}

	## slice the matrix with this sample 
	slice <- gisticObj@data.by.sample[,sample]
	# get a vector of these events, for this sample
	names(slice) <- rownames(gisticObj@data.by.sample)

	## now determine which of these is present...
	dels <- list()
	amps <- list()
	for (loc in unique(locations)) {
		# for all genes at this location, take the max score
		genes.this.loc <- as.character(names(locations[which(locations == loc)]))
		# amps / deletions to -2 or +2 calls
		has.amp <- any(slice[genes.this.loc] > threshold)
		has.del <- any(slice[genes.this.loc] < -threshold)

		# if we get contradictory evidence, skip this location entirely
		if (has.amp & has.del) {
			next
		}
		else if (has.amp) {
			amps[[as.character(loc)]] <- names(slice[genes.this.loc][which(slice[genes.this.loc] > threshold)])
		}
		else if (has.del) {
			dels[[as.character(loc)]] <- names(slice[genes.this.loc][which(slice[genes.this.loc] < -threshold)])
		}
	}

	list(amps=amps, dels=dels)
}

# Register the method to this class here
mapGenes2Band <- function(gisticObj, genes) {
	UseMethod("mapGenes2Band", gisticObj)
}
# class method definition: 
#' diggit.interactions: z-scores 
mapGenes2Band.gisticResult <- function(gisticObj, genes) {
	bands <- unique(unlist(lapply(genes, function(gene) gisticObj@genes.to.band[[gene]])))
	bands
}
	# apply over each TF/MR:
# Register the method to this class here
cnvScoreStouffer <- function(gisticObj, diggit.interaction, cytoband=TRUE, from.p=FALSE, pos.nes.only=TRUE) {
	UseMethod("cnvScoreStouffer", gisticObj)
}

# class method definition: 
#' diggit.interactions: list indexed by MR/TF name in Entrez Space
#' 	each points to a named vector of NES / z-scores associated with entrez IDs for each interacting event. 
cnvScoreStouffer.gisticResult <- function(gisticObj, diggit.interactions, cytoband=TRUE, from.p=FALSE, pos.nes.only=TRUE) {

	if (cytoband) {
		mapping <- 'cytoband'
	} else {
		mapping <- 'location'
	}

	mapped.diggit <- mapScores(gisticObj, diggit.interactions, mapping=mapping, from.p=FALSE, pos.nes.only=TRUE)

	# apply over each TF/MR: sum the named vector of scores 
	integrated.z.scores <- unlist(lapply(mapped.diggit, function(scores) {

		# now integrate scores with Stouffer's method
		integrated.z <- sum(scores)/sqrt(length(scores))
                if (is.null(integrated.z)) { integrated.z <- 0 }
                if (is.nan(integrated.z)) { integrated.z <- 0 }
		integrated.z
	}))

	names(integrated.z.scores) <- names(diggit.interactions)
	integrated.z.scores
}

#' Parse binary matrix into a dataframe. Then 
#' @param : file to the 'all_thresholded_by_gene.txt' from gistic 
#' includes multinomial scores (-2,-1, 0, 1, 2) for high/low level amp
#' and deletion events. 
#' @param : gene.filter is a vector in HUGO space, remove any gene level data not with these ids
parse.gistic <- function(events.file, sig.amp.file, sig.del.file, gene.filter=NULL) {

	##
	## Parse thresholded data matrix, set metadata 
	##
	thresh.by.gene <- read.table(events.file, header=T, sep='\t', row.names=1, check.names=F)
	print (dim(thresh.by.gene))
	# the first two columns are metadata, just process the data here
	short.sample.ids <- unlist(lapply(colnames(thresh.by.gene)[3:ncol(thresh.by.gene)], 
		function(sample.id) 
		paste(strsplit(sample.id, '-')[[1]][1:4], collapse='-')))

	if (length(short.sample.ids) != length(unique(short.sample.ids)) ) {
		print ("Error: short sample ids not unique!")
		return (NA)
	}
	# one more step: remove the A/B/C at the end, hopefully this is unique...otherwise we'll have errors
	# later on
	short.sample.ids <- unlist(lapply(short.sample.ids, function(x) substr(x, 1, 15)))

	# add the two metadata columns and fix the sample ids to be 16 digits
	colnames(thresh.by.gene) <- c(colnames(thresh.by.gene)[1:2], short.sample.ids)

	##
	## map to entrez IDs
	##
	rows.entrez <- map.hugo(rownames(thresh.by.gene))
	rownames(thresh.by.gene)[which(!is.na(rows.entrez))] <- rows.entrez[which(!is.na(rows.entrez))]
	# filter by fCNV
	if (!is.null(gene.filter)) {
		thresh.by.gene <- thresh.by.gene[na.omit(match(gene.filter, rownames(thresh.by.gene))),]
	}

	locus <- as.numeric(thresh.by.gene[,1])
	names(locus) <- rownames(thresh.by.gene)

	cytoband <- as.character(thresh.by.gene[,2])
	names(cytoband) <- rownames(thresh.by.gene)

	##
	## Parse amp and deletions files, find significant
	##
	del.mat <- read.table(sig.del.file, sep='\t', header=T)
	dels.sig <- melt(del.mat)
	amp.mat <- read.table(sig.amp.file, sep='\t', header=T)
	amps.sig <- melt(amp.mat)

	##
	## Collapse down to band-based representation   
	##
	unique.bands <- unique(as.character(thresh.by.gene[,2]))
	same.vals <- unlist(lapply(unique.bands, function(band) {
		dat <- thresh.by.gene[which(as.character(thresh.by.gene[,2]) == unique.bands[band]),]
		dat <- dat[,3:ncol(dat)]
		res <- apply(dat, 2, function(x) length(unique(x)) )
		(length(unique(res))==1)
	}))
	if (!all(same.vals)) {
		print ("Error: different gene-level calls for one or more locations!")
		q();
	}

	# data in band format: calls by band, since they're all the same	
	# this removes most of the space

	# this removes most of the space
	data.by.band <- matrix(unlist(lapply(1:length(unique.bands), function(band) {
      		row <- thresh.by.gene[which(as.character(thresh.by.gene[,2]) == unique.bands[band])[1],3:ncol(thresh.by.gene)]
		row
	})), byrow=T, nrow=length(unique.bands))
	rownames(data.by.band) <- unique.bands
	colnames(data.by.band) <- colnames(thresh.by.gene)[3:ncol(thresh.by.gene)]
	# and a list mapping band names to gene names
	band.to.genes <- lapply(1:length(unique.bands), function(band) {
      		rownames(thresh.by.gene)[which(as.character(thresh.by.gene[,2]) == unique.bands[band])]
	}) 
	names(band.to.genes) <- unique.bands   

	# and a list mapping genes to band ID
	genes.to.band <- lapply(1:nrow(thresh.by.gene), function(i) { as.character(thresh.by.gene[i,2]) })
	names(genes.to.band) <- as.character(rownames(thresh.by.gene))

	res <- gisticResult(data.by.sample=thresh.by.gene[,3:ncol(thresh.by.gene)], 
		data.by.band=data.by.band[,3:ncol(data.by.band)],
		band.to.genes=band.to.genes,
		genes.to.band=genes.to.band,
		locus=locus, 
		sig.dels=dels.sig,
		sig.amps=amps.sig,
		cytoband=cytoband)
}

## FIXME: need acessors to get datamatrix, subset based on confidence or not, collapse based on cytoband. 
