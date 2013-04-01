#!/usr/bin/env Rscript

# SPAN.R: 
# Runs the Fast-Heinz Heuristic algorithm to connect a set of 
# nodes over a network. A minimum spanning tree will be constructed using
# an optimal subset of the nodes provided and any other required 'linker' nodes
# in the network. 
#
# This is the computationally cheap version of the PCST algorithm, which is not
# guaranteed to converge in a reasonable amount of time (or even give an estimate of
# the time required to converge).
#
# Required Packages: BioNet
#
# Author: Evan Paull 
# Date: Oct, 2011

options(warn = -1)
library('getopt')

opt = getopt(matrix(c(
    'help' , 'h', 1, "character",
    'network' , 'n', 1, "character",
    'activities' , 'a', 1, "character",
    'pcst' , 'p', 1, "character",
    'fdr' , 'f', 1, "character"
	),ncol=4,byrow=TRUE));

if (!is.null(opt$help) || is.null(opt$network) || is.null(opt$activities)) {
	self = commandArgs()[1];
	#print a friendly message and exit with a non-zero error code
	cat(paste("Usage: ",self,"  --network <network file> --activities <activities (.lst)> file \n"))
	q();
}

library(BioNet)
library(igraph)

# PCST_RUN_FAST - Take a igraph network and a set of node scores and find
# a given number of suboptimal solutions using the PCST Fast Heinz heuristic 
# from BioNet package R libraries (http://bionet.bioapps.biozentrum.uni-wuerzburg.de/)
# 
# INPUT: 
# 	- subnet : an igraph network of named nodes
# 	- scores : a named vector of node scores corresponding to a subset of the nodes in (subnet)
#	- no_solutions : the number of suboptimal solutions to find (0 for only the optimal solution)
#
# RETURNS:
# 	- modules : a list object containing each solution as a graphNEL object
pcstRunFast <- function(subnet, scores, no_solutions) {

	modules <- list()
	runFastHeinz(subnet, scores)
	modules <- list()
	for (k in c(1:no_solutions)) {
		modules[[k]] <- runFastHeinz(subnet, scores)
	}

	return (modules)
}

# PCST_RUN - Take a igraph network and a set of node scores and find
# a given number of suboptimal solutions using the PCST Heinz algorithm 
# implemented by Ivana Ljubic (http://homepage.univie.ac.at/ivana.ljubic/research/pcstp/) 
# along with the BioNet package R libraries (http://bionet.bioapps.biozentrum.uni-wuerzburg.de/)
# 
# INPUT: 
# 	- subnet : an igraph network of named nodes
# 	- scores : a named vector of node scores corresponding to a subset of the nodes in (subnet)
#	- no_solutions : the number of suboptimal solutions to find (0 for only the optimal solution)
#
# RETURNS:
# 	- modules : a list object containing each solution as a graphNEL object
pcstRun <- function(subnet, scores, no_solutions) {

	# create temp directories for heinz files
	dir.create("heinztmp")
	dir.create("tmp")
	heinz.node.file = "tmp/heinz.n.file.txt"
	heinz.edge.file = "tmp/heinz.e.file.txt"


	writeHeinzEdges(network=subnet, file=heinz.edge.file, use.score=TRUE)
	writeHeinzNodes(network=subnet, file=heinz.node.file, node.scores=scores)

	runHeinz(heinz.folder="/projects/sysbio/apps/x86_64/heinz_1.63", heinz.e.file=heinz.edge.file, heinz.n.file = heinz.node.file, N=TRUE, E=FALSE, diff=-1, n=no_solutions)
	out.node.file <- paste(heinz.node.file,".0.hnz",collapse="",sep="")

	solutions <- list.files("./tmp",pattern="heinz.n.*.hnz")
	modules <- list()
	k <- 1
	for (s in sort(solutions)) {
		print(paste("reading solution: ", s))
		modules[[k]] <- readHeinzGraph(node.file = paste("./tmp",s,collapse="",sep="/"), subnet, format="igraph")
		k <- k + 1 
	}

	return (modules)
}

# MAP_NAMES_2_INTERNAL_IDS - convert a set of Ds
# into gene names, using a map matrix object
#
# INPUT: 
#	- names : a vector of gene ids 
#	- map : a 2-column matrix with gene ids
#	in column 1 paired with gene names in column 2 
#
# RETURNS: 
#	- a vector of the corresponding gene names (or
#	just the original IDs if none match)
mapModuleNames2Ids <- function(names, map) {

	ids <- vector(length=length(names))

	k <- 1
    for (i in c(1:length(names))) {

		# id + ":" type string to prepend to the gene name
		name <- names[i]

		# should get 3 values: the second is the gene name
		id <- map[which(map[,2] == name), 1][1]

		ids[k] <- id
		k <- k + 1
    }

	return (ids)
}

# SCORE: Score a network using the source and target
# genes used to generate it. 
#
# S = Source Genes ; Sr = Recovered Source Genes
# I = Total number of genes in the discovered network
# 
# SCORE = f_B(Sr/S) - 0.1*f_C(I-Sr)/S)
scoreSubnet <- function(activities, subnet) {

	net_genes <- V(subnet)$name
	Sr <- length(intersect(net_genes, activities))
	S <- length(activities)	
	I <- length(net_genes)

	#print(paste("Sr: ", Sr, collapse=""))
	#print(paste("S: ", S, collapse=""))
	#print(paste("(I-Sr)/S: ", (I-Sr)/S, collapse=""))
	return (Sr/S - 0.1*(I-Sr)/S)
}

# PARSE_ACTIVITIES_PVALS_FILE
#
# INPUT: 
#	- file : the relative path to a pvals activities file
#	<name> <ids> <type> <pvals>
#
# RETURNS:
#	- matrix of the ids, gene names and pvals 
#	for each tuple that has a p-value of less than 0.5
parseActivitiesPvalsFile <- function(file) {

	m <- as.matrix(read.delim(file, header=FALSE, sep="\t"))
	m <- m[which(as.numeric(c(m[,2])) < 0.49999),1:2]

	#gene_names <- c(m[,1])

	for (row in c(1:dim(m)[1])) {
		m[row,2] = gsub(" +","#", m[row,2])
		m[row,1] = gsub(" +","#", m[row,1])
		#m[row,1] = paste(gsub(" +","", m[row,1]),collapse="",sep=":")
	}
	ids <- c(m[,1])
	pvals <- as.numeric(c(m[,2]))

	m <- cbind(ids, pvals)

	return (m)
}

# PARSE_NETWORK_FILE - Read in a tab-separated network file list and 
# return an igraph object
#
#	rows: <ID_A> <ID_B> 
#   i.e.: HAP1  HGS   
# 
# INPUT:
#	- network_file : relative file path on the system
#	- weight_col : the (optional) column in the network file that contains
#		the edge weights (or zeros in the unweighted case)
#
# RETURNS:
#	- and igraph object using IDS:type
parseNetworkFile <- function(network_file, weight_col) {
	# input: a 2 column matrix of the edge lists
	m <- as.matrix(read.delim(network_file, sep="\t", header=FALSE))
	# change the names to add the type -- edges should be unique now
	for (row in c(1:dim(m)[1])) {
		# 4th column is the type for the second row
		#m[row,2] = paste(gsub(" +","#", m[row,2]), collapse="",sep=":")
		m[row,3] = gsub(" +","_", m[row,3])
		m[row,1] = gsub(" +","_", m[row,1])
	}
    g <- graph.edgelist(cbind(m[,1],m[,3]), directed=FALSE)
	# assign the weight vector: Heinz will automatically use the edge $score attribute
	# to assign weights
	if (weight_col > 0) {
		E(g)$score <- m[,weight_col]
	} else {
		E(g)$score <- rep.int(1, dim(m)[1])
	}

	return (g)
}

# GET_ROWS (UTILITY FUNCTION) - get a slice of one of the data columns
# of a matrix.
#
# INPUT: 
#	- names : the set of unique ids in the lookup column that are used
#	to select the subset of rows of the matrix
#	- m : the matrix object
# 	- lookup_col : the lookup column number
#	- data_col : the number of the data/values column to return as a 
#	vector
#
# RETURNS:
#	- a subset of data_col in the matrix, as a vector.
getRows <- function(names, m, lookup_col, data_col) {

	values <- vector(length=length(names))
	k <- 1
	for (name in names) {

		if (length(m[which(m[,lookup_col] == name),data_col]) > 0) {
			values[k] <- m[which(m[,lookup_col] == name),data_col]
		} else {
			values[k] = NA
		}
		k <- k + 1	

	}

	return (values)
}

# GET_PVALS_FROM_IDS - Take a set of gene IDS
# and return the associated set of experimental p-values
# 
# INPUT: 
#	- ids : a vector of Gene IDS
#	- pvals : a vector of p-values indexed by name
# 	of the associated  gene IDS
#
# RETURNS: 
#	- a named vector of p-values for the supplied subset
#	(or 1 for each id that was not found in the experiment pvals vector)
getPvalsFromIDS <- function(ids, pvals) {

	k <- 1

	pval_subset <- rep.int(1, length(ids))
	n <- names(pvals)

	for (i in c(1:length(ids))) {
		id <- ids[i]
		for (j in c(1:length(pvals))) {
			if (n[j] == id) {
				pval_subset[k] <- pvals[j]	
				k <- k + 1
				break
			}
		}
	}

	names(pval_subset) <- ids	

	return (pval_subset)
}

findModules <- function(network, activities_data) {

	pvals <- as.numeric(activities_data[,2])
	# assign gene ids
	names(pvals) <- activities_data[,1]
	
	
	# fit beta-uniform model and run pcst
	fb <- fitBumModel(pvals, plot = FALSE)
	# FIXME: how to test if this fit actually found signal?
	if (fb$negLL > -10) {
		print("ERROR: couldn't fit Beta-Uniform model!")
		q();
	}
	scores <- scoreNodes(network, fb, fdr = opt$fdr)
	scores[is.na(scores)] <- rep.int(0, length(is.na(scores)))
	if (is.null(opt$pcst)) {
		modules <- pcstRunFast(network, scores, 1)
	} else {
		modules <- pcstRun(network, scores, 1)
	}

	# return the single module
	return (modules[[1]])
}

plotPCSTModule <- function (network, layout = layout.fruchterman.reingold, labels = NULL, 
    diff.expr = NULL, scores = NULL, main = NULL, shapes = NULL, colors = NULL, ...) 
{
    if (is(network, "graphNEL")) {
        network <- igraph.from.graphNEL(network)
    }
    if (is.null(V(network)$name)) {
        V(network)$name <- as.character(V(network))
    }
    if (is.null(labels)) {
        if ("geneSymbol" %in% list.vertex.attributes(network)) {
            labels <- V(network)$geneSymbol
        }
        else {
            labels <- V(network)$name
        }
    }
    #shapes <- rep("circle", length(V(network)))
    names(shapes) <- V(network)$name
    #if (!is.null(scores) && !is.null(names(scores))) {
    #    shapes[intersect(names(which(scores < 0)), V(network)$name)] <- "csquare"
    #}
    #if (is.null(scores) && "score" %in% list.vertex.attributes(network)) {
    #    scores <- V(network)$score
    #    names(scores) <- V(network)$name
    #    shapes[names(which(scores < 0))] <- "csquare"
    #}
    #if (!is.null(diff.expr) && !is.null(names(diff.expr))) {
    #    coloring <- .node.color(network, diff.expr)
    #}

	coloring <- colors
	#diff.expr <- diff.expr / max(abs(as.numeric(diff.expr)))
	# go in order of the named nodes to build the coloring
	#for (j in c(1:length(V(network)$name))) {
	#	this_node <- V(network)$name[j]
	##	for (i in c(1:length(diff.expr))) {
	#		if (this_node == names(diff.expr)[i]) {
	#			print (this_node)
	#			print(diff.expr[i])
	#			if (diff.expr[i] > 0) {
    #				coloring <- c(coloring, "Pink")
	#			} else {
    ###				coloring <- c(coloring, "LightBlue")
	#			}
	#		}
	#	}
	#}

    if (is.null(diff.expr) && "diff.expr" %in% list.vertex.attributes(network)) {
        diff.exprs = V(network)$diff.expr
        names(diff.exprs) <- V(network)$name
        coloring <- .node.color(network, diff.exprs)
    }
    max.labels <- max(nchar(labels))
    network.size = length(V(network))
    cex = 0.6
    if (network.size < 50) {
        if (max.labels > 2) {
            vertex.size <- 8
            labels.dist <- 0.5
        }
        else {
            vertex.size <- 15
            labels.dist <- 0
        }
    }
    if (network.size < 100 && network.size >= 50) {
        vertex.size <- 8
        if (max.labels > 2) {
            labels.dist <- 0.5
        }
        else {
            labels.dist <- 0
        }
    }
    if (network.size >= 100) {
        vertex.size <- 8
        if (max.labels > 3) {
            labels.dist <- 0.5
        }
        else {
            labels.dist <- 0
        }
    }
    plot(network, layout = layout, vertex.size = vertex.size, 
        vertex.label = labels, vertex.label.cex = cex, vertex.label.dist = labels.dist, 
        vertex.color = coloring, vertex.label.family = "sans", 
        vertex.shape = shapes, main = main, ...)
}

createInternalMap <- function(names) {

	m <- matrix(ncol=2, nrow=length(names))
	m[,2] <- names
	m[,1] <- seq(1,length(names))

	return (m)
}

if (is.null(opt$fdr)) {
	opt$fdr = 0.05
} else {
	opt$fdr = as.numeric(opt$fdr)
}
# parse the PPI network into an undirected graph with uniform weights
write(paste("Parsing network file ", opt$network, "\n",collapse=""), stderr())
network <- parseNetworkFile(opt$network, 0)
write(paste("found network with ",  length(V(network)$name), " nodes and ", length(E(network)), " edges ", "\n", collapse=""), stderr())
# get mutations data matrix: 3 columns, with each tuple having
# (entrez gene ID, gene Name, experimentally dervied p-value of no-mutation NH)
write(paste("Parsing activities file ", opt$activies, "\n", collapse=""), stderr())
activities_data <- parseActivitiesPvalsFile(opt$activities)
#print(activities_data)
write(paste("Attempting to link up to ", dim(activities_data)[1], " genes over the network...\n"), stderr())

igraph_module <- findModules(network, activities_data)
score <- scoreSubnet(activities_data[,1], igraph_module)
#print (paste("Stat: ", "Score ", score, collapse=""))
#print(paste("Stat: ", "CC ",  transitivity(igraph_module, type="global"), collapse=""))
#print(paste("Stat: ", "Nodes ", length(V(igraph_module)$name), collapse=""))
#print(paste("Stat: ", "Edges ", length(E(igraph_module)), collapse=""))
#pdf(opt$pdf)
#cshapes <- rep("circle", length(V(igraph_module)$name))
#ccolors <- rep("White", length(V(igraph_module)$name))
#plotPCSTModule(igraph_module, diff.expr=V(igraph_module)$number, scores=V(igraph_module)$number, colors=ccolors, shapes=cshapes)
if (length(V(igraph_module)$name) == 1) {
	print (V(igraph_module)$name)
} else {
	print(E(igraph_module))
}
