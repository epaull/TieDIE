#!/usr/bin/env	python


from collections import defaultdict

for maf in open('mafs.txt', 'r'):

	data_by_sample = defaultdict(set)

	maf = maf.rstrip()
	tissue, mafFile = maf.split('\t')
	print ("processing "+tissue)

	entrez_gene_ids = set()	
	fh = open(mafFile, 'r')
	header = None
	for line in fh:
		parts = line.rstrip().split('\t')
		if not header:
			header = parts
			continue

		data = dict(zip(header, parts))
		if data["Variant_Classification"] != "Silent":
			barcode = data["Tumor_Sample_Barcode"]
			sample_id = '-'.join(barcode.split('-')[0:4])

			data_by_sample[sample_id].add(data["Entrez_Gene_Id"])
			# add data  
			entrez_gene_ids.add(data["Entrez_Gene_Id"])

	fh.close()

	
	## print out to matrix
	outfile = tissue+'/mutations.tab'
	fh = open(outfile, 'w')	

	sample_order = sorted(data_by_sample.keys())
	print ("number of samples: "+str(len(sample_order)) )
	fh.write('Entrez_ID\t'+'\t'.join(sample_order)+'\n')

	for gene in entrez_gene_ids:

		printstr = gene
		for sample in sample_order:
			if gene in data_by_sample[sample]:
				printstr += '\t1'
			else:
				printstr += '\t0'
		fh.write(printstr+'\n')
	
	fh.close()

