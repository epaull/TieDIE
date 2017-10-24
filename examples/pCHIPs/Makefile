
# make the mutation status matrix
data/ucec.tab:
	python scripts/makeMutMatrix.py

# genomic events with differential frequency in the sample groups
ucec/DIFF_A:
	Rscript scripts/ucec.genomics.R

ucec/upstream.input:
	python scripts/inputHeats.py ucec > $@

# run TieDIE, generate a scaffold network:
ucec/tiedie:
	python ../../bin/tiedie \
		-k pathways/multinet.kernel.obj \
		-n pathways/multinet.sif \
		ucec/upstream.input \
		ucec/downstream.input \
		-s 0.5 \
		--output ucec/TieDIE \
		--skip_nullModel

ucec/n1_network/n1_network.txt:
	python ../../bin/n1_pchip \
		ucec/n1.patient.input \
		ucec/downstream.input \
		ucec/TieDIE/TieDIE.sif \
		--output ucec/n1_network
