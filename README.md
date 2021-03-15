# immunedeconv_RNAseq

CIBERSORTx:
Data: https://cibersortx.stanford.edu/download.php
Required genes: see LM22.txt
Example genes: from CIBERSORT, (https://cibersort.stanford.edu/download.php)


EPIC:
Data: https://github.com/GfellerLab/EPIC
Required genes: see data/TRef.rda
write.csv(TRef[["sigGenes"]], "TRef_sigGenes.csv", row.names = FALSE)
TRef_sigGenes.csv
Example genes: TRef_refProfiles.csv
write.csv(TRef[["refProfiles"]], "TRef_refProfiles.csv", row.names = FALSE)

MCP-counter: 
Data: https://github.com/ebecht/MCPcounter
Required genes: see Signatures/genes.txt
Example genes: ??

quanTIseq:
Data: https://github.com/icbi-lab/quanTIseq
Required genes: see quantiseq/deconvolution/TIL10_signature.txt
Remove genes see quantiseq/deconvolution/TIL10_rmgenes.txt
Example genes: example_gene_tpm.txt
https://icbi.i-med.ac.at/software/quantiseq/doc/sampledata/example_gene_tpm.txt

TIMER:
Data: http://cistrome.org/TIMER/download.html
Required genes: geneMarker.Rdata see  (data file 8: Gene markers to infer immune infiltration)
write.csv(geneMarker[["BRCA"]], "geneMarker_BRCA.csv", row.names = FALSE)
Example genes: see TIMER2.0?

xCell:
Data:
Required genes: see genes_used_by_xCell.txt
Example genes: 

Example genes:
Use TIMER2.0?

HGNC: HGNCbiomart
