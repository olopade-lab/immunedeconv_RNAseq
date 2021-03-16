# Immune Deconvolution RNAseq  

## Gene name conversion:  

### Data:   

#### CIBERSORTx:  
Source: https://cibersortx.stanford.edu/download.php  
Required genes: see **`LM22.txt`**  

#### EPIC:
Source: https://github.com/GfellerLab/EPIC  
Required genes: see **`data/TRef.rda`**  
Note: convert to csv file with `write.csv(TRef[["sigGenes"]], "TRef_sigGenes.csv", row.names = FALSE)`  

#### MCP-counter: 
Source: https://github.com/ebecht/MCPcounter  
Required genes: see **`Signatures/genes.txt`**  

#### quanTIseq:
Source: https://github.com/icbi-lab/quanTIseq  
Required genes: see **`quantiseq/deconvolution/TIL10_signature.txt`**  
Note: a few genes are not recommended for RNA-seq, see **`quantiseq/deconvolution/TIL10_rmgenes.txt`**  

#### TIMER:
Source: http://cistrome.org/TIMER/download.html  
Required genes: see **`geneMarker.Rdata`** (data file 8: Gene markers to infer immune infiltration)  
Note: convert to csv file with `write.csv(geneMarker[["BRCA"]], "geneMarker_BRCA.csv", row.names = FALSE)`

#### xCell:
Source: https://github.com/dviraran/xCell  
Required genes: see **`genes_used_by_xCell.txt`**  

#### HGNCbiomart:
Source: https://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart  
Note: Select "Ensembl gene ID", "Approved symbol", "Previous symbol", "Alias Symbol". Output file is `HGNC_biomart/results.txt`  

### Scripts:

- **`HGNC_biomart/create_biomart_database.py`**: reformats the HGNCbiomart output to **`biomart_database.csv`**  
- **`required_genes/extract_required_genes.py`**: for each algorithm, it extracts the list of required immune genes for deconvolution. For quanTIseq, it also removes unwanted genes.  
- **`create_conversion_tables.py`**: for each algorithm, it creates the specific reference database to convert Ensembl IDs to HGNC symbols. First the algorithm's specific immune genes are converted, using the current symbols, the previous symbols and the alias symbols to 
- **`create_conversion_tables.py`**: for each algorithm
    - `-i` or 
