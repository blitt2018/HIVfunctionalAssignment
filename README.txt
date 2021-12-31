All of the scripts and files needed to blast HIV proteins and get their enrichment scores are present in this folder. The scripts and files to analyze the enriched annotations are in the folder "analysis". Path names may need to be changed to make the pipeline runnable; however, nothing else has been modified about the pipeline that would make it not runnable.  

HIVProteinsFirstHalf: The first half of the HIV proteins. Proteins were split in this way so that two BLAST scripts could be run at once.

HIVProteinsSecondHalf: The second half of the HIV proteins. 

blastMC30_firstHalf_BLASTFIRST.sh: BLASTs first half with Protsub matrix and stores output in rawOutput. Must be run before blastMC30_secondHalf.sh
blastB62_firstHalf.sh: Same but for BLOSUM62 matrix. 

blastMC30_secondHalf.sh: BLASTs second half with Protsub matrix and stores output in rawOutput. 
blastB62_firstHalf.sh: Same but for BLOSUM62 matrix. 

swissprot folder: This folder is the database used when BLASTing and was made using the fasta file swissprot.fasta along with the makeblastdb tool.  

swissprot_ID_GOs.tab: This tsv file contains the GO annotations in swissprot for each protein. It is used as input to the uniIdsToGOsSwissprot.py file.   

makeEvalCut.py: Takes the raw blast output and only keeps hits under the specified evalue. 

uniIDsToGOsSwissprot.py: Takes the evalue-cutoff BLAST output and finds the Gene Ontology terms associated with each output protein.

applyFisherTestSwissprot.py: Takes the annotated proteins as input and outputs the results of applying fishers exact test on each of the annotations. The last column of the output file is the bonferroni-corrected p-value.

Uni_GO_Fisher.sh: A bash file to run all of the python files above in succession. All of the input and output files for each script are listed in this file. 



 

