# protein IBA

  These scripts are an extension of a pipeline to process Anchord Hybrid enrichment (AHE) devolped in Breinholt et al. 2018 (Resolving Relationships among the Megadiverse Butterflies and Moths with a Novel Pipeline for Anchored Phylogenomics. Systematic Biology, Volume 67, Issue 1, 1 January 2018, Pages 78–93, https://doi.org/10.1093/sysbio/syx048). The core set of scripts and for this pipeline are avalible on dryad https://datadryad.org/resource/doi:10.5061/dryad.rf7g5.2. These scripts extend the abillity of IBA to use amino acid instead of nucleotide references sequences and were developed for Wolfe et al. 2018 in review to allow identification of orthologs from a set of more distantly related taxa by using aminio acid for the refrence squnces used as baits for assembly.  

>## genome_getprobe_TBLASTX.py
>  
>This is a modification of the script (genome_getprobe_BLAST.py) from  Espeland et al. 2018. to use tblastx instead of blastn to screen loci from denovo assemblies of capture data against set of reference coding exons loci. For use "python genome_getprobe_TBLASTX.py -h" for script help menu and defination of needed input.
>  

>## protein_IBA.py

>This is a modification of the script (IBA.py) from  Breinholt et al. 2018. to use protein refrence sequences as baits for IBA assembly. For use "protein_IBA.py -h" for script help menu and defination of needed input. Script contains detailed commented header with python, modules and dependancy requierments.
> 

>## protein_dir_fixer.py
 
>This script takes the output of protein_IBA.py (*_finaltable.table and *_finalseqs.fasta) and puts sequences in the same direction of the protein reference and produces a sequnce file with data trimed to blast hit coordinates from the reference 
For use "python protein_dir_fixer.py -h" for script help menu and defination of needed input. Script contains detailed commented header with python, modules and dependancy requierments.
>
