#Scripts for Tn-Seq 

This is a mostly random collection of python and R scripts for performming differential experession analysis on Tn-Seq experiments.  It originally started out for analysis from the Tufts Galaxy server and morphed into something more general.

#Basic Flow

This is the main use of the ```HopCounts.py``` script.  

1. Cleanup the reads using [flexbar](http://sourceforge.net/projects/flexbar/).  Previously I had been using fastx but Flexbar seemed a little more up to date and nicely parallelized.  Using Flexbar we clean up the reads for: quality; adapter sequences and length. Flexbar could be a pain to install because of its dependency on Intel's [Threaded Building Blocks](https://www.threadingbuildingblocks.org/).  These are just template headers so should compile for any server but this is untested.  For local work the binary downloads work out of the box with the appropriate ```LD_LIBRARY_PATH``` or ```DYLD_LIBRARY_PATH``` setup.  

```./flexbar -r ../CPR400bp_S1_L001_R1_001.fastq -f i1.8 --max-uncalled 2 --pre-trim-phred 30 --threads 2```


2. Map the reads using [bowtie2](https://github.com/BenLangmead/bowtie2).  I moved from bowtie to bowtie for no good reason other than we seemed to have longer reads which bowtie2 is supposed to be better at.  This mapping step needs a reference sequence to build an index for.  The reference sequence should be a fasta format and to prevent grief down the road make sure the sequence name matches the full annotation gff file used later.

```bowtie2-build PAO1.fasta PAO1```
```bowtie2 -x PAO1 -q --threads 2 -k 1 -U flexbar_v2.4_macosx/flexbar.fastq  -S test.mapped.sam```

3. Bin the mapped reads into features using [htseq](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html).  htseq is a python package (install with pip) that provides a nice python interface and a solid C backend for speed.   The stranded=reverse was necessary otherwise everything ended up as a __no_feature.  Again note that the reference sequence name (which acts as a chromosome) needs to match beween the gff file and reference used for mapping. 

```htseq-count --stranded=reverse --type=gene -i locus_tag test.mapped.sam PAO1.gff > test.count.tsv```


