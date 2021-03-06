indelGeneralizer
================

Created by Rob Denroche on 2012-10-21
Copyright (C) 2012 The Ontario Institute for Cancer Research

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.




Overview
--------
The genomic position of an indel is ambiguous when the bases being inserted or deleted are in the context of a homopolymer or simple repeat.  Since indel variants are often reported as occuring at a specific genomic position, an ambiguous indel could be called correctly at multiple positions due to differences in pipelines or conventions, and then erroneously considered to be separate mutation events.

The indel generalizer is a perl script capable of determining all the possible ambiguous alignments for a given indel and then generalizing the alignment to a user selected convention (the leftmost position in the reference, for example). It can be run in two modes: interactive mode where a single indel is supplied on the command line, or vcf mode where a vcf formatted file of variants can be processed, converting all ambiguous indels to a single format.


Installation and References
---------------------------
The indelGeneralizer requires a bioperl module in order to index and pull bases from your fasta reference. Bioperl is available here: www.bioperl.org.

Naturally you will also require a fasta reference, many of which can be obtained from UCSC here: http://hgdownload.soe.ucsc.edu/downloads.html.  Due to the current implementation, your fasta reference must have one file per chromosome all placed in the same directory. The format of your chromosome names in your input must match the names of your files (if you use chr12, you'll need a chr12.fa; for 12 you'll need 12.fa).

Once you have bioperl and a proper fasta reference you should be able to run the script as is! 

NOTE: The first time bioperl accesses a fasta file it will create an index file in the same directory as the fasta. As a result, the initial run of the script may take longer than you expect. If you happen to abort the script during index generation, future runs of the tool will attempt to use the incomplete index file and fail with a bioperl error. To resolve this, simply delete the incomplete index file and rerun the indelGeneralizer script.


Usage
-----
Usage is ./indelGeneralizer.pl [-i INDEL -c CHROM -p POS || -v INPUT.VCF] -f /PATH/TO/FASTA/ [options]


### Options are as follows:

**-f /PATH/TO/FASTA/** is the path to your fasta reference, which is split into separate chromosome files.

**-a ANCHOR_POS** (as left, left-1, right, right+1). Optional, default is left-1.  Controls which position is reported for an indel.  For insertions, left and left-1 (or right and right+1) behave the same, reporting the position of the base of the left of the insertion point.  For deletions, left is the leftmost deleted base, and left-1 is the anchor base to the left of the deletion event.  Right is the rightmost base, and right+1 is the base after the deletion.

**-b ANCHOR_BASE** (as left, right, off).  Optional, default is left.  Toggles whether and where to include the anchor base.  To avoid confusion, "-a left-1 -b left" or "-a left -b off" are recommended.  Using right is likely to confuse.

**-s CHUNK_SIZE**  Optional, default is 5 (used as 10^CHUNK_SIZE to set the chunk size in bases).  Controls how many sequential bases are read from the fasta reference at a time (as I/O reads are costly).  Could be tuned lower if making multiple calls to interactive mode or if your VCF is unsorted, or higher if your I/O is slow and you have memory to spare.

**-h** displays this usage message.


### Interactive mode options:

**-i INDEL** (as +ACG or -T).  Enables interactive mode and specifies the inserted (+) or deleted (-) bases.

**-c CHROM** (as chr2 or 2 - whatever matches your fasta files).  Enables interactive mode and specifies the chromosome where the indel occurs.

**-p POS** (as 12345).  Enables interactive mode and specifies the position of the anchor base to the left to the indel event.

**-m ALIGN_MODE** (as left, right, full or all).  Optional, default is left.  Specifies the output of interactive mode.


### VCF mode options:

**-v INPUT.VCF** enables VCF mode and specifies the VCF file to generalize.

**-m ALIGN_MODE** (as left, right).  Optional, default is left.  Specifies whether to generalize indels to their leftmost or rightmost position.

**-o** Preserves existing VCF order.  By default (assuming the input is sorted) indels may be assigned new anchor positions and will be resorted accordingly to generate sorted output.  -o disables this sort check and the output will be in the same order as the input even if this violates the sort.


Examples
--------
```
# The four interactive output modes for a deletion of C at chr14:22298234.
./indelGeneralizer.pl -i -C -c chr14 -p 22298234 -f /references/fasta/ -m left
chr14   22298234        .       TC      T


./indelGeneralizer.pl -i -C -c chr14 -p 22298234 -f /references/fasta/ -m right
chr14   22298237        .       CC      C


./indelGeneralizer.pl -i -C -c chr14 -p 22298234 -f /references/fasta/ -m all
chr14   22298234        .       TC      T
chr14   22298235        .       CC      C
chr14   22298236        .       CC      C
chr14   22298237        .       CC      C


./indelGeneralizer.pl -i -C -c chr14 -p 22298234 -f /references/fasta/ -m full

                                       |22298232
(reference sequence)                   TCTCCCCAGC
chr14   22298234   .   TC   T          TCT-CCCAGC   <-- (original deletion)
chr14   22298235   .   CC   C          TCTC-CCAGC
chr14   22298236   .   CC   C          TCTCC-CAGC
chr14   22298237   .   CC   C          TCTCCC-AGC
(reference sequence)                   TCTCCCCAGC
                                        22298241|
                                        
                                        
# A three base insertion of AGA into chr8:48910891.
./indelGeneralizer.pl -i +AGA -c chr8 -p 48910891 -f /references/fasta/ -m full

                                        |48910880
(aligned consensus)                     CGTAGAAGAAGAAGAAGAGAC
chr8   48910882   .   T   TAGA          CGT***AGAAGAAGAAGAGAC
chr8   48910885   .   A   AAGA          CGTAGA***AGAAGAAGAGAC
chr8   48910888   .   A   AAGA          CGTAGAAGA***AGAAGAGAC
chr8   48910891   .   A   AAGA          CGTAGAAGAAGA***AGAGAC   <-- (original insertion)
chr8   48910894   .   A   AAGA          CGTAGAAGAAGAAGA***GAC
(aligned consensus)                     CGTAGAAGAAGAAGAAGAGAC
                                                    48910897|


# A five base deletion of ATTTA at chr17:41246043.  Note how the order of the deleted bases can change yet still result
# in the same observed sequence.
./indelGeneralizer.pl -i -ATTTA -c chr17 -p 41246043 -f /references/fasta/ -m full

                                           |41246037
(reference sequence)                       CGCTTTAATTTATTT
chr17   41246039   .   CTTTAA   C          CGC-----TTTATTT
chr17   41246040   .   TTTAAT   T          CGCT-----TTATTT
chr17   41246041   .   TTAATT   T          CGCTT-----TATTT
chr17   41246042   .   TAATTT   T          CGCTTT-----ATTT
chr17   41246043   .   AATTTA   A          CGCTTTA-----TTT   <-- (original deletion)
(reference sequence)                       CGCTTTAATTTATTT
                                                 41246051|


```



Compute Resources and Sorting
-----------------------------
When processing a sorted VCF file, the indelGeneralizer is designed to stream over the variants in order to minimize the run time and memory required. To reduce the number of reads operations from the fasta (which we have observed are costly) large chunks of sequence are read into memory at a time. Whenever a base is required by the script and isn't already in memory, the surrounding chunk is loaded. To reclaim memory, once the (presumed sorted) input has moved on to the next chunk, the old chunk is deleted. The algorithm still works on unsorted input, but will be less efficient. Unsorted input is discussed below.

By default the chunk size is 10^5, or 100,000 bases. In our tests, reading a 100,000 base chunk only takes twice as long as reading 1 base, while reading 100,000 bases individually takes 25 times as long as reading 1 base. This may not be the case in your environment, and tuning the chunk size may increase performance.

Because we are likely changing the reported position of some indels, it is possible that the sorted output will have a different order than the sorted input. A similar chunking scheme is used to locally sort the output; variants are placed into blocks once they have been generalized and once the input has moved sufficiently far past a block, its variants are sorted and printed. NOTE: in the unlikely event that a variant is generalized to a position that is further away than the size of the sorting blocks (currently hard coded as 10,000 bp) then it is possible that the output will not be completely sorted. A warning will be written to stderr if this is the case, and you will need to re-sort your vcf.


Unsorted VCF
------------
By default we assume the input is sorted - if this is not the case the -o flag can be used to skip local resorting efforts and preserve the order of the input. If the input is not sorted, reducing the size of the sequence chunks read from the reference may improve run time.


Bugs
----
This script has been tested but, like all code, it is possible bugs still exist. If you notice an issue, please let the author know.
