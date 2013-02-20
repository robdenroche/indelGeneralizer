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



README
======
The genomic position of an indel is ambiguous when the bases being inserted or deleted are in the context of a homopolymer or simple repeat.  Since indel variants are often reported as occuring at a specific genomic position, an ambiguous indel could be called correctly at multiple positions due to differences in pipelines or conventions, and then erroneously considered to be separate mutation events.

The indel generalizer is a perl script capable of determining all the possible ambiguous alignments for a given indel and then generalizing the alignment to a user selected convention (the leftmost position in the reference, for example). It can be run in two modes: interactive mode where a single indel is supplied on the command line, or vcf mode where a vcf formatted file of variants can be processed, converting all ambiguous indels to a single convention.


INSTALLATION AND REFERENCES
The indelGeneralizer requires a bioperl module in order to index and pull bases from your fasta reference. Bioperl is available here: www.bioperl.org

Naturally you will also require a fasta reference, many of which can be obtained from UCSC here: http://hgdownload.soe.ucsc.edu/downloads.html .  Due to the current implementation, your fasta reference must have one file per chromosome all placed in the same directory. The format of your chromosome names in your input must match the names of your files (if you use chr12, you'll need a chr12.fa; for 12 you'll need 12.fa).

Once you have bioperl and a proper fasta reference you should be able to run the script as is! 

NOTE: The first time bioperl accesses a fasta file it will create an index file in the same directory as the fasta. As a result, the initial run of the script may take longer than you expect. If you happen to abort the script during index generation, future runs of the tool will attempt to use the incomplete index file and fail with a bioperl error. To resolve this, simply delete the incomplete index file and rerun the indelGeneralizer script.


USAGE


EXAMPLES


COMPUTE RESOURCES AND SORTING
When processing a sorted VCF file, the indelGeneralizer is designed to stream over the variants in order to minimize the run time and memory required. To reduce the number of reads operations from the fasta (which we have observed are costly) large chunks of sequence are read into memory at a time. Whenever a base is required by the script and isn't already in memory, the surrounding chunk is loaded. To reclaim memory, once the (presumed sorted) input has moved on to the next chunk, the old chunk is deleted. The algorithm still works on unsorted input, but will be less efficient. Unsorted input is discussed below.

By default the chunk size is 10^5, or 100,000 bases. In our tests, reading a 100,000 base chunk only takes twice as long as reading 1 base, while reading 100,000 bases individually takes 25 times as long as reading 1 base. This may not be the case in your environment, and tuning the chunk size may increase performance.

Because we are likely changing the reported position of some indels, it is possible that the sorted output will have a different order than the sorted input. A similar chunking scheme is used to locally sort the output; variants are placed into blocks once they have been generalized and once the input has moved sufficiently far past a block, its variants are sorted and printed. NOTE: in the unlikely event that a variant is generalized to a position that is further away than the size of the sorting blocks (currently hard coded as 10,000 bp) then it is possible that the output will not be completely sorted. A warning will be written to stderr if this is the case, and you will need to re-sort your vcf.


UNSORTED VCF
By default we assume the input is sorted - if this is not the case the -o flag can be used to skip local resorting efforts and preserve the order of the input. If the input is not sorted, reducing the size of the sequence chunks read from the reference may improve run time.


BUGS
This script has been tested but, like all code, it is possible bugs still exist. If you notice an issue, please let the author know.
