#!/usr/bin/perl

# Created by Rob Denroche on 2012-10-21
# Copyright (C) 2012 The Ontario Institute for Cancer Research
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
#      your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.


use strict;
use warnings;

use Getopt::Std;
use vars qw/ %opt /;

use Bio::DB::Fasta;


my $version = "1.0";

# defaults
my $alignMode = "left";
my $anchorPos = "left-1";
my $anchorMode = "left";
my $chunkPower = 5;
my $bufferChunkPower = 5;
my $outputVcfStats = 1;
my $doSort = 1;

sub usage
{
	warn "
 usage is ./indelGeneralizer.pl [-i INDEL -c CHROM -p POS || -v INPUT.VCF] -f /PATH/TO/FASTA/ [options]
 
 Options are as follows:
 	-f /PATH/TO/FASTA is the path to your fasta reference, which is split into separate chromosome files.
 	-a ANCHOR_POS (as left, left-1, right, right+1). Optional, default is left-1.  Controls which position is reported for an indel.
 		For insertions, left and left-1 (or right and right+1) behave the same, reporting the position of the base of the left of the insertion point.
 		For deletions, left is the leftmost deleted base, and left-1 is the anchor base to the left of the deletion event.  Right is the rightmost base,
 			and right+1 is the base after the deletion.
 	-b ANCHOR_BASE (as left, right, off).  Optional, default is left.  Toggles whether and where to include the anchor base.
 		To avoid confusion, \"-a left-1 -b left\" or \"-a left -b off\" are recommended.  Using right is likely to confuse.
 	-s CHUNK_SIZE (as 5).  Optional, default is 5 (used as 10^CHUNK_SIZE to set the chunk size in bases).  Controls how many sequential bases are read from 
 		the fasta reference at a time (as I/O reads are costly).  Could be tuned lower if making multiple calls to interactive mode or if your VCF is unsorted, or higher if
 		your I/O is slow and you have memory to spare.
 	-h displays this usage message.
 		
 Interactive mode options:
 	-i INDEL (as +ACG or -T).  Enables interactive mode and specifies the inserted (+) or deleted (-) bases.
 	-c CHROM (as chr2 or 2 - whatever matches your fasta files).  Enables interactive mode and specifies the chromosome
 		where the indel occurs.
 	-p POS (as 12345).  Enables interactive mode and specifies the position of the anchor base to the left to the indel event.
 	-m ALIGN_MODE (as left, right, full or all).  Optional, default is left.  Specifies the output of interactive mode.
 	
 VCF mode options:
 	-v INPUT.VCF enables VCF mode and specifies the VCF file to generalize.
 	-m ALIGN_MODE (as left, right).  Optional, default is left.  Specifies whether to generalize indels to their leftmost or rightmost position.
	-o Preserves existing VCF order.  By default (assuming the input is sorted) indels may be assigned new anchor positions and will be resorted accordingly to generate
		sorted output.  -o disables this sort check and the output will be in the same order as the input even if this violates the sort.
		
 	
 This is version $version of the indelGeneralizer.\n";

die "\n @_\n\n";
}


my $interactiveMode = 0;		# set to 1 if we're in interactive mode
my %validModes = ("left" => 1, "right" => 1, "full" => 1, "all" => 1);
my %validVcfModes = ("left" => 1, "right" => 1);
my %validAnchorPos = ("left" => 1, "right" => 1, "left-1" => 1, "right+1" => 1);
my %validAnchorBase = ("left" => 1, "right" => 1, "off" => 1);


my $type;
my $indel;
my $chr;
my $pos;

my %fastaHandles;

my $vcfPath;


###########################################################
# process command line options
###########################################################
my $opt_string = "f:s:i:c:p:m:v:a:b:oh";
getopts ($opt_string, \%opt) or usage("Incorrect arguments.");

if (exists $opt{h})
{
	usage("Help requested.");
}

if (exists $opt{f})
{
	if (-d $opt{f})
	{
		%fastaHandles = ("path" => $opt{f});
	}
	else
	{
		usage("$opt{f} is not a valid directory.");
	}
}
else
{
	usage("You must specify the path to your fasta reference files with -f.");
}

if (exists $opt{a})
{
	if (exists $validAnchorPos{$opt{a}})
	{
		$anchorPos = $opt{a};
	}
	else
	{
		usage("-a $opt{a} is not a valid ANCHOR_POS.");
	}
}

if (exists $opt{b})
{
	if (exists $validAnchorBase{$opt{b}})
	{
		$anchorMode = $opt{b};
	}
	else
	{
		usage("-b $opt{b} is not a valid ANCHOR_BASE.");
	}
}

if (exists $opt{s})
{
	if ($opt{s} =~ /^\d+$/)
	{
		$chunkPower = $opt{s};
	}
	else
	{
		usage("-s CHUNK_SIZE must be an integer.");
	}
}

if (exists $opt{m})
{
	if (exists $validModes{$opt{m}})
	{
		$alignMode = $opt{m};
	}
	else
	{
		usage("-m $opt{m} is not a valid alignment mode.");
	}
}

# interactive options
if (exists $opt{i})
{
	if ($opt{i} =~ /([+-])([ACGT]*)/)	# in the format +A or -TAG
	{
		$type = $1;
		$indel = $2;
		$interactiveMode = 1;
	}
	else
	{
		usage("Couldn't parse indel format, should be in the style of +A or -TAG.");
	}
}
if (exists $opt{c})
{
	$chr = $opt{c};
	$interactiveMode = 1;
}
if (exists $opt{p})
{
	$pos = $opt{p};
	$interactiveMode = 1;
}

if ($interactiveMode == 1)
{
	unless ((exists $opt{i}) and (exists $opt{c}) and (exists $opt{p}))
	{
		usage("-i, -c and -p must all be supplied for interactive mode.");
	}
}

# VCF options
if (exists $opt{v})
{
	if ($interactiveMode == 1)
	{
		usage("Please specifiy either a single indel on the command line or a vcf file to process, not both.");
	}
	elsif (-e $opt{v})
	{
		$vcfPath = $opt{v};
	}
	else
	{
		usage("-v $opt{v} does not seem to be a valid file.");
	}
}
if (exists $opt{o})
{
	$doSort = 0;
}




###########################################################
# process indel(s)
###########################################################
my $candidatePos;
my $chunkSize = 10 ** $chunkPower;
my %referenceHash = ("chunkSize" => $chunkSize);		# chunk size is how many bases will be pulled from the fasta reference at a time
														# referenceHash will also contain the bases in the form $hash{chromosome}{chunk}{position}
my %unsortedOutputChunks;		# $unsortedOutputChunks{$chr}{$chunk}{$line}

# counts for vcf stat output
my $snpCount = 0;
my $insCount = 0;
my $delCount = 0;
my $commaCount = 0;
my $noMatchCount = 0;

my %ambigCount;
my %leftCount;
my %rightCount;

my $ref;
my $var;
my $dbsnp;

if ($interactiveMode == 1)
{
	$candidatePos = findCandidatePositions($type, $indel, $chr, $pos, \%referenceHash, \%fastaHandles);
	
	if (exists $candidatePos->{"deletionNotPossible"})
	{
		die "Exiting.\n";
	}
	
	# format and output
	my @sortedPositions = sort { $a <=> $b } keys %$candidatePos;
	my $outputPos;
	my $outputIndel;
	my $outputAnchor;

	my $outputLine;
	
	if ($alignMode eq "left")
	{
		$outputPos = $sortedPositions[0];
		$outputIndel = $candidatePos->{$outputPos};
		
		$outputAnchor = getAnchorBase($chr, $outputPos, $indel, $anchorMode, \%referenceHash, \%fastaHandles);
		$outputPos = adjustAnchorPos($outputPos, $anchorPos, $type, $outputIndel);
		$outputLine = printIndel($chr, $outputPos, ".", $outputAnchor, $anchorMode, $type, $outputIndel, "\n", "\t");
		print $outputLine;
	}
	elsif ($alignMode eq "right")
	{
		$outputPos = $sortedPositions[$#sortedPositions];
		$outputIndel = $candidatePos->{$outputPos};

		$outputAnchor = getAnchorBase($chr, $outputPos, $indel, $anchorMode, \%referenceHash, \%fastaHandles);
		$outputPos = adjustAnchorPos($outputPos, $anchorPos, $type, $outputIndel);
		$outputLine = printIndel($chr, $outputPos, ".", $outputAnchor, $anchorMode, $type, $outputIndel, "\n", "\t");
		print $outputLine;
	}
	elsif ($alignMode eq "all")
	{
		for $outputPos (@sortedPositions)
		{
			$outputIndel = $candidatePos->{$outputPos};
			
			$outputAnchor = getAnchorBase($chr, $outputPos, $indel, $anchorMode, \%referenceHash, \%fastaHandles);
			$outputPos = adjustAnchorPos($outputPos, $anchorPos, $type, $outputIndel);
			$outputLine = printIndel($chr, $outputPos, ".", $outputAnchor, $anchorMode, $type, $outputIndel, "\n", "\t");			
			print $outputLine;
		}
	}
	elsif ($alignMode eq "full")
	{
		my $refPad = 3;
		doFullOutput($chr, $pos, $type, $candidatePos, $anchorMode, $anchorPos, \%referenceHash, \%fastaHandles);#
	}

}
else		# process vcf file
{
	open(VCFFILE, $vcfPath) or die "Couldn't open $vcfPath for reading.\n";

	my $line;
	my @fields;

	while ($line = <VCFFILE>)
	{
		# read line from file
		if ($line =~ /^#/)
		{
			print $line;
		}
		else
		{
			chomp $line;
			@fields = split(/\t/, $line);
	
			$chr = $fields[0];
			$pos = $fields[1];
			$dbsnp = $fields[2];
			$ref = $fields[3];
			$var = $fields[4];

			# tests to make sure we've got a workable indel
			if (length($ref) == length($var))
			{
				$snpCount++;
				if ($doSort == 0)
				{
					print $line . "\n";
				}
				else
				{
					doSortedOutput($chr, $pos, $line, \%unsortedOutputChunks);
				}
			}
			elsif (($ref =~ /,/) or ($var =~ /,/))
			{
				$commaCount++;
				if ($doSort == 0)
				{
					print $line . "\n";
				}
				else
				{
					doSortedOutput($chr, $pos, $line, \%unsortedOutputChunks);
				}
			}
			else
			{
				if (length($ref) > length($var))
				{
					$type = "-";	# deletion
					$indel = substr($ref, 1);		# 1 to drop the anchor base
					$delCount++;
				}
				else
				{
					$type = "+";	# insertion
					$indel = substr($var, 1);		# 1 to drop the anchor base
					$insCount++;
				}
				
				# attempt to free mem (assuming sorted, once we're > 1 chunk away we can free it)
				freeMemChunks($chr, $pos, \%referenceHash);
				
				# find ambiguous alignments
				$candidatePos = findCandidatePositions($type, $indel, $chr, $pos, \%referenceHash, \%fastaHandles);

				if (exists $candidatePos->{"deletionNotPossible"})
				{
					if ($doSort == 0)
					{
						print $line . "\n";
					}
					else
					{
						doSortedOutput($chr, $pos, $line, \%unsortedOutputChunks);
					}
					$noMatchCount++;
				}
				else
				{
					my @sortedPositions = sort { $a <=> $b } keys %$candidatePos;
					my $outputPos;
					my $outputIndel;
					my $outputAnchor;
					my $outputLine;
	
					if ((scalar keys %$candidatePos) > 1)
					{
						$ambigCount{$type}++;
						if ($sortedPositions[0] == $pos)
						{
							$leftCount{$type}++;
						}
						elsif ($sortedPositions[$#sortedPositions] == $pos)
						{
							$rightCount{$type}++;
						}
					}
	
					# format appropriately and output
					
					my $restOfLine = $fields[5];
					for (my $i = 6; $i < scalar @fields; $i++)
					{
						$restOfLine .= "\t$fields[$i]";
					}
					if ($alignMode eq "left")
					{
						$outputPos = $sortedPositions[0];
						$outputIndel = $candidatePos->{$outputPos};
	
						
						$outputAnchor = getAnchorBase($chr, $outputPos, $indel, $anchorMode, \%referenceHash, \%fastaHandles);
						$outputPos = adjustAnchorPos($outputPos, $anchorPos, $type, $outputIndel);
						if ($doSort == 0)
						{
							$outputLine = printIndel($chr, $outputPos, $dbsnp, $outputAnchor, $anchorMode, $type, $outputIndel, "$restOfLine\n", "\t");
							print $outputLine;
						}
						else
						{
							$outputLine = printIndel($chr, $outputPos, $dbsnp, $outputAnchor, $anchorMode, $type, $outputIndel, "$restOfLine\n", "\t");
							doSortedOutput($chr, $outputPos, $outputLine, \%unsortedOutputChunks);
						}
					}
					elsif ($alignMode eq "right")
					{
						$outputPos = $sortedPositions[$#sortedPositions];
						$outputIndel = $candidatePos->{$outputPos};
				
						$outputAnchor = getAnchorBase($chr, $outputPos, $indel, $anchorMode, \%referenceHash, \%fastaHandles);
						$outputPos = adjustAnchorPos($outputPos, $anchorPos, $type, $outputIndel);
						if ($doSort == 0)
						{
							$outputLine = printIndel($chr, $outputPos, $dbsnp, $outputAnchor, $anchorMode, $type, $outputIndel, "$restOfLine\n", "\t");
							print $outputLine;
						}
						else
						{
							$outputLine = printIndel($chr, $outputPos, $dbsnp, $outputAnchor, $anchorMode, $type, $outputIndel, "$restOfLine\n", "\t");
							doSortedOutput($chr, $outputPos, $outputLine, \%unsortedOutputChunks);
						}
					}
				}
			}
		}
	}
	close VCFFILE;
	# output stats
	if ($outputVcfStats == 1)
	{
		for $type (qw(+ -))
		{
			unless (exists $ambigCount{$type})
			{
				$ambigCount{$type} = 0;
			}
			unless (exists $leftCount{$type})
			{
				$leftCount{$type} = 0;
			}
			unless (exists $rightCount{$type})
			{
				$rightCount{$type} = 0;
			}
		}

		warn "\n";
		warn "# snps:           $snpCount\n";
		warn "# insertions:     $insCount\n";
		warn "# deletions:      $delCount\n";
		warn "# no match dels:  $noMatchCount\n";
		warn "# comma vars:     $commaCount\n";
		warn "\n";
		warn "(the following is before generalization)\n";
		warn "# ambiguous insertions:     $ambigCount{'+'}\n";
		warn "# left aligned insertions:  $leftCount{'+'}\n";
		warn "# right aligned insertions: $rightCount{'+'}\n";
		warn "# other aligned insertions: " . ($ambigCount{'+'} - $leftCount{'+'} - $rightCount{'+'}) . "\n";
		warn "\n";
		warn "# ambiguous deletions:      $ambigCount{'-'}\n";
		warn "# left aligned deletions:   $leftCount{'-'}\n";
		warn "# right aligned deletions:  $rightCount{'-'}\n";
		warn "# other aligned deletions:  " . ($ambigCount{'-'} - $leftCount{'-'} - $rightCount{'-'}) . "\n";
		warn "\n";
	}
}


# getAnchorBase considers the position of the indel and current anchor mode in order to return the anchor base
# input is the chromosome, position and content of the indel (in order to determine length for right-most anchor modes) as well as the reference and fasta handles
# output is the anchor base which is returned
sub getAnchorBase
{
	my $chr = $_[0];
	my $pos = $_[1];
	my $indel = $_[2];
	my $anchorMode = $_[3];
	my $reference = $_[4];
	my $fastaHandles = $_[5];

	my $outputAnchor;
	
	if ($anchorMode eq "left")
	{
		$outputAnchor = getBase($chr, $pos, $reference, $fastaHandles);
	}
	elsif ($anchorMode eq "right")
	{
		$outputAnchor = getBase($chr, $pos + length($indel) + 1, $reference, $fastaHandles);
	}
	else	# anchorMode eq "off"
	{
		$outputAnchor = "-";
	}
	
	return $outputAnchor;
}


# adjustAnchorPos consideres the type of indel and current anchor mode in order to adjust the position of the indel
# input is the position, current anchor position mode, indel type and indel context (in order to determine indel length for right-most anchors)
# output is the adjusted anchor pos
sub adjustAnchorPos
{
	my $pos = $_[0];
	my $anchorPos = $_[1];
	my $type = $_[2];
	my $indel = $_[3];
	
	if ($anchorPos eq "left-1")
	{
		return $pos;
	}
	elsif ($anchorPos eq "left")
	{
		if ($type eq "+")
		{
			return $pos;
		}
		else	# type is -
		{
			return $pos + 1;
		}
	}
	elsif ($anchorPos eq "right")
	{
		if ($type eq "+")
		{
			return $pos + length($indel) + 1;
		}
		else	# type is -
		{
			return $pos + length($indel);
		}
	}
	else	# mode is right+1
	{
		return $pos + length($indel) + 1;
	}
}

# printIndel takes an indel and some meta data and formats the output line
# input is the chromosome and position of the indel, any string for the id field (e.g. dbsnp), the anchor base and where to place it, the type of indel and its change
# 	and finally the rest of the output line and the delimitor to use to separate the fields
# output is the formatted string which is returned
sub printIndel
{
	my $chr = $_[0];
	my $pos = $_[1];
	my $id = $_[2];
	my $anchorBase = $_[3];
	my $anchorMode = $_[4];
	my $type = $_[5];
	my $indel = $_[6];
	my $line = $_[7];
	my $delim = $_[8];
	
	my $insertionRef;		
	my $deletionRef;
	
	if ($anchorMode eq "left")
	{
		$insertionRef = $anchorBase;
		$deletionRef = $anchorBase . $indel;
	}
	elsif ($anchorMode eq "right")
	{
		$insertionRef = $anchorBase;
		$deletionRef = $indel . $anchorBase;			
	}
	else	# anchor mode is off
	{
		$insertionRef = $anchorBase;
		$deletionRef = $indel;
	}

	if ($type eq "+")	# insertion
	{
		return "$chr$delim$pos$delim$id$delim$insertionRef$delim$deletionRef$delim$line";
	}
	else	# deletion
	{
		return "$chr$delim$pos$delim$id$delim$deletionRef$delim$insertionRef$delim$line";
	}
}

# doFullOutput takes the indel of interest and a list of candidate positions and prints 'full mode' output to stdout
# input is ($chr, $pos, $type, $candidatePos, $anchorMode, $anchorPos, \%referenceHash, \%fastaHandles)
# outputs to stdout and doesn't return anything
sub doFullOutput
{
	my $chr = $_[0];
	my $pos = $_[1];
	my $type = $_[2];
	my $candidatePos = $_[3];
	my $anchorMode = $_[4];
	my $anchorPos = $_[5];
	my $referenceHash = $_[6];
	my $fastaHandles = $_[7];


	my $outputPos;
	my $outputIndel;
	my $outputAnchor;
	my $fullIndelDisplay;
	my $outputLine;

	my @sortedPositions = sort { $a <=> $b } keys %$candidatePos;
	my $indelLength = length($candidatePos->{$sortedPositions[0]});

	my $gapChars = "";
	for (my $i = 0; $i < $indelLength; $i++)
	{
		if ($type eq "+")
		{
			$gapChars .= "*";
		}
		elsif ($type eq "-")
		{
			$gapChars .= "-";
		}
	}


	my $leadingPad = 15;
	my $refPad = 3;
	my $tabPad = "   ";

	my $leadingStart = $sortedPositions[0] - $refPad + 1;
	my $trailingEnd = $sortedPositions[$#sortedPositions] + $refPad;
	if ($type eq "-")
	{
		$trailingEnd += $indelLength;
	}

	my $fullPad = spaces(length($chr) + $leadingPad + (length($tabPad) * 5) + $indelLength + 2);		# 2 is for the dbSNP "." and the anchor base
	unless ($anchorMode eq "off")
	{
		$fullPad .= " ";
	}

	print "\n";
	print "$fullPad|$leadingStart\n";


	if ($type eq "-")
	{
		print "(reference sequence)";
		print spaces(length($fullPad) - length('(reference sequence)'));
		print getRange($chr, $leadingStart, $trailingEnd, \%referenceHash, \%fastaHandles) . "\n";
	}
	else
	{
		print "(aligned consensus)";
		print spaces(length($fullPad) - length('(aligned consensus)'));
		print getRange($chr, $leadingStart, $pos, \%referenceHash, \%fastaHandles);
		print $indel;
		print getRange($chr, $pos + 1, $trailingEnd, \%referenceHash, \%fastaHandles) . "\n";
	}

	for $outputPos (@sortedPositions)
	{
		$outputIndel = $candidatePos->{$outputPos};
		

		$fullIndelDisplay = "";
		$fullIndelDisplay .= spaces($leadingPad - length($outputPos));

		if ($type eq "-")
		{
			$fullIndelDisplay .= getRange($chr, $leadingStart, $outputPos, \%referenceHash, \%fastaHandles);
			$fullIndelDisplay .= $gapChars;
			$fullIndelDisplay .= getRange($chr, $outputPos + $indelLength + 1, $trailingEnd, \%referenceHash, \%fastaHandles);
		}
		else # type is insertion
		{
			$fullIndelDisplay .= getRange($chr, $leadingStart, $outputPos, \%referenceHash, \%fastaHandles);
			$fullIndelDisplay .= $gapChars;
			$fullIndelDisplay .= getRange($chr, $outputPos + 1, $trailingEnd, \%referenceHash, \%fastaHandles);
		}


		if ($outputPos == $pos)
		{
			if ($type eq "+")
			{
				$fullIndelDisplay .= "$tabPad<-- (original insertion)";
			}
			else
			{
				$fullIndelDisplay .= "$tabPad<-- (original deletion)";
			}
		}
		$outputAnchor = getAnchorBase($chr, $outputPos, $indel, $anchorMode, \%referenceHash, \%fastaHandles);
		$outputPos = adjustAnchorPos($outputPos, $anchorPos, $type, $outputIndel);
		$outputLine = printIndel($chr, $outputPos, ".", $outputAnchor, $anchorMode, $type, $outputIndel, "$fullIndelDisplay\n", $tabPad);
		print $outputLine;
	}

	if ($type eq "-")
	{
		print "(reference sequence)";
		print spaces(length($fullPad) - length('(reference sequence)'));
		print getRange($chr, $leadingStart, $trailingEnd, \%referenceHash, \%fastaHandles) . "\n";
	}
	else
	{
		print "(aligned consensus)";
		print spaces(length($fullPad) - length('(aligned consensus)'));
		print getRange($chr, $leadingStart, $pos, \%referenceHash, \%fastaHandles);
		print $indel;
		print getRange($chr, $pos + 1, $trailingEnd, \%referenceHash, \%fastaHandles) . "\n";
	}

	print spaces(length($chr) + $leadingPad + (length($tabPad) * 5) + $indelLength + 2 + ($trailingEnd - $leadingStart) - length($trailingEnd));
	unless ($anchorMode eq "off")
	{
		print " ";
	}
	if ($type eq "+")
	{
		print spaces($indelLength);
	}
	print "$trailingEnd|\n";
	print "\n";
}

# doSortedOutput takes a line that needs to be outputted and stores it in a hash which is broken up by chunks.  Once sufficiently far from a chunk (assuming sorted input)
# 	the lines in the chunk are sorted and output to stdout.  It's possible (however unlikely) that an indel could be generalized to a position that is far enough away that
# 	its chunk has already been outputted and freed from memory, in which case there is a warning that the output vcf is not sorted
# input is the chromosome, position, output line and handle for the hash of unsorted output chunks
# output is sorted chunks to stdout, doesn't return anything
sub doSortedOutput
{
	my $chr = $_[0];
	my $pos = $_[1];
	my $line = $_[2];
	my $unsortedOutputChunks = $_[3];

	my $chunkSize = 10000;		# could adjust this up if indels are being generalized more than 10kb away from their reported positions (there is a warning for this)
	my $currentChunk = int(($pos - 1) / $chunkSize) * $chunkSize + 1;

	# add line to hash
	$unsortedOutputChunks->{$chr}{$currentChunk}{$pos}{$line}++;

	# check if we're on a new chromosome - sort and output old chromosomes
	for my $refChr (keys %$unsortedOutputChunks)
	{
		if ($refChr ne $chr)
		{
			for my $refChunk (sort { $a <=> $b } keys %{ $unsortedOutputChunks->{$refChr} })
			{
				printChunk($refChr, $refChunk, $unsortedOutputChunks);
			}
			delete $unsortedOutputChunks->{$refChr};
		}
	}

	
	# check if we're far enough away from existing chunks - sort and output chunks that are no longer needed
	for my $chunkPos (sort { $a <=> $b } keys %{ $unsortedOutputChunks->{$chr} })
	{
		if ($chunkPos < ($pos - (2 * $chunkSize)))
		{
			printChunk($chr, $chunkPos, $unsortedOutputChunks);
			delete $unsortedOutputChunks->{$chr}{$chunkPos};
		}

	}

	return;
}

# printChunk takes a chunk and a reference to the unsorted chunk hash as input, sorts the chunks by the position and prints them to stdout
sub printChunk
{
	my $chr = $_[0];
	my $chunk = $_[1];
	my $unsortedOutputChunks = $_[2];

	for my $pos (sort { $a <=> $b } keys %{ $unsortedOutputChunks->{$chr}{$chunk} })
	{
		for my $line (keys %{ $unsortedOutputChunks->{$chr}{$chunk}{$pos} })
		{
			for (my $i = 0; $i < $unsortedOutputChunks->{$chr}{$chunk}{$pos}{$line}; $i++)	# print multiple lines multiple times
			{
				print $line;
			}
		}
	}
	return;
}


# returns $_[0] spaces
sub spaces
{
	my $num = $_[0];
	my $output = "";
	for (my $i = 0; $i < $num; $i++)
	{
		$output .= " ";
	}
	return $output;
}


# getBase returns the base at a specific position in the reference.  It will initialize fasta handles and pull new chunks of reference if necessary
# input is chromosome, position, reference hash and fasta handles
# output is a single base (and the reference hash and fasta handles may be modified)
sub getBase
{
	my $chr = $_[0];
	my $pos = $_[1];
	my $reference = $_[2];
	my $fastaHandles = $_[3];

	my $chunkStart = int(($pos - 1) / $reference->{"chunkSize"}) * $reference->{"chunkSize"} + 1;		# +1 because the first base in the reference is 1, $pos - 1 so that multiples of chunk size resolve to the correct chunk
	my $chunkEnd = $chunkStart + $reference->{"chunkSize"} - 1;
	
	unless (exists $reference->{$chr}{$chunkStart}{$pos})		# if the position isn't in our hash, we need to get a new chunk from the reference
	{
		unless (exists ($fastaHandles->{$chr}))		# create a handle for the chromosome fasta, if it doesn't exist
		{
#			warn "Creating fasta handle for $chr\n";
			$fastaHandles->{$chr} = Bio::DB::Fasta->new("$fastaHandles{path}/$chr.fa");
		}
		
#		warn "Pulling $chr:$chunkStart-$chunkEnd from fasta\n";
		my $newChunk = uc($fastaHandles{$chr}->seq($chr, $chunkStart, $chunkEnd));
		my $i = $chunkStart;
		for my $base (split("", $newChunk))
		{
			$reference->{$chr}{$chunkStart}{$i} = $base;
			$i++;
		}
	}
#	warn "returning $reference->{$chr}{$chunkStart}{$pos}\n";
	return $reference->{$chr}{$chunkStart}{$pos};
}


# getRange returns a string of bases from the reference in the specified range by calling getBase
# input is chromosome, start pos, end pos, reference hash and fasta handles
# output is a string of bases
sub getRange	
{
	my $chr = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	my $reference = $_[3];
	my $fastaHandles = $_[4];

	my $seq = "";

	for (my $p = $start; $p <= $end; $p++)
	{
		$seq .= getBase($chr, $p, $reference, $fastaHandles);
#		warn "Got base: $chr:$p\t$seq\n";
	}

	return $seq;
}

# freeMemChunks tests if the next indel to be processed is on a different chromosome or more than a chunk away from the reference sequences currently in memory
# 	if there exist chunks that we won't need again (assuming the input is sorted) the chunks will be removed from the reference hash
# input is the chromosome and position of the current indel, and the reference hash
# there is no output
sub freeMemChunks
{
	my $chr = $_[0];
	my $pos = $_[1];
	my $reference = $_[2];
	
	# delete chunks from non-current chromosomes
	for my $refChr (keys %$reference)
	{
		if (($refChr ne $chr) and ($refChr ne "chunkSize"))
		{
#			warn "deleting all chunks for $refChr.\n";
			delete $reference->{$refChr};
		}
	}
	
	# delete chunks if they are more than 1.5 chunks away from the current indel
	# 1.5 so that we are at least in the middle of the current chunk before dropping the previous one
	for my $chunkPos (keys %{ $reference->{$chr} })
	{
		if ($chunkPos < ($pos - (1.5 * $reference->{"chunkSize"})))
		{
#			warn "deleting $chr:$chunkPos chunk.\n";
			delete $reference->{$chr}{$chunkPos};
		}
	}
	
	return;
}


# findCandidatePositions compares the indel sequence to the surrounding reference bases and builds a list of analogous indels
# input is the type and sequence of the indel, its chromosome and position, and references to the persistant reference hash and fasta handle hash
# output is a hash of candidate positions and resulting indel sequences
sub findCandidatePositions
{
	my $type = $_[0];
	my $indel = $_[1];
	my $chr = $_[2];
	my $pos = $_[3];
	my $reference = $_[4];
	my $fastaHandles = $_[5];

	my %candidatePos;
	my $candidateIndel;
	my $smallestInsertionSubunit;
	my $sortedIndel = join('',sort(split('',$indel)));
	my $indelLength = length($indel);

	if ($type eq "+")   # won't be able to find the insertion in reference seqeuence
	{
		$candidatePos{$pos} = $indel;
		$smallestInsertionSubunit = findSmallestInsertionSubunit($indel);
	}
	else
	{
		# check if deletion is possible
		
		$candidateIndel = getRange($chr, $pos + 1, $pos + 1 + $indelLength - 1, \%referenceHash, \%fastaHandles);		# pos + 1 because pos is left of the deletion

		if ($indel eq $candidateIndel)
		{
			$candidatePos{$pos} = $indel;
		}
		else
		{
			warn "Deletion of $indel anchored at $chr:$pos is not possible (the reference is $candidateIndel), skipping generalization.\n";
			%candidatePos = ("deletionNotPossible" => 1);
			return \%candidatePos;
		}
	}

	
	# search left
	my $done = 0;
	my $searchPos = $pos;

	my $lastFoundInsertion = $pos;


	while ($done == 0)
	{
		if ($type eq "+")	# insertion
		{
			$searchPos -= length($smallestInsertionSubunit);
			$candidateIndel = getRange($chr, $searchPos + 1, $searchPos + 1 + length($smallestInsertionSubunit) - 1, \%referenceHash, \%fastaHandles);  # check the bases to the right of the potential insertion for a subunit match

			if ($candidateIndel eq $smallestInsertionSubunit)
			{
				$candidatePos{$searchPos} = $indel;
			}
			else
			{
				$done = 1;
			}
		}
		else	# deletion
		{
			$candidateIndel = getRange($chr, $searchPos, $searchPos + $indelLength - 1, \%referenceHash, \%fastaHandles);  # old code: substr($reference, $searchPos, $indelLength);
			if (join('',sort(split('',$candidateIndel))) eq $sortedIndel)
			{
				$candidatePos{$searchPos - 1} = $candidateIndel;	# searchPos - 1 so that we're returning the anchor base
			}
			else
			{
				$done = 1;
			}
			$searchPos--;
		}
	}


	# search right
	$searchPos = $pos;
	if ($type eq "+")
	{
		$searchPos -= length($smallestInsertionSubunit);
	}
	$done = 0;

	$lastFoundInsertion = $pos;

	while ($done == 0)
	{

		if ($type eq "+")
		{
			$searchPos += length($smallestInsertionSubunit);
			$candidateIndel = getRange($chr, $searchPos + 1, $searchPos + 1 + length($smallestInsertionSubunit) - 1, \%referenceHash, \%fastaHandles);  # search the bases to the right of the position base for a subunit match

			if ($candidateIndel eq $smallestInsertionSubunit)
			{
				$candidatePos{$searchPos} = $indel;
				$candidatePos{$searchPos + length($smallestInsertionSubunit)} = $indel;     # can insert on either side of a repeated motif

			}
			else
			{
				$done = 1;
			}
		}
		else    # deletion
		{
			$searchPos++;
			$candidateIndel = getRange($chr, $searchPos, $searchPos + $indelLength - 1, \%referenceHash, \%fastaHandles);  # old code: substr($reference, $searchPos, $indelLength);
			if (join('',sort(split('',$candidateIndel))) eq $sortedIndel)
			{
				$candidatePos{$searchPos - 1} = $candidateIndel; # searchPos - 1 so that we're returning the anchor base
			}
			else
			{
				$done = 1;
			}
		}
	}
	
	
	return \%candidatePos;
	
}

# findSmallestInsertionSubunit returns the simplist subunit inside an insertion (e.g. for an insertion of +AAA, A is the simplist subunit)
# input is the insertion string
# output is the subunit string
sub findSmallestInsertionSubunit
{
	my $insertion = $_[0];

	my $sim;
	my $simFound = 0;

	my @subunits;

	for (my $i = 1; $i * 2 <= length($insertion); $i++)
	{
		if ((length($insertion) % $i) == 0)
		{
			$sim = substr($insertion, 0, $i);
			@subunits = ();

			for (my $j = length($sim); $j < length($insertion); $j += length($sim))
			{
				push(@subunits, substr($insertion, $j, length($sim)));
			}

			$simFound = 1;	  # innocent until guilty
			for (my $j = 0; $j < scalar(@subunits); $j++)
			{
				unless ($sim eq $subunits[$j])
				{
					$simFound = 0;
				}
			}

			if ($simFound == 1)
			{
				return $sim;
			}
		}
	}

	return $insertion;
}


