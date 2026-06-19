#!/usr/bin/perl

use POSIX;


my $scriptName = "study.quality_of_mismatches.for_GitHub.pl";

my %arg;
&parseCommandLine;
#
# required
#
my $prefix = $arg{-prefix};
my $bam = $arg{-bam};
my $output_directory = $arg{-output_directory};
#
# optional
#
my $dump_details = $arg{-dump_details};

#
# fill defaults if necessary
#
# SET TO YES TO GET DETAILS ABOUT READS WITH INDEL MISMATCHES
#   (usually some single-nucleotide mismatches get output as collateral)
#
if (exists($arg{-dump_details}) && ($arg{-dump_details} eq "yes" || $arg{-dump_details} eq "y")) {
    $dump_details = $arg{-dump_details};
} else {
	$dump_details = "no";
}

print STDOUT "\nRunning $scriptName with the following arguments:\n";
print STDOUT "\e[1m-prefix\e[0m\t$prefix\n";
print STDOUT "\e[1m-bam\e[0m\t$bam\n";
print STDOUT "\e[1m-output_directory\e[0m\t$output_directory\n";
print STDOUT "\e[1m-dump_details\e[0m\t$dump_details\n";
#
#################################################################################

my $dumpReadDetails = 0;
if ($dump_details eq "yes" || $dump_details eq "y") {
	$dumpReadDetails = 1;
}


# for CIGAR string analysis
#
my %cigarHash = ("M",1,"I",1,"D",1,"N",1,"S",1,"H",1,"P",1,"X",1,"=",1);
my %matchHash = ("A",1,"C",1,"T",1,"G",1,"^",1);
#
my %cigarCountHash = ();

my $headerTag = "\@";

my @alleleArray = ("A","C","G","T");


my $indelbyqual = "$output_directory/indel.summary.$prefix.txt";
open(INDBYQUAL,">$indelbyqual") || die "can't open site by qual indel output: $indelbyqual\n";
print INDBYQUAL "type\tlength\tavgBQ\tinstances\n";

my $snmbyqual = "$output_directory/single-nucleotide.summary.$prefix.txt";
open(SNMBYQUAL,">$snmbyqual") || die "can't open snm by qual output: $snmbyqual\n";
print SNMBYQUAL "ref\talt\tflanks\tBQ\tinstances\n";

my %snmHash = ();

my $logfile = "$output_directory/logfile.$prefix.txt";
open(LOG,">$logfile") || die "can't open logfile: $logfile\n";

my $ratesfile = "$output_directory/raw_mismatch_rates.summary.$prefix.txt";
open(RATES,">$ratesfile") || die "can't open raw mismatch rates file: $ratesfile\n";

my $ARLfile = "$output_directory/total_adjusted_read_length.$prefix.txt";
open(ARL,">$ARLfile") || die "can't open ARL file: $ARLfile\n";

my $rlfile = "$output_directory/read_length_distribution.$prefix.txt";
open(RL,">$rlfile") || die "can't open RL dist file: $rlfile\n";
print RL "read_length\tcounts\n";
my $isfile = "$output_directory/insert_size_distribution.$prefix.txt";
open(IS,">$isfile") || die "can't open insert size dist file: $isfile\n";
print IS "insert_size\tcounts\n";

my $snmbed = "$output_directory/single-nucleotide.details.$prefix.bed";
open(SNMBED,">$snmbed") || die "can't open snm bed file: $snmbed\n";
print SNMBED "chrom\tstart\tstop\tquality\tfrac_pos\tchange\n";

my $insbed = "$output_directory/insertion.details.$prefix.bed";
open(INSBED,">$insbed") || die "can't open insertions bed file: $insbed\n";
print INSBED "chrom\tstart\tstop\tquality\tchange\n";

my $delbed = "$output_directory/deletion.details.$prefix.bed";
open(DELBED,">$delbed") || die "can't open deletions bed file: $delbed\n";
print DELBED "chrom\tstart\tstop\tquality\tchange\n";


my $totalRL = 0;
my $totalIS = 0;
my $consideredRLReads = 0;
my $consideredISReads = 0;

my $rlisfile = "$output_directory/mean_and_median.read_length_and_insert_size.$prefix.txt";
open(RLIS,">$rlisfile") || die "can't open mean RL insert size file: $rlisfile\n";
print RLIS "mean_RL\tmedian_RL\tmean_IS\tmedian_IS\n";


my %allMatchHash = ();
my $readCount=0;

my $delTag = "del";
my $insTag = "ins";

my $strandFlag = 0x10; # bit 4 defines strand

my %rlHash = ();
my %isHash = ();

my %indelsByInstanceHash = ();

my $totalSNMs = 0;
my $totalInsertions = 0;
my $totalDeletions = 0;

print LOG "samtools view $bam  \n";
open(BAM,"samtools view $bam |") || die "can't open pipe to bam file: $bam\n";
while(<BAM>) {
	my $line = "$_";
	chomp($line);
	my @lineArray = split(/\t/,$line);
	if (substr($lineArray[0],0,1) ne $headerTag) {

		my $id = $lineArray[0];
		my $readName = $id;
		$id =~ s/\//\_/g;

		my $toChrom = $lineArray[2];
		my $mappedTo = $lineArray[3];
		
		$readCount++;
		if ($readCount%1000000 == 0) {
		#if ($readCount%10000 == 0) {
			print STDOUT "analyzed $readCount $reads\n";
			#last;
		}
		
		
		my $read_proceed = 1; # archaic
		
		if ($read_proceed) {
		
			my $cigarStr = $lineArray[5];
		
			my $bitwiseflag = $lineArray[1];
		
			my $strand = "+";
			if ($bitwiseflag & $strandFlag) {
				$strand = "-";
			} 
		
			my $mdStr = "";
			#
			# analyze the CIGAR string
			#
			my $integerStr = "";
			my $softclipSkip = 0;
			my $softclipRemove = 0;
			my $hardclipSkip = 0;
			my $hardclipRemove = 0;
			my %cigarPosHash = ();
			my %matchPosHash = ();
			my $readPos = 0;
			my $firstS = 1;
			my $firstH = 1;
			my %insertHash = ();
			my %insertHashToKeep = ();
			my $indelType = "NA";
		
			my $adjustedReadLength = 0;
			my $mappedAlleles = 0;
			my $mismatches = 0;
			my $softClips = 0;
			my $hardClips = 0;
			my $insertions = 0;
			my $deletions = 0;
			my %combinedChangesHash = ();
			my %posCoveredThisRead = ();
		
			my $mapQ = $lineArray[4];

			$rlHash{length($lineArray[9])}++;
			if ($lineArray[8] > 0) {			
				$isHash{$lineArray[8]}++;
			}

			$totalRL+= length($lineArray[9]);
			$consideredRLReads++;
			if ($lineArray[8] > 0) {			
				$totalIS+=$lineArray[8];
				$consideredISReads++;
			}		
		
	        for (my $i = 0; $i < length($cigarStr); $i++) {
	            my $character =  substr($cigarStr,$i,1);
				$cigarCountHash{$character}++;
	            if (!exists($cigarHash{$character})) {
	                $integerStr = "$integerStr$character";
	            } else {
					if ($character eq "S" && $firstS) {
						$softclipSkip = $integerStr;
						#print LOG "soft $softclipSkip\n";
						$softClips+=$integerStr;
					} elsif ($character eq "S" && !$firstS) {
						# need to remove this section of soft-clipping from read.
						$softclipRemove = $integerStr;
						$softClips+=$integerStr;
					} elsif ($character eq "H" && $firstH) { # is any of this necessary? 
						$hardclipSkip = $integerStr;    #  not really. adding the "H" in next elsif if adequate
						$hardClips+=$integerStr;
					} elsif ($character eq "H" && !$firstH) {
						$hardclipRemove = $integerStr;
						$hardClips+=$integerStr;
					} elsif (exists($cigarHash{$character}) && $character ne "S" && $character ne "H") {
						# this just accumulates a record of the cigar entries
						if ($character eq "I") {
							$insertHash{($readPos+1)} = $integerStr;
							$insertHashToKeep{($readPos+1)} = $integerStr;
							$cigarPosHash{($readPos+1)}{$character} = $integerStr;
							$readPos+=$integerStr;
							$insertions+=$integerStr;
						} elsif ($character eq "D") {  
							$deletions+=$integerStr;
						} else {
							$readPos+=$integerStr;
							$cigarPosHash{$readPos}{$character} = $integerStr;
							$mappedAlleles+=$integerStr;
						}
						$firstS = 0;
						$firstH = 0;
					}
	                $integerStr = "";
	            }
	        }
		
			#
			# matching string analysis
			#
			$integerStr = "";
			$readPos = 0;
			my $foundMD=0;
			my $mdStr = "";
			for (my $i = 11; $i<@lineArray; $i++) {
				if (substr($lineArray[$i],0,2) eq "MD") { 
					$foundMD=1;
					my $afterCaret = 0;
					my $delStr = "";
					my $lastIntegerStr = 0;
					my $lastReadPos = 0;
					(my $tag, my $z, my $matchStr) = split(/\:/,$lineArray[$i]);
					$mdStr = $matchStr;
				
		            for (my $i = 0; $i < length($matchStr); $i++) {
		                my $character =  substr($matchStr,$i,1);
						$allMatchHash{$character}++;
		                if (!exists($matchHash{$character})) {
		                    $integerStr = "$integerStr$character";
							if ($afterCaret) {
								$matchPosHash{$readPos}{$delStr} = $lastIntegerStr;
								$delStr = "";
							}
							$afterCaret = 0;
		                } else {

							if ($character ne "^" && !$afterCaret) {
								$readPos+=$integerStr+1;
								#
								# need to make a correction for each instertion found
								#
								for (my $k=0; $k<=$readPos; $k++) {
									if (exists($insertHash{$k})){
										$readPos+=$insertHash{$k};
										delete $insertHash{$k}; # remove after it's been used once
									}
								}
								$matchPosHash{$readPos}{$character} = $integerStr;

							} else {
							
								# an example of an MD string: 5T60^A28^T14
								if ($afterCaret){
									$delStr = "$delStr$character";
								} else {
									$delStr = $delTag;
								}
								$afterCaret = 1;
								$readPos+=$integerStr;
								$lastIntegerStr = $integerStr;
							}
							$integerStr = "";
						}
						$lastReadPos = $readPos;
					}
					last;
				}
			}

			# adjust read and qual if soft-clipped
			#
			my $scread = substr($lineArray[9],$softclipSkip,length($lineArray[9])-$softclipSkip-$softclipRemove);
			my $scqualStr = substr($lineArray[10],$softclipSkip,length($lineArray[10])-$softclipSkip-$softclipRemove);

			$adjustedReadLength = length($lineArray[9]) - $softClips;
			$totalAdjustedReadLength+=$adjustedReadLength; # this is the denominator for all mismatch calculations

	
			#
			# THIS WILL DUMP INDEL MISMATCH READS
			#
			if ($dumpReadDetails && ($cigarStr =~ /I/ ||  $cigarStr =~ /D/)){ 
			#if ($dumpReadDetails){  # you can go this route if you want ALL mismatch details
				print LOG "\n>>>>>>>>>>>>>>>>>>>>>>>\n";
				print LOG "locale:\t$toChrom\:".($mappedTo-1)."\-".($mappedTo+length($scread))."\n";
				print LOG "readID:\t$id\n";
				print LOG " CIGAR:\t$cigarStr\n";
				print LOG " MATCH:\t$mdStr\n";
				print LOG "  MapQ:\t$mapQ\n\n";
				print LOG "  read:\t$scread\n";
				print LOG " quals:\t$scqualStr\n";
				print LOG " index:\t";
				my $index = 1;
				#
				# this just adds a ruler for easier checking of correlations
				#
				for (my $k=1;$k<=length($scread); $k++) {
					if ($k%10 == 0) {
						print LOG "|";
						$index=1;
					} else {
						print LOG "$index";
						$index++;
					}

				}
				print LOG "\n";
				#
				# COMBINE EVERYTHING INTO ONE
				#
				print LOG "\nevent_type\tpos\tref\talt\tflanks\tBQ\t(finalQ\tleftOrRight)\n";
			}
		
			foreach my $pos (sort {$a<=>$b} keys(%insertHashToKeep)) {
				my $insertStr = substr($scread,($pos - 1),$insertHashToKeep{$pos});
				$combinedChangesHash{$pos} = "ins$insertStr";
			}
			foreach my $pos (sort {$a<=>$b} keys(%matchPosHash)) {
				foreach my $char (sort keys(%{$matchPosHash{$pos}})) {
					if (exists($combinedChangesHash{$pos}) && $combinedChangesHash{$pos} =~ /ins/) {
						my $insStr =$combinedChangesHash{$pos};
						$insStr =~ s/ins//;
						if ($dumpReadDetails) {
							print LOG "OVERLAPPING INDELs: $pos\t$combinedChangesHash{$pos}\t$char\t$matchPosHash{$pos}{$char}\n";
						}
						$combinedChangesHash{$pos+length($insStr)} = $char;
					} else {
						$combinedChangesHash{$pos} = $char;
					}
				}
			}
			my %eventHash = ();
			foreach my $pos (sort {$a<=>$b} keys(%combinedChangesHash)) {
				
				my $event_proceed = 1; # archaic
	
				if ($event_proceed) {

					my $ref = $combinedChangesHash{$pos};
					my $alt = substr($scread,($pos-1),1);

					my $snmQ = ord(substr($scqualStr,($pos-1),1))-33;

					my $left_flank = substr($scread,($pos - 2),1);
					if (($pos - 2) < 0) {
						$left_flank = "X";
					}

					my $right_flank = substr($scread,($pos),1);
					if ($pos > (length($scread)-1)) {
						$right_flank = "X";
					}
					
					my $flanks = "$left_flank\_$right_flank";
					
					if ($combinedChangesHash{$pos} !~ /$delTag/ && $combinedChangesHash{$pos} !~ /$insTag/){
						
						$snmHash{$ref}{$alt}{$flanks}{$snmQ}++;		
						if ($dumpReadDetails && ($cigarStr =~ /I/ ||  $cigarStr =~ /D/)){ # just dump reads with indels					
							print LOG "mismatch\t$pos\t$ref\t$alt\t$flanks\t$snmQ\n";
						}
						$totalSNMs++;
						
						
						my $read_length = length($scread);
						my $string_pos = $pos;
						if ($strand eq "-") {
							$string_pos = $read_length - ($string_pos - 1);
						}
						
						my $string_pos_frac = "NA";
						if ($read_length > 0) {
							$string_pos_frac = sprintf("%6.4f",$string_pos/$read_length);
						}
						
						my $chrom_pos = $mappedTo+$softclipSkip+$pos;
						print SNMBED "$toChrom\t".($chrom_pos-1)."\t$chrom_pos\t$snmQ\t$string_pos_frac\t$ref>$alt\n";
						
					} else {
						
						my $change = $combinedChangesHash{$pos};

						$changeStr = $change;
						$changeStr =~ s/$insTag//;
						$changeStr =~ s/$delTag//;
		
						my $changeLength = length($changeStr);
						my $InsOrDel = substr($change,0,3);

					    my $character =  substr($scread,$pos,1);
			            my $qcharacter =  substr($scqualStr,$pos,1);
			
						my $intQ = ord($qcharacter)-33;
						my $averageIntQ = "NA";
						if ($InsOrDel eq "del"){
							my $qDel1 = ord(substr($scqualStr,$pos-1,1))-33;
							my $qDel2 = ord(substr($scqualStr,$pos,1))-33;
							$averageIntQ = &round_number(($qDel1 + $qDel2)/2);
							$changeQualStr = "$qDel1,$qDel2";
							$ref = $changeStr;
							$alt = "_";
							if ($dumpReadDetails) {
								print LOG "deletion\t$pos\t$ref\t$alt\t$flanks\t$changeQualStr\t$averageIntQ\n";
							}
							$totalDeletions++;
							
							my $chrom_pos = $mappedTo+$softclipSkip+$pos;
							print DELBED "$toChrom\t".($chrom_pos-1)."\t$chrom_pos\t$averageIntQ\t$ref>$alt\n";
							
						} else {
							my $localReadPos = $pos;
							$changeQualStr = "";
							$origchangeQualStr = "";
							my $origqSum = 0;
							for (my $i = $pos; $i<($pos+$changeLength); $i++) {
								my $qVal = ord(substr($scqualStr,($i-1),1))-33;
								$origqSum += $qVal;
								$origchangeQualStr = "$origchangeQualStr$qVal,";
							}
							chop($origchangeQualStr);
							$changeQualStr = "$origchangeQualStr"."_or_";
							$alt = $changeStr;
							$ref = "_";
							
							#
							# this bit checks for differences in right or left justifications
							#
							my $rightJustified = substr($scread,$pos,$changeLength);
							my $k = $pos;
							until ($rightJustified ne $changeStr) {
								$rightJustified = substr($scread,$k++,$changeLength);
								$localReadPos+=1;
							}
							## back it up one
							$k--;
					
							my $qSum = 0;
							my $qSumStr = "TAKING RIGHT-JUSTIFIED QSUM";
							for (my $qi=($k); $qi<($k+$changeLength); $qi++) {
								my $qVal = ord(substr($scqualStr,($qi-1),1))-33;
								$qSum+=$qVal;
								$changeQualStr = "$changeQualStr$qVal,";
								$right_flank = substr($scread,($qi+length($changeStr)),1);
							}
							if ($origqSum < $qSum) {
								$qSum = $origqSum;
								$qSumStr = "TAKING LEFT-JUSTIFIED QSUM\n";
								$right_flank = substr($scread,($pos+length($changeStr)-1),1);
							}
							chop($changeQualStr);
							$flanks = "$left_flank\_$right_flank";
					
							$averageIntQ = "-";
							if ($changeLength > 0) {
								$averageIntQ = &round_number(($qSum)/$changeLength);
							}
							if ($dumpReadDetails) { 				
								print LOG "insertion\t$pos\t$ref\t$alt\t$flanks\t$changeQualStr\t$averageIntQ\t$qSumStr\n";
							}
							$totalInsertions++;	
																			
							my $chrom_pos = $mappedTo+$softclipSkip+$pos;
							print INSBED "$toChrom\t".($chrom_pos-1)."\t".($chrom_pos+length($changeStr))."\t$averageIntQ\t$ref>$alt\n";
						}						
						
						$indelsByInstanceHash{$InsOrDel}{$changeLength}{$averageIntQ}++;
					} 
					
				}
			}
		}
	}

}
close BAM;
close IDRDS;

close SNMBED;
close INSBED;
close DELBED;


my $medianRL = "NA";
my $rlSum = 0;
foreach my $rl (sort {$a<=>$b} keys(%rlHash)) {
	print RL "$rl\t$rlHash{$rl}\n";
	$rlSum+=$rlHash{$rl};
	if ($rlSum >= ($consideredRLReads/2) && $medianRL eq "NA") {
		$medianRL = $rl;
	}
}
my $medianIS = "NA";
my $isSum = 0;
foreach my $is (sort {$a<=>$b} keys(%isHash)) {
	print IS "$is\t$isHash{$is}\n";
	$isSum+=$isHash{$is};
	if ($isSum >= ($consideredISReads/2) && $medianIS eq "NA") {
		$medianIS = $is;
	}
}

close RL;
close IS;

my $meanRL = "NA";
if ($consideredRLReads > 0) {
	$meanRL = sprintf("%6.1f",$totalRL/$consideredRLReads);
}
my $meanIS = "NA";
if ($consideredISReads > 0) {
	$meanIS = sprintf("%6.1f",$totalIS/$consideredISReads);
}

print RLIS "$meanRL\t$medianRL\t$meanIS\t$medianIS\n";
close RLIS;


print LOG "\n----------------FINAL SUMMARY-------------------\n";

print LOG "Total Adjusted Read Length (bp): $totalAdjustedReadLength\n";
my $snmRate = "NA";
my $insRate = "NA";
my $delRate = "NA";
if ($totalAdjustedReadLength > 0) {
	$snmRate = sprintf("%8.6f",1000*$totalSNMs/$totalAdjustedReadLength);
	$insRate = sprintf("%8.6f",1000*$totalInsertions/$totalAdjustedReadLength);
	$delRate = sprintf("%8.6f",1000*$totalDeletions/$totalAdjustedReadLength);
}
print LOG "mismatch type (raw)\ttotal\trate per kbp\n";
print LOG "single-nucleotide:\t$totalSNMs\t$snmRate\n";
print LOG "insertion:\t$totalInsertions\t$insRate\n";
print LOG "deletion:\t$totalDeletions\t$delRate\n";

print RATES "mismatch type (raw)\ttotal\trate per kbp\n";
print RATES "single-nucleotide:\t$totalSNMs\t$snmRate\n";
print RATES "insertion:\t$totalInsertions\t$insRate\n";
print RATES "deletion:\t$totalDeletions\t$delRate\n";


print ARL "$totalAdjustedReadLength\n";
close ARL;


foreach my $type (keys(%indelsByInstanceHash)) {
	foreach my $length (sort {$a<=>$b} keys(%{$indelsByInstanceHash{$type}})) {
		foreach my $intQ (sort {$a<=>$b} keys(%{$indelsByInstanceHash{$type}{$length}})) {
			print INDBYQUAL "$type\t$length\t$intQ\t$indelsByInstanceHash{$type}{$length}{$intQ}\n";
		}		
	}
}

foreach my $ref (@alleleArray) {
	foreach my $alt (@alleleArray,"N") {
		foreach my $flanks (sort {$a<=>$b} keys(%{$snmHash{$ref}{$alt}})) {
			foreach my $snmQ (sort {$a<=>$b} keys(%{$snmHash{$ref}{$alt}{$flanks}})) {
				print SNMBYQUAL "$ref\t$alt\t$flanks\t$snmQ\t$snmHash{$ref}{$alt}{$flanks}{$snmQ}\n";
			}		
		}
	}
}
close SNMBYQUAL;

close LOG;
close RATES;

sub parseCommandLine
{
	my $useage = "
    You appear to be missing some required arguments in your submission.
    
    \t-prefix\t\e[1m[REQUIRED]\e[0m unique identifier for use in naming files and directories
    
    \t-bam\t\e[1m[REQUIRED]\e[0m full path to indexed, sorted bam file
        
    \t-output_directory\t\e[1m[REQUIRED]\e[0m where you would like all the output files and logs located
    
    \t-complexity_regions\t[OPTIONAL] NOTE: ONLY FOR USE WITH DATA MAPPED TO GRCh38! options: 'all','low','no-low'; default = 'all'
        
    \t-dump_details\t[OPTIONAL] 'yes' or 'y' if you'd like read details; default = 'no'
    
    ";
    
    for (my $i = 0; $i <= $#ARGV; $i++)
    {
        if ($ARGV[$i] =~ /^-/)
        {
            $arg{$ARGV[$i]} = $ARGV[$i+1];
        }
    }
	die($useage) if (!($arg{-prefix}));
	die($useage) if (!($arg{-bam}));
	die($useage) if (!($arg{-output_directory}));
}



sub round_number {
    my ($num) = @_;
    if ($num >= 0) {
        return floor($num + 0.5);
    } else {
        return ceil($num - 0.5);
    }
}
