#!/usr/bin/perl/ -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use List::Util qw ( sum );
use Storable;

######################################################################################################################################################
#
#	Description
#		This is a perl script to discover sequence motifs in a set of fasta alignments in a folder using dreme;
#
#	Input
#		--inDir=					a path contains all fastaSeq;
#		--refFastaPath=				the path pf the reference fasta;
#		--maxk=						maxk in dreme; default = 10;
#		--mink=						mink in dreme; default = 6;
#		--minE=						minE in dreme; default = 0.05;
#		--halfHalfSearch=			"yes" or "no"; search for motifs by comparing the 1st and 2nd halves; default = no
#		--fullShuffleSearch=		"yes" or "no"; search for motifs by comparing the alignment by shuffling itself; default = no, but will be changed to yes if none of the other modes is yes;
#		--posSpecificSearch=		"\d+,\d+:\d+,\d+" or no; position specific search, e.g. "160,200:0,40"; #---positive seq rng base 160-200 on left, negetive seq rng 0-40 on right; default = no
#		--randGenomicSearch=		"yes" or "no"; search for motifs by comparing the alignment to random genomic fragments; default = no
#		--fileFilterTag=			will use only files containning this tag in the inDir; default = "";
#		--maxSeqNum= 				maximum number of sequence to be used, in order to strict running time; default = 999999;
#		--preDefinedDremeXMLPath=	only scan for occurence of these predefined motif and will not search for new motifs usinf dreme; use "no" to skip; default = no;
#		--targetRngStart=			target range start for search the occurence of motif in scanForMotifOccurence subroutine; default = 0, i.e. off
#		--targetRngEnd=				target range end for search the occurence of motif in scanForMotifOccurence subroutine; default = 9999999, i.e. off
#		--outDir=					directory for output;
#
#	Output
#
#	Usage
#		perl dremeMotifSearcher_v0.1.pl --inDir=/Volumes/B_MPro2TB/NGS/results/EHI_polyA_pairEndOnly/finalSAM/trck5_up10.6_down10.6NR.no_Ext.10_Len18.999_NM.2_NH.50/test/fastaSeq/ --refFastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/EHI_v13.fa
#		
#
#	Assumption
#		1. In the gff file, for each feature, "gene" must be appear at first, mRNA/tRNA/ncRNA etc must be appearsed as the second, and exon must be the last, which is the default setting;
#		2. The sam file is sorted;
#		3. The sam file must be edited using SAMMultiHitAndUnalignLabeller.pl, so that the NH:i attribute (number of hit) is added as the last attribute;
#		4. All reads are unspiced, i.e. with the cigar string xxM, while xx is the length
#
#	History:
#		
#		#---built based on polyATailDiscoverer

#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#
#----------Read parameters ----------#
use vars qw ($inDir $refFastaPath $maxk $mink $minE $halfHalfSearch $fullShuffleSearch $posSpecificSearch $randGenomicSearch $fileFilterTag $maxSeqNum $preDefinedDremeXMLPath $targetRngStart $targetRngEnd $outDir);
my ($inDir, $refFastaPath, $maxk, $mink, $minE, $halfHalfSearch, $fullShuffleSearch, $posSpecificSearch, $randGenomicSearch, $fileFilterTag, $maxSeqNum, $preDefinedDremeXMLPath, $targetRngStart, $targetRngEnd, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

#----------Read the fasta file----------#
my ($cntgSeqSSHsh_ref, $cntgLenHsh_ref) = readAndRevComContig($refFastaPath);

my ($seqForMEMEHsh_ref, $alignmentLength) = readInDirForFasta($inDir);

MEMEPolyASeq($seqForMEMEHsh_ref, $cntgSeqSSHsh_ref, $maxk, $mink, $minE, $fullShuffleSearch, $randGenomicSearch, $halfHalfSearch, $posSpecificSearch, $preDefinedDremeXMLPath, $alignmentLength);

exit;
#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	$outDir = "./dremeMotifSearcher/";
	$maxk = 10;
	$mink = 6;
	$minE = 0.05;
	$fullShuffleSearch = "no";
	$halfHalfSearch = "no";
	$posSpecificSearch = "no";
	$randGenomicSearch = "no";
	$fileFilterTag = "";
	$maxSeqNum = 9999999;
	$preDefinedDremeXMLPath = "no";
	$targetRngStart = 0;
	$targetRngEnd = 9999999;
	
	foreach my $param (@ARGV) {
		if ($param =~ m/--inDir=/) {$inDir = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--refFastaPath=/) {$refFastaPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--maxk=/) {$maxk = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--mink=/) {$mink = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--minE=/) {$minE = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--halfHalfSearch=/) {$halfHalfSearch = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--fullShuffleSearch=/) {$fullShuffleSearch = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--posSpecificSearch=/) {$posSpecificSearch = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--randGenomicSearch=/) {$randGenomicSearch = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--fileFilterTag=/) {$fileFilterTag = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--maxSeqNum=/) {$maxSeqNum = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--preDefinedDremeXMLPath=/) {$preDefinedDremeXMLPath = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--targetRngStart=/) {$targetRngStart = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--targetRngEnd=/) {$targetRngEnd = substr ($param, index ($param, "=")+1);} 
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);} 
	}

	if (($fullShuffleSearch eq "no") and ($halfHalfSearch eq "no") and ($posSpecificSearch eq "no") and ($randGenomicSearch eq "no")) {
		print "No search mode was specified. fullShuffleSearch is activated as the default search mode.\n";
		$fullShuffleSearch = "yes";
	}
	
	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	
	system "mkdir -p -m 777 $outDir/";
	system "mkdir -p -m 777 $outDir/meme/";

	open (CMDLOG, ">$outDir/CMD.log.txt");
	print CMDLOG join "", ((join "\t", ("inDir", "refFastaPath", "maxk", "mink", "minE", "halfHalfSearch", "fullShuffleSearch", "posSpecificSearch", "randGenomicSearch", "fileFilterTag", "maxSeqNum", "preDefinedDremeXMLPath", "targetRngStart", "targetRngEnd", "outDir")), "\n");
	print CMDLOG join "", ((join "\t", ($inDir, $refFastaPath, $maxk, $mink, $minE, $halfHalfSearch, $fullShuffleSearch, $posSpecificSearch, $randGenomicSearch, $fileFilterTag, $maxSeqNum, $preDefinedDremeXMLPath, $targetRngStart, $targetRngEnd, $outDir)), "\n");
	close CMDLOG;

	return ($inDir, $refFastaPath, $maxk, $mink, $minE, $halfHalfSearch, $fullShuffleSearch, $posSpecificSearch, $randGenomicSearch, $fileFilterTag, $maxSeqNum, $preDefinedDremeXMLPath, $targetRngStart, $targetRngEnd, $outDir);
}
########################################################################## readAndRevComContig
sub readAndRevComContig {

	my $refFastaPath = $_[0];
	
	my %cntgLenHsh;
	my $cntgSeqHsh_ref = readMultiFasta($refFastaPath);
	my %cntgSeqHsh = %$cntgSeqHsh_ref;
	my %cntgSeqSSHsh; #---strndSpecific

	print "Generating strnd specifc refFasta.\n";	
	foreach my $cntg (keys %cntgSeqHsh) {
		$cntgLenHsh{$cntg} = length ($cntgSeqHsh{$cntg});
		my $seqPlus = $cntgSeqHsh{$cntg};
		my $seqMinus = reverse $seqPlus;
		$seqMinus =~ tr/ACGTacgt/TGCAtgca/;
		${$cntgSeqSSHsh{$cntg}}{"+"} = $seqPlus;
		${$cntgSeqSSHsh{$cntg}}{"-"} = $seqMinus;
		delete $cntgSeqHsh{$cntg};
	}

	return (\%cntgSeqSSHsh, \%cntgLenHsh);
}
########################################################################## readMultiFasta
sub readMultiFasta {

	my $refFastaPath = $_[0];
	my ($seq, $seqName, %fastaHsh);
	my $i = 0;
	print "Reading refFasta into a hash.\n";
	open (INFILE, $refFastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$fastaHsh{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq = $seq.$nextLine;
			$fastaHsh{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}
	close INFILE;
	
	return (\%fastaHsh);
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## GNUPlotXYScatterWithLines
sub GNUPlotXYScatterWithLines {

	my %XYHsh = %{$_[0]};
	my $plotFilePath = $_[1];
	my $plotDataPath = $_[2];
	my $xlable = $_[3];
	my $ylable = $_[4];
	my $xscale = $_[5];
	my $yscale = $_[6];
	my $title = $_[7];
	
	my $GNULogXCmd = "";
	$GNULogXCmd = "set logscale x" if ($xscale eq "log");
	my $GNULogYCmd = "";
	$GNULogYCmd = "set logscale y" if ($yscale eq "log");

	$plotFilePath .= ".pdf" if ($plotFilePath !~ m/\.pdf$/);

	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];

	print "Running GNUPlotXYScatterWithLines for $fileName.\n";
	
	#---creat a tmp file
	open (TMPFILE, ">$plotDataPath");
	for my $x (sort {$a <=> $b} keys %XYHsh) {
		print TMPFILE $x."\t".$XYHsh{$x}."\n";
	}
	close TMPFILE;
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	set output "| ps2pdf - $plotFilePath 2>/dev/null";
	unset logscale x; 
	unset logscale y; 
	$GNULogXCmd;
	$GNULogYCmd;
	set xlabel "$xlable";
	set ylabel "$ylable";
	set title "$title";
	set nokey;
   	plot '$plotDataPath' using 1:2 with lines;
EOPLOT
	close(GNUPLOT);
	#rmtree(['tmp.dat'], 0, 1); #---non-verbose removal of tmp file
}
########################################################################## calculateEntropyAndPlotLogo
sub calculateEntropyAndPlotLogo {
	
	#---depends on subroutine runAndCheckSerialTask
	#calculateEntropyAndPlotLogo($seqAry_ref, $fileTag, $outFolder);
	
	my @seqAry = @{$_[0]};
	my $fileTag = $_[1];
	my $outFolder = $_[2];

	print "Running Weblogo for $fileTag.\n";
	
	my $entropyPath = "$outFolder/$fileTag.entropy.txt";
	my $logoBitPath = "$outFolder/$fileTag.weblogo.bit.pdf";
	my $logoProbPath = "$outFolder/$fileTag.weblogo.prob.pdf";
	my $seqPath = "$outFolder/$fileTag.weblogo.fasta";

	open (WEBLOGOFASTA, ">$seqPath");
	foreach my $seq (@seqAry) {
		print WEBLOGOFASTA ">seq\n";
		print WEBLOGOFASTA "$seq\n";
	}
	close WEBLOGOFASTA;
	
	my $numSeq = @seqAry;
	my $title = $fileTag."_n=$numSeq";
	
	my $weblogoEntropyCMD = "weblogo -f $seqPath -F txt -s large --title $title -A rna -c classic >$entropyPath";
	my $grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
	runAndCheckSerialTask($grepCMD, "weblogo", $weblogoEntropyCMD, "$outFolder/error.log.txt");

	$weblogoEntropyCMD = "weblogo -f $seqPath -F pdf -s large --title $title -A rna -c classic >$logoBitPath";
	$grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
	runAndCheckSerialTask($grepCMD, "weblogo", $weblogoEntropyCMD, "$outFolder/error.log.txt");

	$weblogoEntropyCMD = "weblogo -f $seqPath -F pdf -s large --title $title -A rna --units probability -c classic >$logoProbPath";
	$grepCMD = "ps -ef | grep weblogo | grep -v grep | grep -v perl";
	runAndCheckSerialTask($grepCMD, "weblogo", $weblogoEntropyCMD, "$outFolder/error.log.txt");

}
########################################################################## runAndCheckSerialTask
sub runAndCheckSerialTask {

	my $grepCmd = $_[0];
	my $grepStr = $_[1];
	my $cmd = $_[2];
	my $errorLogPath = $_[3];

	system (qq|$cmd 2>>$errorLogPath &|);
	my $sdout = $grepStr;
	while ($sdout =~ m/$grepStr/) {
		$sdout = `$grepCmd`;
		sleep (0.001);
	}
}
########################################################################## MEMEPolyASeq
sub MEMEPolyASeq {
	
	#---if all search mode is turn off, only base composition will be plot
	
	my %seqForMEMEHsh = %{$_[0]};
	my %cntgSeqSSHsh = %{$_[1]};
	my $maxk = $_[2];
	my $mink = $_[3];
	my $minE = $_[4];
	my $fullShuffleSearch = $_[5];
	my $randGenomicSearch = $_[6];
	my $halfHalfSearch = $_[7];
	my $posSpecificSearch = $_[8];
	my $preDefinedDremeXMLPathToRead = $_[9];
	my $alignmentLength = $_[10];
	
	my ($posSeqRngStart, $posSeqRngEnd, $posRngTag);
	my ($negSeqRngStart, $negSeqRngEnd, $negRngTag);
	my $posSpecificSearchRandGenomicFastaPath = "$outDir/meme/posSpecific.rand.fasta";

	#---generate only once
	if ($posSpecificSearch ne "no") {
		my @posSpecificSearchSplt = split /:/, $posSpecificSearch;
		($posSeqRngStart, $posSeqRngEnd) = split /,/, $posSpecificSearchSplt[0];
		($negSeqRngStart, $negSeqRngEnd) = split /,/, $posSpecificSearchSplt[1];
		$posRngTag = "pos.$posSeqRngStart.$posSeqRngEnd";
		$negRngTag = "neg.$negSeqRngStart.$negSeqRngEnd";
		my $randLen = $posSeqRngEnd - $posSeqRngStart;
		randomGenDNASeq(\%cntgSeqSSHsh, $randLen, 1000, $posSpecificSearchRandGenomicFastaPath);
	}
	
	#----scan for the same thing in random Genomic dataset
	if ($preDefinedDremeXMLPathToRead ne "no") {
		system "mkdir -p -m 777 $outDir/meme/randGenomic/";
		my $fullSeqForDremePath = "$outDir/meme/randGenomic/full.randGenomic.fasta";
		my $dirForMotifOccurPath = "$outDir/meme/randGenomic/preDefinedDremeXML/";
		randomGenDNASeq(\%cntgSeqSSHsh, $alignmentLength, 1000, $fullSeqForDremePath);
		plotBaseComposition("$outDir/meme/randGenomic/", $fullSeqForDremePath, "randGenomic");
		plotDremeMotif($dirForMotifOccurPath, $preDefinedDremeXMLPathToRead, $fullSeqForDremePath, "randGenomic", "randGenomic preDefinedDremeXML");
	}
	
	foreach my $fileTag (keys %seqForMEMEHsh) {
		system "mkdir -p -m 777 $outDir/meme/$fileTag/";
		my $fullSeqForDremePath = "$outDir/meme/$fileTag/full.$fileTag.fasta";
		open (SEQDREME, ">$fullSeqForDremePath");
		foreach my $gene (keys %{$seqForMEMEHsh{$fileTag}}) {
			print SEQDREME ">".$gene."\n";
			print SEQDREME ${$seqForMEMEHsh{$fileTag}}{$gene}."\n";
		}
		close SEQDREME;
		plotBaseComposition("$outDir/meme/$fileTag/", $fullSeqForDremePath, $fileTag);

		if ($preDefinedDremeXMLPathToRead eq "no") {
			
			if ($fullShuffleSearch ne "no") {
				#---using shuffle as negetive
				print "Running dreme for $fileTag.\n";
				my $dirForDremePath = "$outDir/meme/$fileTag/shuffle/";
				system "mkdir -p -m 777 $dirForDremePath";
			
				my $dremeCMD = "dreme -oc $dirForDremePath -p $fullSeqForDremePath -maxk $maxk -mink $mink -e $minE";
				my $grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
				my $grepStr = "dreme -oc $dirForDremePath";
				runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $fullSeqForDremePath, $fileTag, "$fileTag fullShuffleSearch");
			}
			
			if ($posSpecificSearch ne "no") {
			
				#---e.g. 40,70:0,30
				my $posSeqForDremePath = "$outDir/meme/$fileTag/$posRngTag.$fileTag.fasta";
				my $negSeqForDremePath = "$outDir/meme/$fileTag/$negRngTag.$fileTag.fasta";
				open (POSSEQ, ">$posSeqForDremePath");
				open (NEGSEQ, ">$negSeqForDremePath");
				foreach my $gene (keys %{$seqForMEMEHsh{$fileTag}}) {
					my $fullSeq = ${$seqForMEMEHsh{$fileTag}}{$gene};
					my $posSeq = substr $fullSeq, $posSeqRngStart, ($posSeqRngEnd-$posSeqRngStart);
					my $negSeq = substr $fullSeq, $negSeqRngStart, ($negSeqRngEnd-$negSeqRngStart);
					print POSSEQ ">".$gene."\n";
					print POSSEQ $posSeq."\n";
					print NEGSEQ ">".$gene."\n";
					print NEGSEQ $negSeq."\n";
				}
				close POSSEQ;
				close NEGSEQ;
				my $rngTag = "pos.$posSeqRngStart.$posSeqRngEnd.neg.$negSeqRngStart.$negSeqRngEnd";
				my $dirForDremePath = "$outDir/meme/$fileTag/pSrchNeg/";
				system "mkdir -p -m 777 $dirForDremePath";
				my $dremeCMD = "dreme -oc $dirForDremePath -p $posSeqForDremePath -n $negSeqForDremePath -maxk $maxk -mink $mink -e $minE";
				my $grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
				my $grepStr = "dreme -oc $dirForDremePath";
				runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $posSeqForDremePath, "pos.".$fileTag, "positive $fileTag neg $rngTag");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $negSeqForDremePath, "neg.".$fileTag, "negative $fileTag neg $rngTag");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $fullSeqForDremePath, "full.".$fileTag, "full $fileTag neg $rngTag");
	
				$dirForDremePath = "$outDir/meme/$fileTag/pSrchRand/";
				system "mkdir -p -m 777 $dirForDremePath";
				$dremeCMD = "dreme -oc $dirForDremePath -p $posSeqForDremePath -n $posSpecificSearchRandGenomicFastaPath -maxk $maxk -mink $mink -e $minE";
				$grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
				$grepStr = "dreme -oc $dirForDremePath";
				runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $posSeqForDremePath, "pos.".$fileTag, "positive $fileTag rand $rngTag");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $posSpecificSearchRandGenomicFastaPath, "rand.".$fileTag, "rand $fileTag rand $rngTag");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $fullSeqForDremePath, "full.".$fileTag, "full $fileTag rand $rngTag");
	
				$dirForDremePath = "$outDir/meme/$fileTag/pSrchShuffle/";
				system "mkdir -p -m 777 $dirForDremePath";
				$dremeCMD = "dreme -oc $dirForDremePath -p $posSeqForDremePath -maxk $maxk -mink $mink -e $minE";
				$grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
				$grepStr = "dreme -oc $dirForDremePath";
				runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $posSeqForDremePath, "pos.".$fileTag, "positive $fileTag shuffle $rngTag");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $fullSeqForDremePath, "full.".$fileTag, "full $fileTag shuffle $rngTag");
			}
			
			#---using randomGenomicDNA as negetive, if any
			if ($randGenomicSearch ne "no") {
				
				#----- find motifs
				my $randNum = 1000;
				my $randGenomicSeqFastaPath = "$outDir/meme/rand_num$randNum"."_len$alignmentLength.fasta";
	
				randomGenDNASeq($cntgSeqSSHsh_ref, $alignmentLength, $randNum, $randGenomicSeqFastaPath);
				
				my $dirForDremePath = "$outDir/meme/$fileTag/randGnmc/";
				system "mkdir -p -m 777 $dirForDremePath";
				open (SEQDREME, ">$fullSeqForDremePath");
				foreach my $gene (keys %{$seqForMEMEHsh{$fileTag}}) {
					print SEQDREME ">".$gene."\n";
					print SEQDREME ${$seqForMEMEHsh{$fileTag}}{$gene}."\n";
				}
				close SEQDREME;
				
				my $dremeCMD = "dreme -oc $dirForDremePath -p $fullSeqForDremePath -n $randGenomicSeqFastaPath -maxk $maxk -mink $mink -e $minE";
				my $grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
				my $grepStr = "dreme -oc $dirForDremePath";
				runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $fullSeqForDremePath, "full.".$fileTag, "full $fileTag randomGenomic");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $randGenomicSeqFastaPath, "rand.".$fileTag, "rand $fileTag randomGenomic");
			}
			
			if ($halfHalfSearch ne "no") {
				my $firstHalfSeqForDremePath = "$outDir/meme/$fileTag/1stHf.$fileTag.fasta";
				my $secondHalfSeqForDremePath = "$outDir/meme/$fileTag/2ndHf.$fileTag.fasta";
				open (FISRTHALF, ">$firstHalfSeqForDremePath");
				open (SECONDHALF, ">$secondHalfSeqForDremePath");
				foreach my $gene (keys %{$seqForMEMEHsh{$fileTag}}) {
					my $fullSeq = ${$seqForMEMEHsh{$fileTag}}{$gene};
					my $seqLenHalf = int ((length $fullSeq)/2);
					my $firstHalfSeq = substr $fullSeq, 0, $seqLenHalf;
					my $secondHalfSeq = substr $fullSeq, $seqLenHalf;
					print FISRTHALF ">".$gene."\n";
					print FISRTHALF $firstHalfSeq."\n";
					print SECONDHALF ">".$gene."\n";
					print SECONDHALF $secondHalfSeq."\n";
				}
				close FISRTHALF;
				close SECONDHALF;
				
				my $dirForDremePath = "$outDir/meme/$fileTag/1stHf/";
				system "mkdir -p -m 777 $dirForDremePath";
				my $dremeCMD = "dreme -oc $dirForDremePath -p $firstHalfSeqForDremePath -n $secondHalfSeqForDremePath -maxk $maxk -mink $mink -e $minE";
				my $grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
				my $grepStr = "dreme -oc $dirForDremePath";
				runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $firstHalfSeqForDremePath, "1stHf.".$fileTag, "1stHf $fileTag 1stHf vs 2ndHf");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $secondHalfSeqForDremePath, "2ndHf.".$fileTag, "2ndHf $fileTag 1stHf vs 2ndHf");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $fullSeqForDremePath, "full.".$fileTag, "full $fileTag 1stHf vs 2ndHf");
	
				$dirForDremePath = "$outDir/meme/$fileTag/2ndHf/";
				system "mkdir -p -m 777 $dirForDremePath";
				$dremeCMD = "dreme -oc $dirForDremePath -p $secondHalfSeqForDremePath -n $firstHalfSeqForDremePath -maxk $maxk -mink $mink -e $minE";
				$grepCMD = "ps -ef | grep dreme | grep -v grep | grep -v perl";
				$grepStr = "dreme -oc $dirForDremePath";
				runAndCheckSerialTask($grepCMD, $grepStr, $dremeCMD, "$outDir/error.log.txt");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $firstHalfSeqForDremePath, "1stHf.".$fileTag, "1stHf $fileTag 2ndHf vs 1stHf");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $secondHalfSeqForDremePath, "2ndHf.".$fileTag, "2ndHf $fileTag 2ndHf vs 1stHf");
				plotDremeMotif($dirForDremePath, $dirForDremePath."./dreme.xml", $fullSeqForDremePath, "full.".$fileTag, "full $fileTag 2ndHf vs 1stHf");
			}
		} else {

			my $dirForMotifOccurPath = "$outDir/meme/$fileTag/preDefinedDremeXML/";
			plotDremeMotif($dirForMotifOccurPath, $preDefinedDremeXMLPathToRead, $fullSeqForDremePath, $fileTag, "$fileTag preDefinedDremeXML");

		}
	}
}
########################################################################## randomGenDNASeq
sub randomGenDNASeq {
	
	#---my $randSeqHsh_ref = randomGenDNASeq(\%refSeqHsh, $seqLen, $seqNum, $fastaPath);
	my %refSeqHsh = %{$_[0]};
	my $seqLen = $_[1];
	my $seqNum = $_[2];
	my $fastaPath = $_[3];#---input the path for printing or use "no" to swicth off;
	
	print "Generating $seqNum random sequences of $seqLen in length.\n";
	
	open (FASTA, ">$fastaPath") if ($fastaPath ne "no");
	
	my %randSeqHsh;
	my $conCatSeq = "";
	
	foreach my $seqName (keys %refSeqHsh) {
		my $seq = ${$refSeqHsh{$seqName}}{"+"};
		$conCatSeq .= $seq;
	}
	
	my $totalConCatLength = length $conCatSeq;
	
	for my $round (1..$seqNum) {
		my $validSeq = "no";
		my $strnd = "+";
		$strnd = "-" if (rand() <= 0.5);
		while ($validSeq eq "no") {#---loop until not out of ctng rng

			my $seqStartPos = int (rand $totalConCatLength); #---rand from 0 to $cntgLen-1;
			next if (($seqStartPos + $seqLen) > $totalConCatLength);
			my $extractSeq = substr $conCatSeq, $seqStartPos, $seqLen;
		 	next if ($extractSeq =~ m/[^ATGCatgc]/);
		 	$validSeq = "yes";
		 	
		 	if ($strnd eq "-") {
		 		$extractSeq = reverse $extractSeq;
		 		$extractSeq =~ tr/ACGTacgt/TGCAtgca/;
		 	}
		 	
		 	$randSeqHsh{"randSeq".$round} = $extractSeq;
		 	
		 	if ($fastaPath ne "no") {
			 	print FASTA ">"."randSeq".$round."\n";
			 	print FASTA $extractSeq."\n";
			}
		}
	}
	
	return \%randSeqHsh;
}
########################################################################## randomGenDNASeq
sub getDremeMotif {

	my $dremeXMLToRead = $_[0];
	
	my %motifInfoHsh;
	my ($motifSeq, $evalue);
	open (DREMEXML, "$dremeXMLToRead");
	while (my $theLine = <DREMEXML>) {
		chomp $theLine;
		if ($theLine =~ m/<motif id=/) {
			my @theLineSplt = split / |\>/, $theLine;
			foreach my $arg (@theLineSplt) {
				$arg =~ s/\"//g;
				if ($arg =~ m/^seq=/) {$motifSeq = substr ($arg, index ($arg, "=")+1);}
				elsif ($arg =~ m/^evalue=/) {$evalue = substr ($arg, index ($arg, "=")+1);}
			}
			${$motifInfoHsh{$motifSeq}}{"evalue"} = $evalue;
		}
		if ($theLine =~ m/<match seq=/) {
			my @theLineSplt = split / |\>/, $theLine;
			foreach my $arg (@theLineSplt) {
				if ($arg =~ m/^seq=/) {
					$arg =~ s/\"//g;
					my $matchSeq = substr ($arg, index ($arg, "=")+1); 
					push @{${$motifInfoHsh{$motifSeq}}{"matchSeq"}}, $matchSeq;
				}
			}
		}
	}
	close (DREMEXML);
	
	foreach my $motifSeq (keys %motifInfoHsh) {
		my %tmpPermutatedSeqHsh;
		foreach my $matchSeq (@{${$motifInfoHsh{$motifSeq}}{"matchSeq"}}) {
			#---permutation 1bp mismatch
			my @matchSeqSplt = split //, $matchSeq;
			foreach my $pos (0..$#matchSeqSplt) {
				foreach my $mutant (("A", "T", "G", "C")) {
					my @permutatedSeqAry = @matchSeqSplt;
					$permutatedSeqAry[$pos] = $mutant;
					my $permutatedSeqStr = join "", @permutatedSeqAry;
					push @{${$motifInfoHsh{$motifSeq}}{"permutateSeq"}}, $permutatedSeqStr if (not exists $tmpPermutatedSeqHsh{$permutatedSeqStr});
					$tmpPermutatedSeqHsh{$permutatedSeqStr}++;
				}
			}
		}
	}

	
	my $motifNum = keys %motifInfoHsh;
	
	print "$motifNum motif stored.\n";
	
	return \%motifInfoHsh;

}
########################################################################## plotDremeMotif
sub plotDremeMotif {
	
	#---subroutine dependency: getDremeMotif, readMultiFasta, GNUPlotMultipleColumnXYLines, getSeqComposition, sortoutUnmatchedSeq
	#---in/out: plotDremeMotif($dremeXMLToRead, $fastaSeq);
	
	my $subOutDir = $_[0];
	my $dremeXMLToRead = $_[1];
	my $fastaSeq = $_[2];
	my $fileTag = $_[3];
	my $title = $_[4];
	
	system "mkdir -p -m 777 $subOutDir/";

	my $motifInfoHsh_ref = getDremeMotif($dremeXMLToRead);
	my %motifInfoHsh = %{$motifInfoHsh_ref};
	if ((keys %motifInfoHsh) >= 1) {
		my $fastaHsh_ref = readMultiFasta($fastaSeq);
		my $motifOccurenceFilePath = $subOutDir."/$fileTag.motifOccurence.permutateSeq.txt";
		my $motifOccurencePlotFilePath = $subOutDir."/$fileTag.motifOccurence.permutateSeq.pdf";
		my ($pmtHitMotifCountHsh_ref, $pmtSeqHitInfoHsh_ref) = scanForMotifOccurence($motifInfoHsh_ref, $fastaHsh_ref, $motifOccurenceFilePath, "permutateSeq", "yes", $targetRngStart, $targetRngEnd, $subOutDir);
		GNUPlotMultipleColumnXYLines($motifOccurencePlotFilePath, $motifOccurenceFilePath, $title, "percentage");
		$motifOccurenceFilePath = $subOutDir."/$fileTag.motifOccurence.matchSeq.txt";
		$motifOccurencePlotFilePath = $subOutDir."/$fileTag.motifOccurence.matchSeq.pdf";
		my ($mchHitMotifCountHsh_ref, $mchSeqHitInfoHsh_ref) =  scanForMotifOccurence($motifInfoHsh_ref, $fastaHsh_ref, $motifOccurenceFilePath, "matchSeq", "yes", $targetRngStart, $targetRngEnd, $subOutDir);
		GNUPlotMultipleColumnXYLines($motifOccurencePlotFilePath, $motifOccurenceFilePath, $title, "percentage");
		
		#---sort of the matched and non-matched sequences
		
		
	} else {
		open TMP, ">$subOutDir/$fileTag.has.no.significant.motifs.txt"; close TMP;
	}
	
	return;
}
########################################################################## plotDremeMotif
sub plotBaseComposition {
	
	#---subroutine dependency: getDremeMotif, readMultiFasta, GNUPlotMultipleColumnXYLines, getSeqComposition
	#---in/out: plotDremeMotif($dremeXMLToRead, $fastaSeq);
	
	my $subOutDir = $_[0];
	my $fastaSeq = $_[1];
	my $fileTag = $_[2];
	
	system "mkdir -p -m 777 $subOutDir/";

	my $fastaHsh_ref = readMultiFasta($fastaSeq);
	my $seqNum = keys %{$fastaHsh_ref};
	my $seqComPath = $subOutDir."/$fileTag.seqComposition.txt";
	getSeqComposition($fastaHsh_ref, $seqComPath);
	my $baseComPlotFilePath = $subOutDir."/$fileTag.seqComposition.pdf";
	GNUPlotMultipleColumnXYLines($baseComPlotFilePath, $seqComPath, $fileTag." n=$seqNum", "percentage");
	return;
}
########################################################################## getSeqComposition
sub getSeqComposition {
	
	#---subroutine dependency:
	#---in/out: getSeqComposition($fastaHsh);
	
	my %fastaHsh = %{$_[0]};
	my $outputPath = $_[1]; #---use no to switch off
	
	my %posCountPctHsh;
	my %posTotalHsh;
	my %allResHsh;
	
	
	foreach my $seqName (keys %fastaHsh) {
		my $seq = $fastaHsh{$seqName};
		my @seqSplt = split //, $seq;
		my $pos = 0;
		foreach my $res (@seqSplt) {
			$allResHsh{$res}++;
			$pos++;
			${${$posCountPctHsh{$pos}}{$res}}{"count"}++;
			$posTotalHsh{$pos}++;
		}
	}
	
	foreach my $pos (keys %posTotalHsh) {
		my $totalCount = $posTotalHsh{$pos};
		foreach my $res (sort {$a cmp $b} keys %allResHsh) {
			my $pct = 0;
			if (exists ${${$posCountPctHsh{$pos}}{$res}}{"count"}) {
				$pct = sprintf "%.06f", 100*${${$posCountPctHsh{$pos}}{$res}}{"count"}/$totalCount;
			}
			${${$posCountPctHsh{$pos}}{$res}}{"pct"} = $pct;
		}
	}
	
	if ($outputPath ne "no") {
		open (OUTSEQCOM, ">$outputPath");
		print OUTSEQCOM "pos";
		foreach my $res (sort {$a cmp $b} keys %allResHsh) {
			print OUTSEQCOM "\t".$res;
		}
		print OUTSEQCOM "\n";

		foreach my $pos (sort {$a <=> $b} keys %posCountPctHsh) {
			print OUTSEQCOM $pos;
			foreach my $res (sort {$a cmp $b} keys %{$posCountPctHsh{$pos}}) {
				print OUTSEQCOM "\t".${${$posCountPctHsh{$pos}}{$res}}{"pct"};
			}
			print OUTSEQCOM "\n";
		}
		close OUTSEQCOM;
	}
	
	return (\%posCountPctHsh, \%posTotalHsh, \%allResHsh);
}
########################################################################## GNUPlotMultipleColumnXYLines
sub GNUPlotMultipleColumnXYLines {

	#---GNUPlotMultipleColumnXYLines($plotFilePath, $plotDataPath, $title, $ylabel);
	#The file must look like this
	#X	Y1	Y2	Y3
	#1	45	872	68
	#2	45	87	53
	#3	68	97	1
	#4	10	2	60
	
	my $plotFilePath = $_[0];
	my $plotDataPath = $_[1];
	my $title = $_[2];
	my $ylabel = $_[3];
	my $extraCmd = $_[4];
	
	$extraCmd = "" if (not defined $extraCmd);
	
	$plotFilePath .= ".pdf" if ($plotFilePath !~ m/\.pdf$/);

	my @filePathSplt = split /\//, $plotFilePath;
	my $fileName = $filePathSplt[-1];

	my @allPlotCmdAry;
	my @indivPlotCmdAry;
	my $xlabel;
	open PLOTDATA, "$plotDataPath" ;
	while (my $theLine = <PLOTDATA>) {
		chomp $theLine;
		my @theLineSplt = split /\t/, $theLine;
		$xlabel = $theLineSplt[0];
		foreach my $i (1..$#theLineSplt) {
			my $legend = $theLineSplt[$i];
			my $YColumn = $i+1;
			push @allPlotCmdAry, "\'$plotDataPath\' using 1:$YColumn with lines title \"$legend\" ";
			push @indivPlotCmdAry, "plot \'$plotDataPath\' using 1:$YColumn with lines title \"$legend\"; ";
		}
		last;
	}
	close PLOTDATA;
	
	my $allPlotCmdstr = join ",", @allPlotCmdAry;
	my $indivPlotCmdStr = join " ", @indivPlotCmdAry;

	$allPlotCmdstr = "plot ".$allPlotCmdstr.";";
	
	print "Running GNUPlotMultipleColumnXYLines for $fileName.\n";
	
	#---do the GNUPLOT
	open (GNUPLOT, "|gnuplot");
	print GNUPLOT <<EOPLOT;
	set terminal postscript color solid
	set output "| ps2pdf - $plotFilePath 2>/dev/null";
	unset logscale x;
	unset logscale y;
	$extraCmd;
	set xlabel "$xlabel";
	set ylabel "$ylabel";
	set title "$title";
	$allPlotCmdstr;
	$indivPlotCmdStr;
EOPLOT
	close(GNUPLOT);
}
########################################################################## scanForMotifOccurence
sub scanForMotifOccurence {
	
	#---subroutine dependency: 
	#---in/out: scanForMotifOccurence($var);
	
	my %motifInfoHsh = %{$_[0]}; #---push @{${$motifInfoHsh{$motifSeq}}{"matchSeq"}}, $matchSeq;
	my %fastaHsh = %{$_[1]};
	my $outFilePath = $_[2];
	my $permutateOrMatchSeq = $_[3];
	my $fullRngHit = $_[4]; #---yes or no, if yes, will count full rng hit rather than start pos only
	my $rngStart = $_[5];
	my $rngEnd = $_[6];
	my $subOutDir = $_[7];
	
	my %motifHitHsh;
	my %seqHitInfoHsh;
	my $totalSeqNum = keys %fastaHsh;
	my %allPosHsh; #---just to store the position;
	my %hitMotifCountHsh;
	
	foreach my $motifSeq (sort {${$motifInfoHsh{$a}}{"evalue"} <=> ${$motifInfoHsh{$b}}{"evalue"}} keys %motifInfoHsh) {
		print "Scanning motif occurence for $motifSeq\r";
		%allPosHsh = ();
		foreach my $seqName (keys %fastaHsh) {
			my %seqHitPosHsh = ();
			my $seq = $fastaHsh{$seqName};
			$seq =~ tr/Uu/Tt/;
			
			for my $pos (0..((length $seq)-1)) {
				${$seqHitPosHsh{$seqName}}{$pos} = 0;
				$allPosHsh{$pos}++;
				${${$motifHitHsh{$motifSeq}}{$pos}}{"count"} = 0 if (not exists ${${$motifHitHsh{$motifSeq}}{$pos}}{"count"}); #---build the motifHitHsh base
			}

			foreach my $seqToSearch (@{${$motifInfoHsh{$motifSeq}}{$permutateOrMatchSeq}}) {
				my $motifLen = length $motifSeq;
				my @startPosAry;
				my $offset = 0;
				my $startPos = index(uc ($seq), $seqToSearch, $offset);#---UC for upper case as dreme always output uppercase
				while ($startPos != -1) {
					push @startPosAry, $startPos;
					push @{${$seqHitInfoHsh{$seqName}}{$motifSeq}}, ($seqToSearch, $startPos);
					push @{${${$hitMotifCountHsh{"full"}}{$motifSeq}}{$seqName}}, $seqToSearch;
					push @{${${$hitMotifCountHsh{"withinRng"}}{$motifSeq}}{$seqName}}, $seqToSearch if (($startPos >= $rngStart) and ($startPos <= $rngEnd));

					$offset = $startPos + 1;
					$startPos = index(uc ($seq), $seqToSearch, $offset);
				}
				
				foreach my $startPos (@startPosAry) {
					if ($fullRngHit eq "yes") {
						foreach my $hitPos ($startPos..($startPos+$motifLen-1)) {
							${$seqHitPosHsh{$seqName}}{$hitPos}++;
						}
					} else {
						${$seqHitPosHsh{$seqName}}{$startPos}++;
					}
				}
			} #---end foreach my $matchSeq (@{${$motifInfoHsh{$motifSeq}}{"matchSeq"}}) {
			
			foreach my $hitPos (keys %{$seqHitPosHsh{$seqName}}) {
				${${$motifHitHsh{$motifSeq}}{$hitPos}}{"count"}++ if (${$seqHitPosHsh{$seqName}}{$hitPos} > 0);
			}
		}#--- end of foreach my $seqName (keys %fastaHsh) {
		
		foreach my $pos (keys %{$motifHitHsh{$motifSeq}}) {
			my $count = ${${$motifHitHsh{$motifSeq}}{$pos}}{"count"};
			my $pct = sprintf "%.06f", 100*$count/$allPosHsh{$pos};
			${${$motifHitHsh{$motifSeq}}{$pos}}{"pct"} = $pct;
		}
	}#---end of foreach my $motifSeq (sort {${$motifInfoHsh{$a}}{"evalue"} <=> ${$motifInfoHsh{$b}}{"evalue"}} keys %motifInfoHsh) {
	print "......finished\n";

	open (MOTIFOCC, ">$outFilePath");
	print MOTIFOCC "pos";
	foreach my $motifSeq (sort {${$motifInfoHsh{$a}}{"evalue"} <=> ${$motifInfoHsh{$b}}{"evalue"}} keys %motifInfoHsh) {
		my $withinRngHit = keys %{${$hitMotifCountHsh{"withinRng"}}{$motifSeq}};
		my $fullHit = keys %{${$hitMotifCountHsh{"full"}}{$motifSeq}};
		print MOTIFOCC "\t".$motifSeq."_full:$fullHit/$totalSeqNum"."_rng:$withinRngHit/$totalSeqNum"."_".${$motifInfoHsh{$motifSeq}}{"evalue"};
	}
	print MOTIFOCC "\n";
	foreach my $pos (sort {$a <=> $b} keys %allPosHsh) {
		print MOTIFOCC $pos;
		foreach my $motifSeq (sort {${$motifInfoHsh{$a}}{"evalue"} <=> ${$motifInfoHsh{$b}}{"evalue"}} keys %motifInfoHsh) {
			print MOTIFOCC "\t".${${$motifHitHsh{$motifSeq}}{$pos}}{"pct"};
		}
		print MOTIFOCC "\n";
	}
	close MOTIFOCC;
	
	#---sortout the macth and unmatched seq
	foreach my $fullOrRng (("full", "withinRng")) {
		foreach my $motifSeq (keys %{$hitMotifCountHsh{$fullOrRng}}) {
			my %tmpMotifCountHsh;
			my %tmpMotifCountTotalHsh;
			system "mkdir -p -m 777 $subOutDir/sortSeq/$permutateOrMatchSeq/$fullOrRng/$motifSeq/";
			system "mkdir -p -m 777 $subOutDir/motifCount/$permutateOrMatchSeq/$fullOrRng/$motifSeq/";
			system "mkdir -p -m 777 $subOutDir/motifSeqFasta/$permutateOrMatchSeq/$fullOrRng/$motifSeq/";
			open (MATCHSEQ, ">$subOutDir/sortSeq/$permutateOrMatchSeq/$fullOrRng/$motifSeq/matched.fasta");
			open (UNMATCHSEQ, ">$subOutDir/sortSeq/$permutateOrMatchSeq/$fullOrRng/$motifSeq/unmatched.fasta");
			open (MOTIFSEQFASTA, ">$subOutDir/motifSeqFasta/$permutateOrMatchSeq/$fullOrRng/$motifSeq/matchedMotifSeq.fasta");
			open (MOTIFSEQCOUNT, ">$subOutDir/motifCount/$permutateOrMatchSeq/$fullOrRng/$motifSeq/motifSeqCount.txt");

			foreach my $seqName (keys %fastaHsh) {
				if (exists ${${$hitMotifCountHsh{$fullOrRng}}{$motifSeq}}{$seqName}) {
					foreach my $seqToSearch (@{${${$hitMotifCountHsh{$fullOrRng}}{$motifSeq}}{$seqName}}) {
						${$tmpMotifCountHsh{$seqToSearch}}{$seqName}++;
						$tmpMotifCountTotalHsh{$seqToSearch}++;
						print MOTIFSEQFASTA ">".$seqName."_".$seqToSearch."_".${$tmpMotifCountHsh{$seqToSearch}}{$seqName}."\n";
						print MOTIFSEQFASTA $seqToSearch."\n";
					}
					print MATCHSEQ ">".$seqName."\n";
					print MATCHSEQ $fastaHsh{$seqName}."\n";
				} else {
					print UNMATCHSEQ ">".$seqName."\n";
					print UNMATCHSEQ $fastaHsh{$seqName}."\n";
				}
			}
			
			foreach my $seqToSearch (sort {$tmpMotifCountTotalHsh{$b} <=> $tmpMotifCountTotalHsh{$a}} keys %tmpMotifCountTotalHsh) {
				my $hitNumSeq = keys %{$tmpMotifCountHsh{$seqToSearch}};
				my $count = 0;
				foreach my $seqName (keys %{$tmpMotifCountHsh{$seqToSearch}}) {
					$count += ${$tmpMotifCountHsh{$seqToSearch}}{$seqName};
				}
				my $pct = sprintf "%.02f", 100*$hitNumSeq/$totalSeqNum;
				my $hitPerSeq = sprintf "%.02f", $count/$hitNumSeq;
				print MOTIFSEQCOUNT $seqToSearch."\t".$count."\t".$pct."\t".$hitPerSeq."\n";
			}
			close MATCHSEQ;
			close UNMATCHSEQ;
			close MOTIFSEQCOUNT;
			close MOTIFSEQFASTA;
		}
	}

	return \%hitMotifCountHsh, \%seqHitInfoHsh;
}
########################################################################## calculateArySDMeanMedianQuartile
sub calculateArySDMeanMedianQuartile {

	my @numAry = @{$_[0]};
	
	my $totalNum = @numAry;

	my ($mean, $stdDev, $max, $pct99, $pct95, $pct90, $pct75, $median, $pct25, $pct10, $pct5, $pct1, $min);

	if ($totalNum >= 10) {
	
		my @sortedNumAry = sort {$a <=> $b} @numAry;

		# Step 1, find the mean of the numbers
		my $total1 = 0;
		foreach my $num (@sortedNumAry) {
			$total1 += $num;
		}
		$mean = $total1 / $totalNum;
	
		# Step 2, find the mean of the squares of the differences
		# between each number and the mean
		my $total2 = 0;
		foreach my $num (@sortedNumAry) {
			$total2 += ($mean-$num)**2;
		}
		my $mean2 = $total2 / $totalNum;
	
		# Step 3, standard deviation is the square root of the
		# above mean
		$stdDev = sqrt($mean2);

		my $maxIndex = $#sortedNumAry; $max = $sortedNumAry[$maxIndex];
		my $pct99Index = int ($totalNum*0.99); $pct99 = $sortedNumAry[$pct99Index];
		my $pct95Index = int ($totalNum*0.95); $pct95 = $sortedNumAry[$pct95Index];
		my $pct90Index = int ($totalNum*0.90); $pct90 = $sortedNumAry[$pct90Index];
		my $pct75Index = int ($totalNum*0.75); $pct75 = $sortedNumAry[$pct75Index];
		my $medianIndex = int ($totalNum*0.5); $median = $sortedNumAry[$medianIndex];
		my $pct25Index = int ($totalNum*0.25); $pct25 = $sortedNumAry[$pct25Index];
		my $pct10Index = int ($totalNum*0.1); $pct10 = $sortedNumAry[$pct10Index];
		my $pct5Index = int ($totalNum*0.05); $pct5 = $sortedNumAry[$pct5Index];
		my $pct1Index = int ($totalNum*0.01); $pct1 = $sortedNumAry[$pct1Index];
		my $minIndex = 0; $min = $sortedNumAry[$minIndex];
	} else {
	 	$totalNum = $mean = $stdDev = $max = $pct99 = $pct95 = $pct90 = $pct75 = $median = $pct25 = $pct10 = $pct5 = $pct1 = $min = 0;
	}
	#0=$totalNum, 1=$mean, 2=$stdDev, 3=$max, 4=$pct99, 5=$pct95, 6=$pct90, 7=$pct75, 8=$median, 9=$pct25, 10=$pct10, 11=$pct5, 12=$pct1, 13=$min);
	return ($totalNum, $mean, $stdDev, $max, $pct99, $pct95, $pct90, $pct75, $median, $pct25, $pct10, $pct5, $pct1, $min);
}
########################################################################## readInDirForFasta
sub readInDirForFasta {
	
	#---subroutine dependency: 
	#---in/out: readInDirForFasta($inDir);
	
	my $inDirToRead = $_[0];
	print "Reading $inDirToRead\n";
	
	opendir (DIR, $inDirToRead);
	my @allFastaAry = grep /\.fasta$/, readdir DIR;
	closedir DIR;
	
	my %seqForMEMEHsh;
	my %aligmentLenHsh;
	my $alignmentLength;
	foreach my $fastaFileName (@allFastaAry) {
		my $storeSeqNum = 0;
		my $fileTag = $fastaFileName;
		$fileTag =~ s/\.fasta$//;
		next if (($fileTag !~ m/$fileFilterTag/) and ($fileFilterTag ne ""));
		print "Reading sequence from $fileTag\n";
		open (FASTA, "$inDirToRead/$fastaFileName");
		my ($gene, $seq);
		my @tmpSeqAry;
		my @tmpGeneAry;
		while (my $theLine = <FASTA>) {
			chomp $theLine;
			if ($theLine =~ m/^>/) {
				$theLine =~ s/^>//;
				$gene = $theLine;
			} elsif ((length $theLine) > 2) {
				$seq = $theLine;
				push @tmpSeqAry, $seq;
				push @tmpGeneAry, $gene;
				my $length = length $seq;
				$alignmentLength = $length;
				$aligmentLenHsh{$length}++;
				$storeSeqNum++;
			}
		}
		close FASTA;

		my $arySize = @tmpSeqAry;
		if ($storeSeqNum >= $maxSeqNum) {
			print "storeSeqNum > $maxSeqNum. Random Sampling Sequence\n";
			foreach (1..$maxSeqNum) {
				$arySize = @tmpSeqAry;
				my $randIndex = int (rand $arySize);
				my $gene = $tmpGeneAry[$randIndex];
				my $seq = $tmpSeqAry[$randIndex];
				${$seqForMEMEHsh{$fileTag}}{$gene} = $seq;
				splice(@tmpGeneAry, $randIndex, 1);
				splice(@tmpSeqAry, $randIndex, 1);
			}
		} else {
			foreach my $index (0..($arySize-1)) {
				my $gene = $tmpGeneAry[$index];
				my $seq = $tmpSeqAry[$index];
				${$seqForMEMEHsh{$fileTag}}{$gene} = $seq;
			}
		}
	}

	my $numLen = keys %aligmentLenHsh;
	if ($numLen > 1) {
		die "Problem with alignment length\n";
	}
	return \%seqForMEMEHsh, $alignmentLength;
}
