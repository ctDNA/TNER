#!/usr/bin/perl
################################################################
# File name: pileup2actg.pl
# This tool processes a pileup file to generate a tab-delimated txt file for the
# TNER main function.
# Last updated: 08/20/2018 
# Contact: shibing.deng {at} pfizer.com & tao.xie {at} pfizer.com or xietao2000 {at} gmail.com
# TNER publication: BMC Bioinformatics (in revision)
# Pfizer Early Clinical Development Biostatistics 
# Pfizer Early Oncology Development and Clinical Research 
#################################################################

#################################################################
# Command line usage: perl pileup2actg.pl test.pileup
#################################################################
#################################################################
# Input:   A pileup file 
# Output:  A tab-delimated txt file with a ".freq" suffix. The fields are: CHR 
#           = chromosome number; POSITION = genomic coordinate position;
#           DEPTH = coverage depth; REF=reference nucleotide;
#           R+ = forward strand coverage for the reference nucleotide; R- = 
#           reverse strand coverage for the reference nucleotide;
#           A+ = forward strand base A count; A- = reverse strand base A count; ... ...
#           Note: only error counts are listed. The reference base count is set to 0.
#           For example, if REF=A, then A+ =0 and A- = 0.
#################################################################
use strict;
use warnings;
if(scalar(@ARGV) != 1){
    print "\nUsage: perl pileup2actg.pl test.pileup\n\n";
    exit(0);
}
open (F1, "< $ARGV[0]" ) || die "Can't find the input file!\n";
open (F2, "> $ARGV[0].freq");

print F2 "CHR\t"."POSITION\t"."DEPTH\t"."REF\t"."R+\t"."R-\t"."A+\t"."A-\t"."C+\t"."C-\t"."T+\t"."T-\t"."G+\t"."g-\n";
while (<F1>) {
        chomp;
        my($chr,$pos,$ref,$depth,$bases,$bq) = split /\s+/;
        $bases=~s/\^.//g;
        $bases=~s/\$//g;
        my @base=split (//,$bases);

        my @f = (); for (my $i=0; $i<10; $i++) {$f[$i]=0;}
        for (my $i=0; $i<@base; $i++) {
                my $s = $base[$i];
                if($s eq '.') {$f[0]++;
                } elsif ($s eq ',') {$f[1]++;
                } elsif ($s eq 'A') {$f[2]++;
                } elsif ($s eq 'a') {$f[3]++;
                } elsif ($s eq 'C') {$f[4]++;
                } elsif ($s eq 'c') {$f[5]++;
                } elsif ($s eq 'T') {$f[6]++;
                } elsif ($s eq 't') {$f[7]++;
                } elsif ($s eq 'G') {$f[8]++;
                } elsif ($s eq 'g') {$f[9]++;
                }
        }
        my $depth_count = 0;
        for(my $i = 0; $i < @f; $i++){
                $depth_count += $f[$i];
        }
        print F2 "$chr\t$pos\t$depth\t$ref\t$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t$f[8]\t$f[9]\n";
}
close(F1);
close(F2);

