#!/usr/bin/perl

#Written by Sarah A. Teichmann on March 6th 1997 and modified on July 2nd 1997.

#Usage: sw_exp_search.pl <query database> <target database> (both in fasta format) <fasta version path>
#Function: to allow a search of a whole family to be run automatically using fasta.

$|=1;

$usage="sw_exp_search.pl <query database> <target database> <fasta version>";

die "\nIncorrect number of arguments.\n\nUsage: $usage\n\ndied" unless ($#ARGV>=1);


my $qdb_main=$ARGV[0];
my $tdb_main=$ARGV[1];
my $fastaver_main=$ARGV[3];

&fasta_kt1_search ($qdb_main, $tdb_main, $fastaver_main);

#------------Subroutine 1: doing a sw search of one database against the other-------#

#_______________________________________________________________________________
# Title     : fasta_kt1_search
# Usage     : &fasta_kt1_search($query_database, $target_database, $fasta_version_to_use
# Function  : to search one database against the other using fasta
#                ktup=1 (default is simply "fasta"). The results are stored in sub dirs
#                which are from the 2 first chars of the query sequence.
# Example   : &fasta_kt1_search ($qdb_main, $tdb_main, $fastaver_main);
# Keywords  : fasta_search, fasta_database_search
# Options   :
# Author    : Sarah A. Teichmann
# Date      : 19th September 1997
# Version   : 1.1
#-------------------------------------------------------------------------------
sub fasta_kt1_search{
    my ($qdb, $tdb, @qdbcont, $fastaver, $gene, $seq,
        @genes, %genes, $genes, $out, $tmp, $sw_score,
        $e_val, @tmpcontent, $i, $dir, @dir);
    $qdb=$_[0];
    open (QDB, "$qdb");
    @qdbcont=<QDB>;
    close QDB;
    $tdb=$_[1];
    open(TDB, "$tdb");
    close TDB;
    if ($_[2]){ $fastaver=$_[2]; }
    else{$fastaver="fasta";}
    for ($i=0; $i<@qdbcont; $i++) {

	my $qdbcont=$qdbcont[$i];
	if ($qdbcont=~/^\>(\S+)/) {
	   $gene=$1;
	   push (@genes, $gene);
	}
	if ($qdbcont=~/^(\w+)/) {
	   $seq=$1;
	   $genes{"$gene"}.="$seq";
       }
	else {next;}
    }
    for ($i=0; $i<@genes; $i++) {
        $genes=$genes[$i];
        @dir=split(//, $genes);
        @dir=splice(@dir, 0, 2);
        $dir=join('', @dir);
        mkdir $dir;
        $out="$dir"."/"."$genes".".sso";
        if (-s $out){next;}  ## -s is better than -e
        $tmp=&tempname;
        open (TMP, ">$tmp");
        print TMP ">$genes\n", $genes{"$genes"}, "\n";
        close TMP;
        $sw_score=0;
        $e_val=10;
        @tmpcontent=`$fastaver -E 0.1 -H -m 10 $tmp $tdb 1`;
        open (OUT, ">$out");
        print OUT "@tmpcontent\n";
        close OUT;
        unlink ("$tmp");
        next;
    }
}


# Function: tempname
#
# Returns a unique temporary filename.
# Reasonably robust but not completely immune to race conditions
# with other processes simultaneously requesting a tempname.
#
sub tempname {
	foreach $suffix (0..99) {
		if (! (-e "tmpxx$suffix")) {
			open(TMP,">tmpxx$suffix"); # Stamp it to reserve it.
			close(TMP);
			return "tmpxx$suffix";
		}
	}
}

