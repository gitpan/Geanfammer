#!/usr/bin/perl

# Written by Sarah A. Teichmann on December 17th 1996.

# Use: creates a summary file from a cluster file.

$usage="make_summ.pl <cluster file>";


#Create summary of cluster sizes and numbers
unless (@ARGV==1){print "$usage" && die;}
my $good_cluster_file_main=shift;

&make_summ($good_cluster_file_main);

#--------Subroutine: prints out summary file from cluster file.----------


#________________________________________________________________________________
# Title     : make_clustering_summary
# Usage     : &make_summ($sorted_cluster_file)
# Function  : to make a summary file of a sorted cluster file
# Example   :
# Keywords  : summary, make_cluster_summary, subclustering summary
# Options   :
# Author    : Sarah A. Teichmann
# Date      : 19th September 1997
# Version   : 1.0
#--------------------------------------------------------------------------------
sub make_clustering_summary{
    my ($good_cluster_file, $summary_file, @filecontent, $i, $filecontent,
        $cluster_size, @cluster_sizes, $cluster_number, $number_of_clusters);
    $good_cluster_file=$_[0];
    open (FILE, "$good_cluster_file");
    $summary_file="$good_cluster_file".".summary";
    open (SUMM, ">$summary_file");
    @filecontent=<FILE>;

    for ($i=0; $i<@filecontent; $i++) {
        $filecontent=$filecontent[$i];
        if ($filecontent=~/^\s*$/){    next;        }
        if ($filecontent=~/^\s*Cluster size (\w+).*/) {
            $cluster_size=$1;
            push (@cluster_sizes, $cluster_size);
            next;
        }
        if ($filecontent=~/^Cluster\s+(\w+).*/) {
            $cluster_number=$1;
            push(@{"cluster_size_$cluster_size"},$cluster_number);
            next;
        }else { next;}
    }
    print SUMM "Cluster size          No. clusters\n";

    for ($i=0; $i<@cluster_sizes; $i++) {
        $number_of_clusters=@{"cluster_size_$cluster_sizes[$i]"};
        print SUMM "$cluster_sizes[$i]          $number_of_clusters\n";
    }

}

