#!/usr/bin/perl
# Lst Updt by /gn0/jong/Perl/update_subroutines.pl: Mon Oct 20 23:45:36 BST 1997
#______________________________________________________________________________
# Title     : convert_bla_to_msp.pl
# Usage     : convert_bla_to_msp.pl *.pbla
#             convert_bla_to_msp.pl *.bpla.gz
#             convert_bla_to_msp.pl .                <--- Give dir name!
#
# Function  : reads in PSI blast output and produces MSP file format.
#              MSP file format is a major format used a lot by people
#              in MRC-LMB genome group.
#             If the input file is gzipped with extension of 'gz', this
#              runs 'gunzip' automatically to unzip it before parsing.
# Example   : convert_bla_to_msp.pl pdb_d1svb_1.pbla.gz -b -o
# Keywords  : bla_to_msp.pl
# Options   :
#
#   $over_write=o                   by o -o ## by o -o means, by char o  or -o  at prompt
#   $pdbd_seq_only=d                by p d -d
#   $nrdb_seq_only=n                by n -n
#   $all_seq=a                      by a -a ## default option, to get all kinds of sequences.
#   $genome_seq_only=g              by g -g
#   $which_iteration=               by i=   # choose which iteration result you want to take
#   $report_only_the_best=b         by b , -b # reports the very first (best) region matched for each hit
#   $accumulative_hits_eval_thresh= by e=
#   $evalue_thresh=                 by E=   # the maximum Eval thresh in reporting any hit
#   $take_only_the_last_iteration=l by l -l
#   $verbose=v                      by v -v
#   $Accumulate_matches=A           by A -A
#   $sort_by_evalue=S               by S
#   $take_last_iter_PSI_BLA=l  by l -l
#   $BLA_entry_list_file=           by L=    ## when one column list of xx.pbla.gz files given
#
# Warning   : do not remove this 'headbox' as it is used by this perl program!
# Author    : Sarah Teichmann, Jong Park, sat@mrc-lmb.cam.ac.uk, jong@salt2.med.harvard.edu
# Version   : 2.0
#--------------------------------------------------------------------------------

exit 0 if fork;         # basic background
#use POSIX qw(setsid);
#setsid();               # disassociate from terminal etc.


push(@INC, "/Bio/App/Perl/perl5/lib") if -d "/Bio/App/Perl/perl5/lib";
my ($i, @files, $evalue_thresh, $accumulative_hits_eval_thresh);


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Important parameters
#________________________________________________________________
$evalue_thresh=30; ## default to be written to MSP file
$accumulative_hits_eval_thresh=$evalue_thresh/1000;
#$which_iteration=1;


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
# parse_arguments is a Jong's standard prompt arguments handler.
# It parses the headerbox in the above and creates global
# variables according to the options section eg) $over_write=o by o -o etc.
#___________________________________________________________________________
@files=@{&parse_arguments(1)};
@files=@{&scramble_array(\@files)};
#@files=reverse(@files);

if(@files < 2 and -d $ARGV[0]){
   print "\n Reading PWD directory or the given one \n";
   @files=@{&read_file_names_only($ARGV[0])};
}elsif(-s $BLA_entry_list_file and @ARGV < 2){
   @files=@{&open_list_file_ARRAY(\$BLA_entry_list_file)};
}elsif($ARGV[0]=~/\S\.list$/ and @ARGV < 2){
   @files=@{&open_list_file_ARRAY($ARGV[0])};
}


#$Accumulate_matches='A'; ## I am setting it as a default
#$pdbd_seq_only='d';


print "\n# $0: \$which_iteration was  $which_iteration with \"@files\" ";

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````
# This is to handle multiple xxxx.bla or xxx.pbla files,
#__________________________________________________________________
for($i=0; $i< @files; $i++){
    local ($j, @keys, $sort_by_evalue, %hash_out_final, $file, $base, $out_msp);
    $file=$files[$i];

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Take only the pbla files
    #___________________________________
    unless($file=~/\w\.[p]*bla/){      next;    }
    $base=${&get_base_names(\$file)};
    $out_msp="$base\.msp";
    print "\n# (i)  \$i=$i $0: $file output msp file  will  be \"$out_msp\"";
    print "\n# (i) \$Accumulate_matches is set to $Accumulate_matches" if $Accumulate_matches;
    if(-f $out_msp and !$over_write){
         #print "# (i) $0: $out_msp already exists, skipping \n";
         next;
    }else{
         open(MSP, ">$out_msp");
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # main function, now bla_to_msp became convert_bla_to_msp to follow the Action Oriented Programming convention
    #______________________________________________________________________________________________________________
    %hash_out_final=%{&convert_bla_to_msp(\$file,
                                          $pdbd_seq_only,
                                          "i=$which_iteration",
                                          "e=$accumulative_hits_eval_thresh",
                                          "E=$evalue_thresh",
                                          $take_only_the_last_iteration,
                                          $report_only_the_best,
                                          $genome_seq_only,
                                          $nrdb_seq_only,
                                          $all_seq,
                                          $take_last_iter_PSI_BLA,
                                          $Accumulate_matches,
                                          $verbose )};

    %hash_out=();
    %accumulative_hits=();
    $sort_by_evalue=1;
    if($sort_by_evalue){
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
         # I use sort_hash_value_by_column to sort the second column, as Eval is in the 2nd column
         #__________________________________________________________________________________________
         @keys=@{&sort_hash_value_by_column(\%hash_out_final, 2, 'n')};
    }else{
         @keys= sort {$hash_out_final{$b}<=>$hash_out_final{$a}} keys %hash_out_final;
    }

    print MSP "# (i) MSP file: $out_msp is created by $0 ", `date`, "\n";
    for($j=0; $j<@keys; $j++){
       #print "$keys[$j]\n==> $hash_out_final{$keys[$j]}\n\n";
         print MSP $hash_out_final{$keys[$j]};
    }
    close (MSP);
    if(-s $out_msp < 120){
       print "\n\a\a\a Something went wrong, $out_msp is too small\n\n\n !!!!!!!!!!\n";
       #exit
    }else{
       print "\n# (i) $0: \"$out_msp\" is created \n\n";
    }
}

&show_options;

print "\n# (i) \$Accumulate_matches was set to $Accumulate_matches\n" if $Accumulate_matches;



#________________________________________________________________________________
# Title     : convert_bla_to_msp
# Usage     : %hash_out_final=%{&convert_bla_to_msp(\$file, [$Lean_output])};
# Function  : reads in PSI blast output and produces MSP file format.
#             Takes all the good hits below certain threshold in multiple iteration
#             Reports the best evalue with a given sequence name
# Example   : %hash_out=%{&convert_bla_to_msp(\$file)};
# Keywords  : pbla_to_msp, blast_to_msp, bla_2_msp, blastp_to_msp_format,
#             blast_to_msp_format, convert_bla_to_msp, convert_bla_to_msp_files
#             bla_to_msp
# Options   :
#   $pdbd_seq_only  d   for getting dxxxx_ like seq names only(pdb40d names for examp)
#   $all_seq  a         for forcing all seq conversion
#   $which_iteration= by i=    # choose which iteration result you want to take
#   $which_iteration   as just a digit
#   $report_only_the_best=b by b -b
#   $take_last_iter_PSI_BLA=l by l
#   $PSI_BLA_ACCUMU_hits_eval_thresh= by e=
#   $genome_seq_only=g      by g
#   $nrdb_seq_only=n        by n
#   $evalue_thresh=         by E=
#   $Accumulate_matches=A   by A -A
#   $Lean_output=L          by L -L  # to remove search output to unclutter
#
# Author    : Sarah Teichmann and Jong Park, jong@salt2.med.harvard.edu
# Version   : 5.2
#--------------------------------------------------------------------------------
sub convert_bla_to_msp{
    my($i, $j, $k, @lines, $match_string_count,  $line_count, $query_string_count,
       $match_length, $Lean_output, $SEQ_NAME, $original_query,
       $duplicated_match_count, $new_sorted_name, $sorted_name, $verbose,
       $pdbd_seq_only, $which_iteration, $report_only_the_best,
       $genome_seq_only, $all_seq, $header_found, $PSI_BLA_ACCUMU_hits_eval_thresh,
       $take_last_iter_PSI_BLA, $nrdb_seq_only, $system_mem_size,
       $get_the_final_iteration, $read_entry_lines, $verbose, $Accumulate_matches,
       $CONVERGED_sign_found, $Evalue_limit, $entry_and_alignment_found, $query);
    my $match_leng_thresh=10;
    ### This localization is critial NOT my, as I use a sub which relies on this
    local(%hash_out, %accumulative_hits, $file, $score, $score_ori, $evalue,
          $evalue_ori, $seq_id, $query_range_start, $query_range_stop,
          $match_string_start, $match_string_stop, $matched, $matched_seq_name,
          $read_point_found, $summary_lines_found, $entry_found, %good_matches_list);
    $duplicated_match_count=0;

    $Evalue_limit=5;
    $PSI_BLA_ACCUMU_hits_eval_thresh=0.0001; ## default eval threshes
    $query=$original_query='query_seq'; ## default query seq name, to avoid blank name
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Processing the input arguments to get file and options etc
    #_____________________________________________________________
    for (@_){
        if(ref $_ eq 'ARRAY'){ @lines =@{$_};
        }elsif( ref $_ eq 'SCALAR' and -s ${$_} ){ $file=${$_};
        }elsif( -s $_ ){            $file=$_;
        }elsif(/^ *d *$/){          $pdbd_seq_only='d'; $all_seq=''; $genome_seq_only='';
            print "\n $0: convert_bla_to_msp,  You set \$pdbd_seq_only option, I will skip others.\n";
        }elsif(/^ *b *$/){          $report_only_the_best='b';
        }elsif(/^ *a *$/){          $all_seq='a'; $genome_seq_only=''; $pdbd_seq_only=''; $nrdb_seq_only='';
        }elsif(/^ *g *$/){          $genome_seq_only='g'; $all_seq=''; $pdbd_seq_only='';$nrdb_seq_only='';
        }elsif(/^ *n *$/){          $nrdb_seq_only='n'; $all_seq=''; $pdbd_seq_only=''; $genome_seq_only=''; }
        if(/^ *l *$/){          $take_last_iter_PSI_BLA='l'; $Accumulate_matches='' }
        if(/^ *v *$/){          $verbose='v'; }
        if(/^ *L *$/){          $Lean_output='L'; }
        if(/^\s*e=(\S+)/){          $PSI_BLA_ACCUMU_hits_eval_thresh=$1; }
        if(/^\s*SEQ_NAME=(\S+)/i){  $query=$original_query=$SEQ_NAME=$1;  }
        if(/^\s*E=(\S+)\s*/){       $Evalue_limit=$1;          }
        if(/^\s*A$/){           $Accumulate_matches='A'; $take_last_iter_PSI_BLA=''; }
        if(/^ *i= *(\d+) *$/){ $which_iteration=$1; }
     }
     unless($which_iteration){  $get_the_final_iteration=1 }

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
     # If the input file is gzipped, uncompress it to text file and then open
     #__________________________________________________________________
     if($file=~/\.gz *$/){
         open(BLA_FILE, "gunzip -c $file|") || die "\n# $0: Failed to open $file\n\n\n";
         if($file=~/^([de]*\d\d*\w\w\w\w\w)\./){         $query=$1;          }
     }else{
         open(BLA_FILE, "$file") || die "\n# !! $0: convert_bla_to_msp : Failed to open \"$file\"\n\n\n";
         if($file=~/^([de]*\d\d*\w\w\w\w\w)\./){         $query=$1;          }
     }
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # UP to NOW is frivalous option handling stuff
     #_______________________________________________________
     $system_mem_size=${&get_total_memory_size_in_linux};
     if($system_mem_size < 128000000){
         ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
         ###                                                                           ###
         ###  (1) Main reading in .pbla file (or any extension)                        ###
         ###   by putting pattern matches which occur most, I can save comparisons     ###
         ###___________________________________________________________________________###
         BLA1: while(<BLA_FILE>){

             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
             # This is just to skip garbage lines
             #_______________________________________________________
             if(/^\s*$/ or /^\s\s\s+Length +\= +\d+ *$/){      next     }  ## skipping some junk lines

             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (1) The most important READING HERE!
             #___________________________________________
             if($entry_and_alignment_found or ($summary_lines_found and /^\>\S/)){
                  ($matched, $entry_and_alignment_found)
                      =&match_seq_entry_and_alignment_block_in_BLAST_output($original_query);
                  $summary_lines_found=0;
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (3) '> xxxx '  New sequence entry, '>' starts
             #__________________________________________________________
             }elsif(/^(\S+)\s.{29,}\s\d+\s+(\S{3,8})\s*$/){   ## mind the size of space!!
                  $name_matched=$1;
                  $evalue=$2;
                  &match_summary_head_lines_in_BLAST_output($name_matched, $evalue);
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~
             # (4) 'Searching......done'  line indicates new search step(iteration)
             #_________________________________________________________________________
             }elsif( /^\s{0,4}Searching\.\.\.+[done]?/i ){
                 $which_iteration=&match_Searching_dot_line_in_BLAST_output($which_iteration);
                 $summary_lines_found=1;
                 next BLA1;
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (5) Check if it converged before the given -j value
             #________________________________________________________
             }elsif(/^\s*CONVERGED/){
                 $CONVERGED_sign_found=1;
                 $which_iteration=$read_point_found;
                 $entry_and_alignment_found=1;
                 next BLA1;
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (6) Extracting query seq name(this is the only place to get it)
             #__________________________________________________________________
             }elsif(/^\s{0,4}Query=\s+(\S+)/){
                  $query=$original_query=$1;
                  next BLA1;
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (6-1) Following is to handle the HTML version of PSI output
             #___________________________________________________________________
             }elsif(/\<\S\> *Query=\<\S\>/i){ $query=$original_query=$SEQ_NAME;
                 next BLA1;  # <b>Query=</b>
             }elsif(eof){
                  @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
                        \%accumulative_hits, $SEQ_NAME,$matched,$evalue, $score, $seq_id,
                        $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                        $match_string_stop, $read_point_found, $accumulative_hits_eval_thresh,
                        $take_last_iter_PSI_BLA, $accumulative_hits_eval_thresh, $Evalue_limit)};
                  %hash_out=         %{$out_from_put_msp_lines[0]};
                  %accumulative_hits=%{$out_from_put_msp_lines[1]};
                  $read_point_found= $out_from_put_msp_lines[2];
                  last;
             }elsif(/^ +\*+ +No hits found +\*+/i){
                  print "\n $_ \n";
                  last;
             }
         }
         close(BLA_FILE);
     }else{
         ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
         ###                                                                           ###
         ###  (1) Main reading in .pbla file (or any extension)                        ###
         ###   by putting pattern matches which occur most, I can save comparisons     ###
         ###___________________________________________________________________________###
         @lines=<BLA_FILE>;
         BLA2: for (@lines){
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
             # This is just to skip garbage lines
             #_______________________________________________________
             if(/^\s*$/ or /^\s\s\s+Length +\= +\d+ *$/){  next BLA2;  }  ## skipping some junk lines

             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (1) The most important READING HERE!
             #___________________________________________
             if($entry_and_alignment_found or ($summary_lines_found and /^\>\S/)){
                  ($matched, $entry_and_alignment_found)
                      =&match_seq_entry_and_alignment_block_in_BLAST_output($original_query);
                  $summary_lines_found=0;
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (3) '> xxxx '  New sequence entry, '>' starts
             #__________________________________________________________
             }elsif(/^(\S+)\s.{29,}\s\d+\s+(\S{3,8})\s*$/){   ## mind the size of space!!
                  $name_matched=$1;
                  $evalue=$2;
                  &match_summary_head_lines_in_BLAST_output($name_matched, $evalue);
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~
             # (4) 'Searching......done'  line indicates new search step(iteration)
             #_________________________________________________________________________
             }elsif( /^\s{0,4}Searching\.\.\.+[done]?/i ){
                 $which_iteration=&match_Searching_dot_line_in_BLAST_output($which_iteration);
                 $summary_lines_found=1;
                 next BLA2;
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (5) Check if it converged before the given -j value
             #________________________________________________________
             }elsif(/^\s*CONVERGED/){
                 $CONVERGED_sign_found=1;
                 $which_iteration=$read_point_found;
                 $entry_and_alignment_found=1;
                 next BLA2;
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (6) Extracting query seq name(this is the only place to get it)
             #__________________________________________________________________
             }elsif(/^\s{0,4}Query=\s+(\S+)/){
                  $query=$original_query=$1;
                  next BLA2;
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # (6-1) Following is to handle the HTML version of PSI output
             #___________________________________________________________________
             }elsif(/\<\S\> *Query=\<\S\>/i){ $query=$original_query=$SEQ_NAME;
                 next BLA2;  # <b>Query=</b>
             }elsif(eof){
                  @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
                        \%accumulative_hits, $SEQ_NAME,$matched,$evalue, $score, $seq_id,
                        $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                        $match_string_stop, $read_point_found, $accumulative_hits_eval_thresh,
                        $take_last_iter_PSI_BLA, $accumulative_hits_eval_thresh, $Evalue_limit)};
                  %hash_out=         %{$out_from_put_msp_lines[0]};
                  %accumulative_hits=%{$out_from_put_msp_lines[1]};
                  $read_point_found= $out_from_put_msp_lines[2];
                  last;
             }elsif(/^ +\*+ +No hits found +\*+/i){
                  print "\n $_ \n";
                  last;
             }
         }
         close(BLA_FILE);
     }

     unless( $take_last_iter_PSI_BLA){
         print "\n# >> ACCUMULATIVE HITS are reported as you did not set \$take_last_iter_PSI_BLA opt!!\n";
         %hash_out=(%hash_out, %accumulative_hits);
     }

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```
     # CLeaning up the BLA file if $Lean_output is set
     #_____________________________________________________
     $gzipped_search_file="$file\.gz";
     if($Lean_output ){ ## If Lean_out opt is set and $file exists and %hash_out is not empty, remove $file
          if(-s $file){
                  unlink($file);  ## removes fam_8_8.pbla etc,
          }elsif(-s $gzipped_search_file){
                  unlink($gzipped_search_file); ## removes fam_8_8.pbla.gz etc,
          }else{
             print "\n# (E) convert_bla_to_msp: tried to remove search out file for \$Lean_output opt,
                   but failed. Something is wrong. Think! or report to jong\@salt2.med.harvard.edu,
                   jong\@mrc-lmb.cam.ac.uk, sat\@mrc-lmb.cam.ac.uk, jong_p\@hotmail.com\n";
                   exit;
          }
     }
     return(\%hash_out);

            sub match_summary_head_lines_in_BLAST_output{
                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
                 # reading the search summary lines. Save time by selecting which match to parse
                 # "d2sga__ 2.35.1.1.3 Protease A [(Streptomyces griseus), strain k1]    273  7e-74"
                 #_________________________________________________________________________________________
                 $matched_seq_name=$_[0];
                 $evalue=$_[1];

                 if($matched_seq_name=~/pdb\|(\S+)\|(\S+)$/i){ $matched_seq_name="$1$2"
                 }elsif($matched_seq_name=~/^gi\|\S*?\|?([^\|]+)$/i
                    or  $matched_seq_name=~/^\S+\|\S*\|([^\|]+)$/){ $matched_seq_name=$1
                 }
                 if($evalue <= $Evalue_limit){
                     if($pdbd_seq_only and ($matched_seq_name=~/^pdb_/
                        or $matched_seq_name=~/^[cde]\d\w{3,6}/)
                        or $matched_seq_name=~/^ds[\d\_]+$/){
                         $good_matches_list{$matched_seq_name}=$matched_seq_name;
                     }elsif(!$pdbd_seq_only){
                         $good_matches_list{$matched_seq_name}=$matched_seq_name;
                     }
                 }
                 return($matched_seq_name, $evalue);
            }

            sub match_seq_entry_and_alignment_block_in_BLAST_output{
                $entry_and_alignment_found=1;
                $original_query=$_[0];

                if(/^\> {0,4}(\S+)/){
                     $temp_match=$1;
                     if($temp_match=~/pdb\|(\S+)\|(\S+)$/i){ $temp_match="$1$2"
                     }elsif($temp_match=~/^gi\|\S*?\|?([^\|]+)$/i
                        or  $temp_match=~/^\S+\|\S*\|([^\|]+)$/){ $temp_match=$1 }

                     unless($good_matches_list{$temp_match}){  next }
                     if($match_string_count){ ## $match_string_count is incremented only by 'Sbjct' line
                           @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
                                  \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                  $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                  $match_string_stop, $read_point_found, $PSI_BLA_ACCUMU_hits_eval_thresh,
                                  $take_last_iter_PSI_BLA, $PSI_BLA_ACCUMU_hits_eval_thresh, $Evalue_limit)};
                           %hash_out=         %{$out_from_put_msp_lines[0]};
                           %accumulative_hits=%{$out_from_put_msp_lines[1]};
                           $read_point_found= $out_from_put_msp_lines[2];
                           $match_string_count=0;
                           $duplicated_match_count=0;
                     }

                     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     #  Only with new seq entry, I count the pair occurrances
                     #__________________________________________________________________
                     $query=$original_query;
                     $query_string_count='';
                     $matched=$temp_match; ## this should be here, after if
                     $sorted_name=join(' ', sort($query, $matched) );
                     $match_string_count=0;
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # (1) Matching  Score =  325 bits (824), Expect = 6e-89           << 2 >>
                #________________________________________________________________________________
                elsif( /^\s*Score\s*\=\s*(\S+)\s*bits +\(\S+\)\,\s*Expect\s*=\s*(\S+)/i
                    or /^\s*Score\s*\=\s*(\S+)\s*bits.+\,\s*Expect\s*=\s*(\S+)/i){
                    $score_ori=$1;  $evalue_ori=$2;
                    if($evalue_ori=~/^e\-\d\d\d/){ $evalue_ori="1".$evalue_ori; } ## bug fix for short eval in blast distribution

                    if($match_string_count){ # $match_string_count is increased when Sbjct word is found
                         if($evalue > $Evalue_limit){ $evalue=$evalue_ori; $score=$score_ori; }
                         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
                         # When Only the first match(best evalue) is required, write msp line and reset $entry_found var
                         #_________________________________________________________________________________________________
                         if($report_only_the_best){
                              #print "      (5)  \$report_only_the_best is set\n" if $verbose;
                              @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
                                      \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                      $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                      $match_string_stop, $read_point_found, $PSI_BLA_ACCUMU_hits_eval_thresh,
                                      $take_last_iter_PSI_BLA, $PSI_BLA_ACCUMU_hits_eval_thresh, $Evalue_limit)};
                              %hash_out=         %{$out_from_put_msp_lines[0]};
                              %accumulative_hits=%{$out_from_put_msp_lines[1]};
                              $read_point_found= $out_from_put_msp_lines[2];
                              #$entry_found=0;
                         }else{
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~`
                              # duplicated match count means, query matched more than one region of a match seq
                              #__________________________________________________________________________________
                              $duplicated_match_count++;
                              $sorted_name="$sorted_name $duplicated_match_count";
                              @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
                                      \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                      $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                      $match_string_stop, $read_point_found, $PSI_BLA_ACCUMU_hits_eval_thresh,
                                      $take_last_iter_PSI_BLA, $PSI_BLA_ACCUMU_hits_eval_thresh, $Evalue_limit)};
                              %hash_out=         %{$out_from_put_msp_lines[0]};
                              %accumulative_hits=%{$out_from_put_msp_lines[1]};
                              $read_point_found= $out_from_put_msp_lines[2];
                         }
                         $score=$score_ori; $evalue=$evalue_ori;
                    }else{
                         $evalue=$evalue_ori; $score=$score_ori;
                    } ## to next line

                    sub reset_all_the_vars{
                         #print "            !!!!  Reseting all the vars !!!!\n" if $verbose;
                         $query_string_count=$score=$evalue=$seq_id=$query_range_stop=$query_range_start='';
                         $match_string_stop=$msp_line=$new_sorted_name=$match_string_count='';
                     }
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
                # (2) Matching   Identities = 158/158 (100%), Positives = 158/158 (100%)    ,
                #____________________________________________________________________________________
                elsif( /^ {0,4}Identities += +\S+\/(\S+) +\( *(\S+) *\%\)/i){
                     $query_string_count=$match_string_count=0;
                     $seq_id=$2/100;
                     $match_length=$1;
                     if($match_length < $match_leng_thresh){  $match_string_count=1; }
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # (3) Matching 'Query: 2 GIRAATSQEINELT..' line    ,
                #_________________________________________________________________
                elsif(/^ {0,4}Query\:?\s+(\d+) +\D+ +(\d+)/){
                     $query_string_count++;
                     $query_line_found=1;
                     if($query_string_count==1){      $query_range_start=$1;   $query_range_stop =$2;
                     }elsif($query_string_count > 1){ $query_range_stop=$2;     }
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # (4) Matching 'Sbjct: 2 GIRAATSQEINELT..' line
                #_________________________________________________________________
                elsif($query_line_found and /^ {0,4}Sbjct\:? +(\d+) +[\w\-]+ +(\d+)/i){
                     $match_string_count++;
                     $subject_line_found=1;
                     if($match_string_count==1){      $match_string_start=$1;
                                                      $match_string_stop =$2;
                     }elsif($match_string_count > 1){ $match_string_stop=$2;      }
                }
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # (5) Matching '   Database: ' line    ,                << END >>
                #_________________________________________________________________
                elsif(/^\s+Database:\s+\S+/ or eof){ # the very last write
                     if($evalue > $Evalue_limit){
                     }else{
                          @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
                                \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                $match_string_stop, $read_point_found, $PSI_BLA_ACCUMU_hits_eval_thresh,
                                $take_last_iter_PSI_BLA, $PSI_BLA_ACCUMU_hits_eval_thresh, $Evalue_limit)};
                          %hash_out=         %{$out_from_put_msp_lines[0]};
                          %accumulative_hits=%{$out_from_put_msp_lines[1]};
                          $read_point_found= $out_from_put_msp_lines[2];
                     }
                }
                return($matched, $entry_and_alignment_found);
             }

             sub match_Searching_dot_line_in_BLAST_output{
                  $which_iteration=$_[0];
                  %good_matches_list=();
                  $read_point_found++;
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
                  #  (3.3) Following is the KEY part for controlling iteration
                  #__________________________________________________________
                  if( $read_point_found < $which_iteration){
                       $match_string_count=$query_string_count=$score=$evalue=$seq_id=$score_ori=$evalue_ori='';
                       $query_range_stop=$query_range_start=$match_string_stop=$msp_line=$new_sorted_name='';
                       $duplicated_match_count=0;
                       if( !$Accumulate_matches){  %hash_out=(); } ## this is to remove any discarded pairs in the iteration
                  }elsif( $read_point_found == $which_iteration){

                  }elsif( $which_iteration and $read_point_found >  $which_iteration){
                       @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
                                  \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                  $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                  $match_string_stop, $read_point_found, $PSI_BLA_ACCUMU_hits_eval_thresh,
                                  $take_last_iter_PSI_BLA, $PSI_BLA_ACCUMU_hits_eval_thresh, $Evalue_limit)};
                       %hash_out=         %{$out_from_put_msp_lines[0]};
                       %accumulative_hits=%{$out_from_put_msp_lines[1]};
                       $read_point_found= $out_from_put_msp_lines[2];
                       last;
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  # If you did not set the which iteration option
                  #_________________________________________________________
                  }elsif(!$which_iteration){
                       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
                       # (3.4) Default situation
                       #____________________________________________________________
                       print "\n# (WARN) You did not set \$which_iteration option \n\n" if $verbose;
                       if($read_point_found > 1){
                             @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
                                    \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                    $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                    $match_string_stop, $read_point_found, $PSI_BLA_ACCUMU_hits_eval_thresh,
                                    $take_last_iter_PSI_BLA, $PSI_BLA_ACCUMU_hits_eval_thresh, $Evalue_limit)};
                             %hash_out=         %{$out_from_put_msp_lines[0]};
                             %accumulative_hits=%{$out_from_put_msp_lines[1]};
                             $read_point_found= $out_from_put_msp_lines[2];
                       }
                       $match_string_count=$query_string_count=$score=$evalue=$seq_id=$score_ori=$evalue_ori='';
                       $query_range_stop=$query_range_start=$match_string_stop=$msp_line=$new_sorted_name='';
                       $entry_found=$duplicated_match_count=0;
                       if( !$Accumulate_matches){  %hash_out=(); $entry_found=0; $duplicated_match_count=0;     }
                  }
                  return($which_iteration);
             }


}





#______________________________________________________________________________
# Title     : put_msp_lines_to_hash_from_bla
# Usage     : @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
#                                        $query,$matched,$evalue, $score, $seq_id,
#                                        $sorted_name, $query_range_start,
#                                        $query_range_stop,$match_string_start,
#                                        $match_string_stop, $read_point_found,
#                                        $accumulative_hits_eval_thresh,
#                                        $take_last_iter_PSI_BLA)};
# Function  :
# Example   :
# Keywords  :
# Options   :
# Author    : jong@salt2.med.harvard.edu,
# Category  :
# Version   : 1.4
#------------------------------------------------------------------------------
sub put_msp_lines_to_hash_from_bla{
    my (@finale_out, $sorted_name, $msp_line, $evalue, $score, $matched,
        $seq_id, $query_range_start,$accumulative_hits_eval_thresh,
        $query_range_stop, $query, $match_string_start, $match_string_stop,
        $read_point_found, %hash_out, %accumulative_hits, $Evalue_thresh);
    %hash_out=%{$_[0]};         %accumulative_hits=%{$_[1]};
    $query=$_[2];               $matched=$_[3];
    $evalue=$_[4];              $score=$_[5];
    $seq_id=$_[6];              $sorted_name=$_[7];
    $query_range_start=$_[8];   $query_range_stop =$_[9];
    $match_string_start=$_[10]; $match_string_stop=$_[11];
    $read_point_found=$_[12];   $accumulative_hits_eval_thresh=$_[13];
    $take_last_iter_PSI_BLA=$_[14];
    $accumulative_hits_eval_thresh=$_[15];
    $Evalue_thresh=$_[16];
    $query  ="$query\_$query_range_start\-$query_range_stop";

    if($matched !~/^\S+\_\d+\-\d+ *$/){         $matched="$matched\_$match_string_start\-$match_string_stop";
    }elsif($matched =~/^(\S+)\_\d+\-\d+ *$/){   $matched="$1\_$match_string_start\-$match_string_stop";     }

    if($score=~/\S/ and $evalue=~/\S/ and $match_string_start=~/\S/ and $Evalue_thresh > $evalue){
        $msp_line=sprintf("%-6s %-9s %-5s %-5s %-5s %-32s %-s\t%-s\t%-s\t%-s\n",
                           $score, $evalue, $seq_id, $query_range_start, $query_range_stop,
                           $query, $match_string_start, $match_string_stop, $matched, $read_point_found);
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # This is where I really put the matches !!!
        #_____________________________________________________
        if($hash_out{$sorted_name}=~/^\S+ +(\S+) +/){
            if($1 >= $evalue){
                print "                    (1) put_msp_lines_to_hash_from_bla: $1 >= $evalue WRITING to hash. 1\n" if $verbose;
                $hash_out{$sorted_name}=$msp_line;
            }else{
                print "                    put_msp_lines_to_hash_from_bla: $1 < $evalue_ori NO write to hash\n" if $verbose;  }
            }else{
                print "                    (2) put_msp_lines_to_hash_from_bla: NO eval >= $evalue WRITING to hash. 2\n" if $verbose;
                $hash_out{$sorted_name}=$msp_line;
            }
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # This part is to rescue the hits dropped by matrix migration
            #_________________________________________________________________
            if(!$take_last_iter_PSI_BLA and $evalue <= $accumulative_hits_eval_thresh ){
                if($accumulative_hits{$sorted_name}){
                    if($accumulative_hits{$sorted_name}=~/^[\t ]*\S+[\t ]+(\S+)[\t ]/){
                         if($evalue < $1){
                                 $accumulative_hits{$sorted_name}=$msp_line;   }   }
                }else{ $accumulative_hits{$sorted_name}=$msp_line;     }
            }
    }else{
    }
    @finale_out=(\%hash_out, \%accumulative_hits, $read_point_found, $query, $matched, $evalue, $score, $seq_id, $sorted_name,
                 $query_range_start, $query_range_stop, $match_string_start, $match_string_stop  );
    return(\@finale_out);
}






#________________________________________________________________________
# Title     : assign_options_to_variables
# Usage     : &assign_options_to_variables(\$input_line);
# Function  : Assigns the values set in head box to the variables used in
#             the programs according to the values given at prompt.
#             This produces global values.
#             When numbers are given at prompt, they go to @num_opt
#              global variable. %vars global option will be made
#
# Example   : When you want to set 'a' char to a variable called '$dummy' in
#             the program, you put a head box commented line
#             '#  $dummy    becomes  a  by  -a '
#             Then, the parse_arguments and this sub routine will read the head
#             box and assigns 'a' to $dummy IF you put an argument of '-a' in
#             the prompt.
# Warning   : This is a global vars generator!!!
# Keywords  :
# Options   : '#' at prompt will make a var  $debug set to 1. This is to
#              print out all the print lines to make debugging easier.
# Returns   : Some globaly used variables according to prompt options.
#             @num_opt,
#
# Argument  : None.
# Version   : 2.6
#--------------------------------------------------------------------
sub assign_options_to_variables{
	my($i, %vars, $j, $op, $z, $n, $symb, $value, $var, %val, @val, $ARG_REG,
	 $option_table_example, @input_options, $first_border_and_title, $sym, @arg);

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#      Defining small variables for option table reading
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my($g)='gets';                my($if)='if';
	my($is)='is';                 my(@input_files);
	my($o)='or';   my(@arguments) = sort @ARGV;

	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#  Assigning global arguments(@num_opt, %vars) variables
	#_______________________________________________________________
	for($i=0; $i< @arguments; $i++){
	 if(($arguments[$i]=~/^(\-?\d+[\.\d+]?)$/)&&   ### it mustn't be a file
		( !(-f $arguments[$i]) ) ){                ### getting NUM opt
		push(@num_opt, $1);
	 }elsif( $arguments[$i]=~/^(\S+)=(\S+)$/){
		$vars{$1}=$2;
	 }
	}

	#""""""""""""""""""""""""""""""""""""""""""""""""""
	#   Some DEFAULT $debug variables for debugging purposes
	#""""""""""""""""""""""""""""""""""""""""""""""""""
	&set_debug_option;

	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#   The main processing of self
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	open(SELF, "$0");    ## opens the program you ran to get the options table.
	while(<SELF>){

	  if( $first_border_and_title > 6 ){  ## This is to make it read only the first headbox.
		  last;                            #  $first_border_and_title is an incremental counter.
	  }elsif( /^ *#[_\*\-]{15,}$/ and /^ *# *[Tt][itle]*[ :]*/ ){
		  $first_border_and_title++;
		  print __LINE__, "# assign_options_to_variables : Title line found\n" if $debug eq 1;
	  }elsif(/^ {0,5}# {1,50}[\$\%\@].+$/){
		  $op = $&;  ## $op is for the whole input option line which has $xxxx, @xxx, %xxxx format
		  $op =~ s/^( *\# *)(\W\w+.+)$/$2/;  ## This is removing '#  ' in the line.
		  $op =~ s/^(\W\w+.+)(\s+\#.*)$/$1/;  ## This is removing any comments in the line.
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 ## matching the following line input format.
			 ## $av_sc_segment     becomes    a  by  a  # To smooth the SC rates. Gets the averages of
			 ## $ARG_REG is for arguments regular expression variable.
			 ##  This reg. exp. matches = 'a or A or E or e' part
			 ##  which represents alternative prompt arguments possibilities. \=$b$g$is$e$set
			 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			 $ARG_REG ='(\S*) *[or=\,]* *(\S*) *[=or\,]* *(\S*) *[=or\,]* *(\S*) *[=or\,]* *(\S*)';
			 if($op=~/^([\$\@\%])([\w\-]+) {0,20}[=|$g|$is] *[\$\@\%]*([\- \w\.\d]+) *[bB]y +$ARG_REG/){
							 ## $sym     $var        becomes          a [$a...]       by       a -a -A
				  my $sym = $1;  #### The symbols like ($, @, %), '$' in the above.
				  my $var = $2;  #### Actual variable name 'var' from $var, 'av_sc_segment' in the above.
				  my $val = $3;  #### The becoming value  first 'a' in the above.
				  my @arg = ($4, $5, $6, $7, $8);  ## The alternative prompt arguments, second 'a' in the above..
			      print "\n $sym $var $val \n" if $debug==1;
			      print "\n \@arg are @arg \n" if $debug==1;

				  #""""""""""""""""""""""""""""""""""""""""""""""""""""
				  #  Going through the PROMPT args.
				  #""""""""""""""""""""""""""""""""""""""""""""""""""""
				  for($z=0; $z < @arguments; $z++){     ## $arguments[$z]  is from @ARGV
					  if($arguments[$z]=~/^\-\w+$/){
						  $arguments[$z] =~ s/\-//;
					  }
					  for ($i=0; $i < @arg; $i ++ ){
						 if( ("$arg[$i]" eq "$arguments[$z]" )&& ($arg[$i] !~ /\=/)
							 && ($sym eq '$') ){
							 ${"$var"}="$val";
							 if($debug == 1){
								 print __LINE__," \$${var} is set to \"$1\"\n";
							 }

						 }#'''''''''''''''' $arg = by s=  syntax ~~~~~~~~~~~~~~~~~~~~~~~~~~~
						 elsif( ( $arg[$i] =~ /^(\w+) *\=/ ) &&
							( $arguments[$z] =~ /^${1}= *([\w\.*\-*]+)$/) &&
							( $sym eq '$') ){
							  ${"$var"}="$1";
							  if($debug eq 1){ print __LINE__,"\$${var} is set to \"$1\"\n";  }
						 }
					  }
				  }
			  }
		}
	}
}

#________________________________________________________________________
# Title     : set_debug_option
# Usage     : &set_debug_option;
# Function  : If you put '#' or  '##' at the prompt of any program which uses
#             this sub you will get verbose printouts for the program if the program
#             has a lot of comments.
# Example   : set_debug_option #    <-- at prompt.
# Warning   :
# Keywords  :
# Options   : #   for 1st level of verbose printouts
#             ##  for even more verbose printouts
# $debug  becomes 1 by '#'  or '_'
# $debug2 becomes 1 by '##'  or '__'
#
# Returns   :  $debug
# Argument  :
# Version   : 1.8
#--------------------------------------------------------------------
sub set_debug_option{
	my($j, $i, $level);
	unless( defined($debug) ){
	 for($j=0; $j < @ARGV; $j ++){
		 if( $ARGV[$j] =~/^(_+)$|^(#+)$/){ # in bash, '#' is a special var, so use '_'
			 print __LINE__," >>>>>>> Debug option is set by $1 <<<<<<<<<\n";
			 $debug=1;
				  print chr(7);
			 print __LINE__," \$debug  is set to ", $debug, "\n";
			 splice(@ARGV,$j,1); $j-- ;
			 $level = length($1)+1;
			 for($i=0; $i < $level; $i++){
				 ${"debug$i"}=1;
				 print __LINE__," \$debug${i} is set to ", ${"debug$i"}, "\n";
			 }
		 }
	 }
	}
}
#________________________________________________________________________
# Title     : default_help
# Usage     : &default_help2;  usually with 'parse_arguments' sub.
# Function  : Prints usage information and others when invoked. You need to have
#             sections like this explanation box in your perl code. When invoked,
#             default_help routine reads the running perl code (SELF READING) and
#             displays what you have typed in this box.
#             After one entry names like # Function :, the following lines without
#             entry name (like this very line) are attached to the previous entry.
#             In this example, to # Function : entry.
# Example   : &default_help2; &default_help2(\$arg_num_limit);   &default_help2( '3' );
#             1 scalar digit for the minimum number of arg (optional),
#             or its ref. If this defined, it will produce exit the program
#             telling the minimum arguments.
# Warning   : this uses format and references
# Keywords  :
# Options   :
# Returns   : formated information
# Argument  :
# Version   : 3.3
#--------------------------------------------------------------------
sub default_help{
	my($i, $perl_dir, $arg_num_limit, $head ,$arg_num_limit, $key_press, $e,
	  @entries, @entries_I_want_write, $option_tb_found, $extension, $logname, $tmp );
	$logname=getlogin();
	my($pwd)=`pwd`;
	my($date)=`date`;
	chomp($date,$pwd);
	my($not_provided)="--- not provided ---\n";
	my($file_to_read) = $0;

	for($i=0; $i < @_; $i ++){
	  if((ref($_[$i]) eq 'SCALAR')&&(${$_[$i]} =~ /^\d$/)){
		  $arg_num_limit = ${$_[$i]};  }
	  elsif( (!(ref($_[$i]))) && ($_[$i] =~ /^\d$/)){
		  $arg_num_limit = $_[$i];     }
	}
	my %entries = %{&read_head_box(\$file_to_read )};
	if($option_tb_found ==1){
	 @option_tb=@{&read_option_table(\$file_to_read)};
	}

	@entries = keys %entries;
	foreach $help_item (@entries){
	  ${$help_item}= $not_provided if( ${$help_item}=~/^[\W]*$/  and  !defined(${$help_item}) );
	}
	#""""""""""""""""""""""""""""""""""""""""
	#########  Writing the format <<<<<<<<<<<
	#""""""""""""""""""""""""""""""""""""""""
	$~ =HEADER_HELP;
	write;   ## <<--  $~ is the selection operator
	$~ =DEFAULT_HELP_FORM;

	@entries_I_want_write=sort keys %entries;

	for( @entries_I_want_write ){  write  }

	print chr(7);  print "_"x72,"\n\n";

	if(@ARGV < $arg_num_limit){ print "\* $0 fataly needs $arg_num_limit arguments\n\n" }

	if(  $option_tb_found == 1){
	 #########  Printing the OPTION table contents <<<<<<<<<<<<
	 print "  Press \"Return\" key to see what options $logname\'s \n\n    \"$0\" take... \n";
				 $key_press=getc();
	 print @option_tb, "\n"x2 if(@option_tb > 0);
	}
format HEADER_HELP  =
_____________________________________________________________________
		  __  __      ______     __          _____
		 /\ \/\ \    /\  ___\   /\ \        /\  _ `\
		 \ \ \_\ \   \ \ \__/   \ \ \       \ \ \L\ \
		  \ \  _  \   \ \  _\    \ \ \       \ \ ,__/
		   \ \ \ \ \   \ \ \/___  \ \ \_____  \ \ \/
		    \ \_\ \_\   \ \_____\  \ \______\  \ \_\
		     \/_/\/_/    \/_____/   \/______/   \/_/ V 3.1`
_____________________________________________________________________
.
format DEFAULT_HELP_FORM =
 @<<<<<<<<<: @*
 $_,        $entries{$_}
.
}


#________________________________________________________________________
# Title     : parse_arguments
# Usage     : &parse_arguments; or  (file1, file2)=@{&parse_arguments};
# Function  : Parse and assign any types of arguments on prompt in UNIX to
#             the various variables inside of the running program.
#             This is more visual than getopt and easier.
#             just change the option table_example below for your own variable
#             setttings. This program reads itself and parse the arguments
#             according to the setting you made in this subroutine or
#             option table in anywhere in the program.
#             It also imports the ENV variables to your program.
#
# Example   : &parse_arguments(1);
#             @files=@{&parse_arguments(1)};
# Warning   : HASH and ARRAY mustn't be like = (1, 2,3) or (1,2 ,3)
# Keywords  :
# Options   : '0'  to specify that there is no argument to sub, use
#              &parse_arguments(0);
#             parse_arguments itself does not have any specific option.
#             '#' at prompt will make a var  $debug set to 1. This is to
#              print out all the print lines to make debugging easier.
#             'e=xxxx' for filtering input files by extension xxxx
#
# Returns   : Filenames in a reference of array
#             and input files in an array (file1, file2)=@{&parse_arguments};
# Argument  : uses @ARGV
# Version   : 1.9
#--------------------------------------------------------------------
sub parse_arguments{
	my( $c, $d, $e, $extension, $f, $arg_num, $option_table_seen, $n, $option_filtered,
		$option_table_example, $tmp, $logname, $input_line, @input_files,
		@extension);

	&import_ENV_vars;

	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#   Checks if there were arguments
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	if( @ARGV < 1 ){ #<-- If Argument is not given at prompt
	  for(@_){
		 if($_ eq '0'){
			 last;
		 }else{
			 print "\n \"$0\" requires at least one Argument, suiciding.\n\n";
			 print chr(7); #<-- This is beeping
			 print "  To get help type \"$0  h\"\n\n\n ";
			 exit;
		 }
	  }
	}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#  Checking some input options like 'e=txt' for extension filtering
	#_____________________________________________________________________
	for($i=0; $i< @_; $i++){
	  if($_[$i]=~/e=(\S+)/){
					push(@extension, $1);
	  }
	}

	#""""""""""""""""""""""""""""""""""""""""""""""""""
	#   Some DEFAULT $debug variables for debugging purposes
	#""""""""""""""""""""""""""""""""""""""""""""""""""
	&set_debug_option;
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#  If there is only one prompt arg. and is [-]*[hH][elp]*, it calls
	#   &default_help and exits
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	if( ( @ARGV == 1 ) && ($ARGV[0] =~ /^[\-]*[hH\?][elp ]*$/) ){
		&default_help;
		exit;
	}
	for($f=0; $f < @ARGV; $f++){
	 if( $ARGV[$f] =~ /\w+[\-\.\w]+$/ and -f $ARGV[$f] ){
		 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		 # When extension is defined, filter files by it
		 #____________________________________________
		 if(@extension > 0){
		     for($e=0; $e < @extension; $e++){
				 $extension=$extension[$e];
				 if($ARGV[$f]=~/\S\.$extension/){
					 push(@input_files, $ARGV[$f] );
				 }else{ next }
			 }
		 }else{
			 push(@input_files, $ARGV[$f] );
			 next;
		 }
	 }
	}

	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	#     Reading the running program script
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""
	&assign_options_to_variables;
	if($HELP == 1){ &default_help }
	return(\@input_files);
}

#________________________________________________________________________________
# Title     : import_ENV_vars
# Usage     :
# Function  : You can use any ENV set variables directly in your
#             program. So, you can say $USER instead of $ENV{'USER'}
# Example   :
# Keywords  : import_Env_vars, import_ENV_variables
# Options   :
# Version   : 1.1
#--------------------------------------------------------------------------------
sub import_ENV_vars{
		my($caller_package, $env_var_name);
		$caller_package=caller;
		foreach  $env_var_name (keys %ENV){
			 ${"${'caller_package'}::${'env_var_name'}"}=$ENV{$env_var_name}; ## ' ' are necessary
		}
		print "\n# import_ENV_vars: ALL the ENV settings are imported to $0 program\n";
}

#__________________________________________________________________________
# Title     : sort_hash_value_by_column
# Usage     : @out=@{&sort_by_column(\%input_line_hash, <column num>)};
# Function  : it sorts values of hash by the given column , small comes top. Unless number is
#             is given, it sorts by the first column.
#             It returnns ARRAY of the keys of the input HASH!!!
#
#             It can handle gzipped file. It called gunzip to open and sort.
#
# Example   : Above will sort the file xxxx.msp by its 3rd column(numerically)
#               small numbers will come to the top.
# Keywords  : sort_by_2nd_column, sort_by_second_column, sort_by_e_values,
#             sort_by_evalues, sort_hash_by_column, sort_value_by_column,
# Options   :
#      s  for sorting stringwise
#      d  for sorting by digit
#      n  for sorting by digit(numerically)
#   numerically  an alias of n
#
# Version   : 1.1
#----------------------------------------------------------------------------
sub sort_hash_value_by_column{
	 my (%in, $i, $col, $sort_numerically, $sort_non_numerically, @keys);
	 $sort_numerically=1;
	 if(@_ < 2  ){ print "\n# FATAL: sort_by_column needs 2 arguments\n"; exit }
	 for (@_){
			if(ref $_ eq 'HASH'){ %in =%{$_}; }
			elsif( ref $_ eq 'SCALAR'){ $col=${$_}; }
			elsif(/^ *\d+ *$/){ $col=$_ }
			elsif(/^ *[nd] *$/i){ $sort_numerically=1; $sort_non_numerically=0; }
			elsif(/^ *n[umerically]* *$/i){ $sort_numerically=1; $sort_non_numerically=0; }
			elsif(/^ *s *$/i){ $sort_non_numerically=1; $sort_numerically=0; }
	 }
	 $col--;
	 @keys= keys %in;
	 if($sort_numerically ){
			 @keys= map {$_->[0]} sort { $a->[1] <=> $b->[1] } map { [$_, ($in{$_}=~/(\S+)/g)[$col] ] } @keys;
	 }else{ # here let's do the sring sort
			 @keys= map {$_->[0]} sort { $a->[1] cmp $b->[1] } map { [$_, ($in{$_}=~/(\S+)/g)[$col] ] } @keys;
	 }
	 return(\@keys);
}

#________________________________________________________________________
# Title     : get_base_names
# Usage     : $base =${&get_base_names(\$file_name)};
#             :   or @bases = &get_base_names(\@files);  # <-- uses `pwd` for abs directory
# Function  : produces the file base name(eg, "evalign"  out of "evalign.pl" ).
#              when xxxx.xx.gz form file is given, it removes gz as well
#
# Example   : $base => 'test'  with 'test.txt' or '/home/dir/of/mine/text.txt'
# Warning   :
# Keywords  : get_base_name{, base_name, file_base_name ,  get_file_base_name
#             get_basename, basename, get_root_name
# Options   :
# Returns   :
# Argument  : handles both ref and non-ref.
# Version   : 1.5
#--------------------------------------------------------------------
sub get_base_names{
	my($x, $pos, $pos1, @out_file, $file_only, $file, @file, $base, @base);
	@file=@{$_[0]} if (ref($_[0]) eq 'ARRAY');
	@file=@_ if !(ref($_[0]) eq 'ARRAY');
	for($x=0; $x < @file; $x ++){
		if( ref($file[$x]) ){
			$file = ${$file[$x]};
			$pos1=rindex($file, "/");
	        $file_only=substr($file, ($pos1+1));
						if($file_only=~/(\S+\.\S+)\.gz$/){
							 $pos = rindex($1, ".");
			   $base= substr($1, 0, $pos);
			}elsif($file_only=~/^[^\.]+$/){ ## when file does not have '.' in its name
							 $base= $file_only;
						}else{
							 $pos = rindex($file_only, ".");
							 $base= substr($file_only, 0, $pos);
						}
		}else{
			$file = $file[$x];
			$pos1=rindex($file, "/");
	        $file_only=substr($file, ($pos1+1));
			if($file_only=~/(\S+\.\S+)\.gz$/){
							 $pos = rindex($1, ".");
			   $base= substr($1, 0, $pos);
						}elsif($file_only=~/^[^\.]+$/){ ## when file does not have '.' in its name
							 $base= $file_only;
			}else{
							 $pos = rindex($file_only, ".");
							 $base= substr($file_only, 0, $pos);
						}
		}
		push(@base, $base);
	}
	if(@base == 1 ){ \$base[0] }else{ \@base }
}
#________________________________________________________________________
# Title     : read_option_table
# Usage     :
# Function  : Reads the option table made by Jong in any perl script. The
#             option table is a box with separators.
# Example   :
# Warning   :
# Keywords  :
# Options   :
# Returns   :
# Argument  :
# Version   : 1.0
#--------------------------------------------------------------------
sub read_option_table{
	my($table_found, @option_tb, $head, );
	 open(SELF, "${$_[0]}");
	 while(<SELF>){
		if( (/^ *#+/) && ( $table_found== 1) ){
		  push (@option_tb, "$_");
		}elsif( ($table_found != 1)&&(/^ *\#+ *[Oo]ption *[Tt]able */) ){
			$table_found=1; $head="############## Option Table  \"$0\"\n"; ##
			push(@option_tb, $head);
		}
		if( ($table_found==1)&&(/^ *###################+ *$/)){  ### to find the end point of reading
			$table_found =0; last; }
	 }
	 return(\@option_tb);
}
#________________________________________________________________________
# Title     : read_head_box
# Usage     : %entries = %{&read_head_box([\$file_to_read, \@BOXED ] )};
# Function  : Reads the introductory header box(the one you see on top of sub routines of
#             Jong's programs.). Make a hash(associative array) to put entries
#             and descriptions of the items. The hash values have new lines '\n' are
#             attached, so that later write_head_box just sorts Title to the top
#             and prints without much calculation.
#             This is similar to read_head_box, but
#             This has one long straight string as value(no \n inside)
#             There are two types of ending line one is Jong's #---------- ...
#             the other is Astrid's  #*************** ...
# Example   : Output is something like
#             ('Title', 'read_head_box', 'Tips', 'Use to parse doc', ...)
# Warning   :
# Keywords  : open_head_box, open_headbox, read_headbox
# Options   : 'b' for remove blank lines. This will remove all the entries
#             with no descriptions
# Returns   : A hash ref.
# Argument  : One or None. If you give an argu. it should be a ref. of an ARRAY
#              or a filename, or ref. of a filename.
#             If no arg is given, it reads SELF, ie. the program itself.
# Version   : 2.7
#--------------------------------------------------------------------
sub read_head_box{
	my($i, $c, $d, $j, $s, $z, @whole_file, $title_found, %Final_out,
	  $variable_string, $TITLE, $title, @keys, $end_found, $line, $entry,
	  $entry_match, $End_line_num, $remove_blank,  $title_entry_null,
	  $end_found, $Enclosed_entry, $Enclosed_var, $blank_counter,
	  $title_entry_exist, $entry_value, $temp_W, $Warning_part, $tmp,
	  $option_tb_found
	);

	if(ref($_[0]) eq 'ARRAY'){ ## When array is given
	  @whole_file = @{$_[0]};
	}elsif(-e ${$_[0]}){       ## When filename is given in a ref
	  open(FILE, "${$_[0]}");
	  @whole_file=(<FILE>);
	}elsif(-e $_[0]){          ## When filename is given
	  open(FILE, "$_[0]");
	  @whole_file=(<FILE>);
	}elsif( $_[0] eq 'b'){          ## When filename is given
	  $remove_blank = 1;
	}elsif( ${$_[0]} eq 'b'){          ## When filename is given
	  $remove_blank = 1;
	}else{
	  open(SELF, "$0");
	  @whole_file=(<SELF>);
	}
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	for($i=0; $i<@whole_file; $i++){
	 $whole_file[$i] =~ tr/\t/ {7}/;  ## This is quite important to some parsing!!!
	 #########################################
	 ##  The first and second line of box 1 ##
	 #########################################
	 if( ($whole_file[$i]=~/^#[_\*\~\-\=]{20,}$/)&&    ##  '#______' is discarded
		 ($whole_file[$i+1]=~/ *\# {0,3}([TitlNam]+e) {0,8}: {1,10}([\w\.:]*) *(Copyright.*)/i) ){
		 $TITLE = $1;      $title = "$2\n";   $Final_out{'Warning'}.="$3\n";
		 $entry_match=$TITLE; ## The very first $entry_match is set to 'Title' to prevent null entry
		 if($TITLE =~ /^Title|Name$/i){   #
			  if( ($title=~/^\s+$/)||( $title eq "\n") ){
				  $title_entry_null =1;  $title = '';  }    }
		 $Final_out{$TITLE}=$title;
		 $title_found ++ ;   $i++;  ## << this is essential to prevent reading the same line again.
		 last if $title_found > 1;    }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ## The first and second line of box 2, #__________ or #**************
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&($whole_file[$i]=~/^#[_\*]{20,}$/)&&
		 ($whole_file[$i+1]=~/^# {1,3}(\w{1,6}\s{0,2}\w+) {0,7}: {1,5}(.*) */i) ){
		 $title_found ++ ;        $i++;
		 $entry_match=$1;       $entry_value=$2;
		 $entry_match =~ s#^\S#(($tmp = $&) =~ tr/[a-z]/[A-Z]/,$tmp)#e;  ## Capitalize words
		 $Final_out{$entry_match}.= "$entry_value\n";
		 last if $title_found > 1;  next;   }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ##  'Enclosed' : section. After this, everything is read without discrimination ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($Enclosed_entry == 1)&&($whole_file[$i] =~/^#{1,} {1,}(.*)$/) ){
		 $Final_out{$Enclosed_var}.= "$1\n";    }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ##  With proper entry 1 : for  'eg)'
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&($title_found==1)&&
		 ($whole_file[$i]=~ /^# {1,12}(eg ?\)) {0,8}(.*)/i)){
		 $entry_match='Example';
		 $Final_out{$entry_match}.= "$2\n";
	 }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ##  With PROPER entry 2 : descriptins like. 'Ussage : ssssssxxjkk  kj'
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&($title_found==1)&&
		 ($whole_file[$i]=~ /^# {1,2}(\w{1,4}\s{0,2}\w{1,7}) {0,8}[:\)] {0,6}(.*) */i)){
		 $entry_match=$1;       $entry_value=$2;
		 $entry_match =~ s#^\S#(($tmp = $&) =~ tr/[a-z]/[A-Z]/,$tmp)#e;
		 $Final_out{$entry_match}.= "$entry_value\n";
		 if($entry_match=~/^(Enclosed?)$/i){
			  $Enclosed_entry = 1;  $Enclosed_var=$1;        }    }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 #  With proper entry 3 : descriptins like. 'Ussage :', But blank description ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&($title_found==1)&&
		 ($whole_file[$i]=~ /^# {1,2}(\w{1,4}\s{0,2}\w{1,7}) {0,8}[:\)]( {0,})$/i)){
		 $entry_match=$1;       $entry_value=$2;
		 $entry_match =~ s#^\S#(($tmp = $&) =~ tr/[a-z]/[A-Z]/,$tmp)#e;
		 $Final_out{$entry_match}.= " $entry_value\n";
		 if($entry_match=~/^(Enclosed?)$/i){
			  $Enclosed_entry = 1;  $Enclosed_var=$1;      }    }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 #  $option variable matching                ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1) && ($title_found==1) &&
		 ($whole_file[$i]=~ /^# {1,15}([\$\@]\w+ +[\w=\>]+ +\S+ \w+ \S+ *.*)/ )){
		 $Final_out{$entry_match} .= "$1\n";  }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ###  all space line matching                 ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&  ##<----- If blank line is matched. Take the line
		 ($title_found==1)&&($whole_file[$i]=~/^# {0,}$/) ){
		 $blank_counter++;
		 if($blank_counter > 2){ $blank_counter--; }
		 else{ $Final_out{$entry_match}.= " \n";  }     }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ###  Anything after 3 space to 12 positions  ##
	 ###  To match 'examples' etc. INC. ':'       ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&
		 ($title_found==1)&&($whole_file[$i]=~/^#( {2,12})(.+)/) ){
		 $Final_out{$entry_match}.= "$1$2\n"; $blank_counter=0; }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ###  Anything after 1 space to 11 positions  ##
	 ###  To match 'examples' etc. EXC. ':'       ##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($end_found != 1)&&
		 ($title_found==1)&&($whole_file[$i]=~/^# {1,12}([^:.]+)/) ){
		 $Final_out{$entry_match}.= "$1\n"; $blank_counter=0;}

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 ###-------End of the read_box reading--------##
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( ($title_found==1)&&
		 ($whole_file[$i]=~ /^#[\~\=\*\-]{15,}/)){  ## to match '#-----..' or '#******..'(Astrid's)
		 $End_line_num = $i;       $end_found++;
		 last;      }

	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 #  <<<<  Check if there is option table >>>>  #
	 #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	 elsif( (/^#{10,} option table of this program   #{10,}/)&&($end_found >=1) &&($title_found==1)){
				 $option_tb_found++; ### This is a global var.
	 }
	} ## < End of for loop


	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	### If title is not there at all     ####
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	@keys=keys %Final_out;
	for(@keys){
	  if(/^Title$/i){    ## No Entry of Title at all??
		  $TITLE =$&;
		  $title_entry_exist = 1;
		  if($Final_out{$_}=~/^ *$/){   ## if Title => Null or just space
			  $title_entry_null = 1;    }  }  }

	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	### When title entry is not there    ####
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	if( $title_entry_exist != 1){
		for($s=$End_line_num+1; $s < $End_line_num+20; $s++){
			if( $whole_file[$s] =~ /^sub {1,5}([\w\.]+) {0,6}\{/){
				$Final_out{'Title'} = "$1\n";   last;       }
			elsif( $whole_file[$s] =~/^#________________________+/){
				#######################################
				## Uses running file name as titile  ##
				#######################################
				$Final_out{'Title'} = "$0";     last;
			}else{
				#######################################
				## Uses running file name as titile  ##
				#######################################
				$Final_out{'Title'} = "$0";
			}
		}
	}
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	### When title is blank              ####
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	elsif($title_entry_null ==1){  ## It looks for 'sub xxxxx{ ' line to get title
		### $End_line_num is the last line read.
		for($s = $End_line_num+1; $s < $End_line_num+20; $s++){
			if( $whole_file[$s] =~ /^sub {1,5}(\w+\.*\w*) {0,7}{/){
				$Final_out{$TITLE} = "$1\n";    last;     }
			elsif( $whole_file[$s] =~/^#________________________+/){
				#######################################
				## Uses running file name as titile  ##
				#######################################
				$Final_out{$TITLE} = "$0";     last;
			}else{
				#######################################
				## Uses running file name as titile  ##
				#######################################
				$Final_out{$TITLE} = "$0";
			}
		}
	}
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	## Error handling, if no head box is found   ####
	#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	if($title_found < 1){ print "\nFatal: No headbox found by read_head_box2 sub.\n";  }
	\%Final_out;
}               ##<<--- ENd of the sub read_head_box

#________________________________________________________________________
# Title     : show_hash
# Usage     : &show_hash(\@input_array);
# Function  : for debugging purpose. Shows any array elem line by line.
#             the line is 60 elements long (uses recursion)
# Example   : Output:      item1
#             Output:      item2
#             Output:      item3
# Warning   : There is a global variable:  $show_hash_option
#             It tries to detect any given sting which is joined by ','
# Keywords  :
# Options   : -s or -S or s or S for spaced output. Eg)
#             seq1       1 1 1 1 1 1 1 1 1 1 1 1
#
#             instead of
#             seq1       111111111111
#
#             -h or -H or h or H for horizontal line of '---------...'
#
# Returns   :
# Argument  :
# Category  :
# Version   : 1.7
#--------------------------------------------------------------------
sub show_hash{
	my($k, $i, $t, @in2, $in, $LEN, %TEM ); ## You should not put $show_hash_option
	my(@in)=@_;                     ## and $horizontal_line  in my !!!
	my($KL)=2; # default keys string length;
	my($VL)=80; # default values string length;
	my($GAP)=2;  # default space between keys and values
	my($horizontal_line, $show_hash_optionXX, $Hash_counter, @line);

	## This is to get the option of 'space' to make spaced output.
	for($t=0; $t < @in; $t++){
	 if($in[$t] =~/^[-]+[sS][pace]*$/){
		 $show_hash_optionXX = 1;
		 splice(@in, $t, 1);
	 }elsif(${in[$t]} =~/^[-]+[sS][pace]*$/){
		 $show_hash_optionXX = 1;
		 splice(@in, $t, 1);
	 }elsif($in[$t] =~/^[-]+[hH][rR]*$/){
		 $horizontal_line = 1;
		 splice(@in, $t, 1);
	 }elsif(${in[$t]} =~/^[-]+[hH][rR]*$/){
		 $horizontal_line = 1;
		 splice(@in, $t, 1);
	 }
	}

	######## Main loop #################
	if($horizontal_line ==1){  ## This puts the delimiter '--------------(  )'
	  $Hash_counter ++;
	  print "\n","-"x78,"(${Hash_counter}th hash)", "\n";
	}

	for($k=0; $k < @in; $k++){
	 if(ref($in[$k]) eq 'ARRAY'){  ### When the hashes were given in array ref.
		 &show_hash(@{$in[$k]}, $show_hash_optionXX, $horizontal_line);
		 print "\n";
	 }
	 elsif(ref($in[$k]) eq 'HASH'){  ### recursion
		 &show_hash(%{$in[$k]});
		 print "\n";
	 }
	 elsif(ref($in[$k+1]) eq 'HASH'){  ### recursion
		 &show_hash(%{$in[$k+1]}); print "\n";
	 }
	 elsif(ref($in[$k]) eq 'SCALAR'){  print ${$_[$k]}, "\n";  }
	 elsif( !ref($in[$k]) ){
		 if( !ref($in[$k+1]) && defined($in[$k+1])  ){
			 if($show_hash_optionXX == 1){  #### space option checking.
				#if($in[$k+1] =~ /\,.+\,/){  #### if the string is joined with ','
				#	 @line = split(/\,/, $_[$k+1]);
				# }else{
				#	 @line = split(//, $_[$k+1]);
				# }
				%TEM = @in;
				$LEN = ${&max_elem_string_array_show_hash(keys %TEM)};
				 if($LEN > $KL){ $KL = $LEN + $GAP +2};
				 printf ("%-${KL}s ", $in[$k]);  $k++;
				 printf ("%-${VL}s\n","@line");
			 }else{                        ### If not option is set, just write
				%TEM = @in;
				$LEN = ${&max_elem_string_array_show_hash( keys %TEM)};
				 if($LEN > $KL){ $KL = $LEN + $GAP +2};
				 printf ("%-${KL}s ", $in[$k]);  $k++; # print $in[$k], "\t";  $k++;
				 printf ("%-${VL}s\n",$in[$k]);        # print $in[$k], "\n";
			 }
		 }
		  #________________________________________________________
		  # Title    : max_elem_string_array_show_hash
		  # Keywords : largest string length of array
		  # Function : gets the largest string length of element of any array of numbers.
		  # Usage    : ($out1, $out2)=@{&max_elem_array(\@array1, \@array2)};
		  #            ($out1)       =${&max_elem_array(\@array1)          };
		  # Argument : numerical arrays
		  # returns  : one or more ref. for scalar numbers.
		  # Version  : 1.1
		  #-------------------------------------------------------
		  sub max_elem_string_array_show_hash{
			 my(@input, $i, $max_elem);
			 @input = @{$_[0]} || @_;
			 for($i=0; $i< @input ; $i++){
					$max_elem = length($input[0]);
					if (length($input[$i]) > $max_elem){
						$max_elem = length($input[$i]);
					}
			 }
			 \$max_elem;
		  }
		  #####################################insert_gaps_in_seq_hash
	 }
	}
}

#________________________________________________________________________
# Title     : show_options
# Usage     : &show_options;  usually with 'parse_arguments' sub.
# Function  :
# Example   :
# Keywords  : display_options, show_help_options, show_argument_options,
#             show_options_in_headbox, show_prompt_options
# Options   :
# Category  :
# Version   : 1.2
#--------------------------------------------------------------------
sub show_options{
	my($i, @keys, $perl_dir, $arg_num_limit, $head ,$arg_num_limit,
	 @entries_I_want_write );
	my($logname)=getlogin();
	my($pwd)=`pwd`;
	my($date)=`date`;
	chomp($date,$pwd);
	my($not_provided)="--- not provided ---\n";
	my($file_to_read) = $0;

	for($i=0; $i < @_; $i ++){
	  if((ref($_[$i]) eq 'SCALAR')&&(${$_[$i]} =~ /^\d$/)){
		  $arg_num_limit = ${$_[$i]};  }
	  elsif( (!(ref($_[$i]))) && ($_[$i] =~ /^\d$/)){
		  $arg_num_limit = $_[$i];     }
	}
	my %entries = %{&read_head_box(\$file_to_read )};
	if($option_tb_found ==1){
	 @option_tb=@{&read_option_table(\$file_to_read)};
	}
	foreach $help_item (keys %entries){
	  ${$help_item}= $not_provided if( (${$help_item}=~/^[\W]*$/)||( !defined(${$help_item})) );
	}
	#""""""""""""""""""""""""""""""""""""""""
	#########  Writing the format <<<<<<<<<<<
	#""""""""""""""""""""""""""""""""""""""""
	$~ =HEADER_HELP;
	write;   ## <<--  $~ is the selection operator
	$~ =DEFAULT_HELP_FORM;

	@entries_I_want_write=('Options');

	for( @entries_I_want_write ){  write  }

	print chr(7);  print "_"x72,"\n";

	if(@ARGV < $arg_num_limit){ print "\* $0 fataly needs $arg_num_limit arguments\n\n" }

	if(  $option_tb_found == 1){
	 #########  Printing the OPTION table contents <<<<<<<<<<<<
	 print "  Press \"Return\" key to see what options $logname\'s \n\n    \"$0\" take... \n";
		 $key_press=getc();
	 print @option_tb, "\n"x2 if(@option_tb > 0);
	}
format HEADER_HELP  =

**---------------------------------------------------------------------
	O P T I O N S  (I am &show_options)
**--------------------------------------------------------------------
.
format DEFAULT_HELP_FORM =
 @<<<<<<<<<: @*
 $_        $entries{$_}
.
}



#___________________________________________________________________
# Title     : scramble_array
# Usage     : @in=@{&scramble_array(\@in)};
# Function  : shuffles the elements of array
# Example   :
# Keywords  : randomise_array, randomize_array, shuffle_array
# Options   :
# Category  :
# Version   : 1.4
#---------------------------------------------------------------
sub scramble_array{
	srand(time()|$$);  # or use srand(time^$$);
	my ($i, @scrambled, @out, @each_array);

	for($i =0; $i< @_; $i++){
	   my @each_array = @{$_[$i]};
	   while (@each_array) {
		   push @scrambled, splice @each_array, int(rand(@each_array)), 1;
	   }
	   push(@out, \@scrambled);
	}
	if(@out > 1){
	   return(@out);
	}else{
	   return($out[0]);
	}
}


#________________________________________________________________________
# Title     : read_file_names_only
# Usage     : @all_files=@{&read_file_names_only(<dir>, [extension])};
# Function  : read any file names and REMOVES the '.', '..' and dir entries.
#             And then put in array.  This checks if anything is a real file.
#             You can use 'txt' as well as '.txt' as extension
#             You can put multiple file extension (txt, doc, ....)
#               and multiple dir path (/usr/Perl, /usr/local/Perl....)
#               It will fetch all files wanted in all the direc specified
#
#             It can handle file glob eg)
#             @all_files=@{&read_file_names_only(\$abs_path_dir_name, 'G1_*.txt')};
#               for all txt files starting with 'G1_'
#
# Example   : @all_files=@{&read_file_names_only(\$abs_path_dir_name, ..)};
#             @all_files=@{&read_file_names_only(\$dir1, '.pl', '.txt')};
#             @all_files=@{&read_file_names_only(\$dir1, '.', \$dir2, \$dir3, 'e=pl')};
#             @all_files=@{&read_file_names_only(\$abs_path_dir_name, 'G1_*.txt')};
#             @all_files=@{&read_file_names_only(\$abs_path_dir_name, \@target_file_names)};
#
# Warning   : This does not report '.', '..'
#             Only file names are reported. Compare with &read_any_dir
#             extension size should be less than 15 char.
#             It sorts the results!
# Keywords  : filename only, filename_only, read_files_only, read files
#             get_file_names_only, get_files_only, read_files_only
# Options   : "extension name". If you put , 'pl' as an option, it will show
#             files only with '.pl' extension.
#  '-p'      for path also included resulting in '/path/path/file.ext'
#              rather than 'file.ext' in output @array
#  '-s'      for sorting the results
#  e='xxx'  for extention xxx
#  '.pl'    for files extended by '.pl'
#  'pl'     for files extended by 'pl', same as above
#  D=       for dir name input
#  d=       for dir name input
#
# Category  :
# Version   : 3.0
#--------------------------------------------------------------------
sub read_file_names_only{
		my($in_dir, $i, $j, $x, $k, $dir, @final_files, @possible_dirs, $sort_opt, $ext, @extensions,
			 $path_include, @in, $glob_given, @files_globed, @in_dir, $pwd, $extension_given,
			 %target_file_names, @target_file_names, @read_files);
		$pwd=`pwd`; chomp($pwd);
		$in_dir=$pwd;
		@in=@_;

		print "\n# read_file_names_only: input are @in" if $verbose;

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#  Directory entry and opts detection
		#_________________________________________
		for($k=0; $k < @in; $k++){
			 if   ( $in[$k] eq '.'){ push(@in_dir,$pwd); splice(@in, $k, 1);  $k--; next }
					if( !(ref($in[$k]))){
	    print "\n# read_file_names_only: $in[$k] is not a reference";
				if($in[$k]=~/D=(\S+)/i){
						print "\n# read_file_names_only : $1 is used as input dir ";
						push(@in_dir, $1); splice(@in, $k, 1);    $k--; next;  }
				if( -d "$in[$k]" ){
						print "\n# read_file_names_only: $in[$k] is a dir";
						if($in[$k]=~/\/\S+$/){
								$path_include=1;  ## If the input dir has '/', I assume path should be added to out file names
								print "\n# read_file_names_only: \$path_include is set to 1";
			}
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# Removes the last slash '/' of input dir name
			#________________________________________________
			if($in[$k]=~/\/$/){   chop($in[$k]);  }
						push(@in_dir, $in[$k]);
						splice(@in, $k, 1);    $k--; next;
		}
		if(!(-f $in[$k]) and $in[$k] =~ /^\-p *$/ ){ ## somehow, ' *' is essential
			$path_include=1;
			splice(@in, $k, 1); $k--;
								}elsif(!(-f $in[$k]) and $in[$k] =~ /^\-s *$/   ){$sort_opt=1; splice(@in, $k, 1); $k--;
								}else{
										 print "\n# (W) read_file_names_only: $in[$k] not a file, nor dir, a file extnsion?\n";
								}
	 }elsif(ref($in[$k])){
				if(ref($in[$k]) eq 'SCALAR'){
					 if( -d ${$in[$k]}){
							 if(${$in[$k]}=~/\/$/){ chop(${$in[$k]}) }
							 push(@in_dir,${$in[$k]});
							 splice(@in, $k, 1);
							 $k--;
					 }elsif(!(-f $in[$k]) and ${$in[$k]} =~ /^\-p$/ ){$path_include=1; splice(@in, $k, 1); $k--;
					 }elsif(!(-f $in[$k]) and ${$in[$k]} =~ /^\-s$/ ){$sort_opt=1; splice(@in, $k, 1); $k--;}
				}elsif(ref($in[$k]) eq 'ARRAY'){
					 @target_file_names=@{$in[$k]}; splice(@in, $k, 1); $k--;
					 for($x=0; $x < @target_file_names; $x++){  # making a hash out of @array
							 $target_file_names{$target_file_names[$x]}=$target_file_names[$x];
					 }
				}
	 }
	}
	if(@in_dir < 1){ push(@in_dir, $pwd) }
	if($verbose){
		 print "\n# read_file_names_only: Final input directories are : @in_dir";
		 print "\n# read_file_names_only: going to \'File name and extension detection\' stage with \@in";
	}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#  File name and extension detection
	#_________________________________________
	for $dir (@in_dir){
        print "\n# read_file_names_only: changing to subdir \'$dir\'" if $verbose;
        chdir($dir);
        print "\n# read_file_names_only: trying to detect extension name from \@in: @in\n" if $verbose;

        for($k=0; $k < @in; $k++){
            if( !(ref($in[$k]))){
                 if($in[$k]=~/\*/){
                     $glob_given=1;
                     #~~~~~~~~~~~~~~~~~~~~~  Reads globbed files and attaches path if opt -p is set
                     if($path_include==1){  @final_files=map{ "$dir/$_" } <$in[$k]>;
                     }else{ @final_files=<$in[$k]>;  }
                     splice(@in, $k, 1); $k--;
                 }elsif(!(-f $in[$k]) and $in[$k] =~/e=\.?(\S+)/){ $extension_given =1; push(@extensions, $1); splice(@in, $k, 1);$k--;
                 }elsif(!(-f $in[$k]) and $in[$k] =~/\.*(\S+)/){
                     print "\n# read_file_names_only: pushing $1 as an extension" if $verbose;
                     $extension_given =1; push(@extensions, $1);
                     splice(@in, $k, 1); $k--;
                 }elsif(!(-f $in[$k]) and $in[$k] =~/^([^\-]{0,8})$/){  ## extension name can not be larger than 8 chars
                     print "\n# read_file_names_only: pushing $1 as an extension" if $verbose;
                     $extension_given =1; push(@extensions, $1);
                     splice(@in, $k, 1); $k--;
                 }
            }elsif(ref($in[$k])){
                  if(ref($in[$k]) eq 'SCALAR'){

                      if(${$in[$k]}=~/\*/){
                          $glob_given=1;
                          if($path_include==1){  @final_files=map{ "$dir/$_" } <${$in[$k]}>;
                          }else{ @final_files=<${$in[$k]}> }
                          splice(@in, $k, 1); $k--;
                      }elsif(!(-f ${$in[$k]}) and ${$in[$k]} =~/e=(\S+)/ ){ $extension_given = 1;
                                                      push(@extensions, $1); splice(@in, $k, 1);  $k--;
                      }elsif(!(-f ${$in[$k]}) and ${$in[$k]} =~/^\.?(\S+)/ ){$extension_given =1;
                                                      push(@extensions, $1);  splice(@in, $k, 1);  $k--;
                      }
                  }
            }
        }
        chdir($pwd);
	}
	if( $glob_given == 1 and  $extension_given !=1 ){  # when glob input is given only(without any extension input!
		 print "\n# read_file_names_only: You used glob for file name, but without extension name\n" if $verbose;
	 return(\@final_files);
	}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	#  Main READING PART
	#_________________________________________________
	print "\n# read_file_names_only: \@in_dir is  @in_dir\n";
	for($k=0; $k< @in_dir; $k++){
		 chdir($in_dir[$k]) or die "\n# read_file_names_only: could not get into $in_dir[$k]\n";
	 opendir(DIR1, ".");
				 @read_files = readdir(DIR1);
	 print "\n# read_file_names_only: content of \@read_files in $in_dir[$k] : @read_files\n" if $verbose;
	 if(@read_files < 1){ print "\n# read_file_names_only: ERROR??, \@read_files is empty\n\n\n"; }
	 for($i=0; $i < @read_files; $i ++){
						#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
						# If the user has specified the target file names
						#____________________________________________________
	    if( @target_file_names > 0){
           if( -f "$read_files[$i]" and -s $target_file_names{$read_files[$i]} ){ ##
                        if($extension_given ==1 ){
                                for $ext (@extensions){
                                           if( $read_files[$i] =~ /\.$ext$/){
                                                    if($path_include==1){
                                                        push(@final_files, "$in_dir[$k]\/$read_files[$i]" );
                                                    }else{
                                                        push(@final_files, "$read_files[$i]" );
                                                    }
                                           }
                                }
                        }else{ ## reading everything !!!
                            push(@final_files, $read_files[$i]);
                        }
               }
           }else{
               if( -f "$read_files[$i]" ){ ##
                    if($extension_given ==1 ){
                               for $ext (@extensions){
                                        if( $read_files[$i] =~ /\.?$ext$/){
                                                   if($path_include==1){
                                                        push(@final_files, "$in_dir[$k]\/$read_files[$i]" );
                                                   }else{
                                                        push(@final_files, "$read_files[$i]" );
                                                   }
                                        }
                               }
                    }else{ ## reading everything !!!
                         push(@final_files, $read_files[$i]);
                    }
               }
           }
	 }
	 chdir($pwd);
   }
   @final_files=sort @final_files if $sort_opt == 1;
   return(\@final_files);
}

#______________________________________________________________________________
# Title     : open_list_file
# Usage     :
# Function  :
# Example   :
# Keywords  : open_list_file_HASH
# Options   :
# Author    : jong@bio.cc,
# Category  :
# Returns   :
# Version   : 1.0
#------------------------------------------------------------------------------
sub open_list_file{
    my($list_file, %list);
    $list_file=${$_[0]} || $_[0];
    open(LIST_FILE, "$list_file") || die "\n Can not open $list_file \n";
    while(<LIST_FILE>){
       if(/(\S+)/){
          $list{$1}=$1;
       }
    }
    close(LIST_FILE);
    return(\%list);
}

#______________________________________________________________________________
# Title     : open_list_file_ARRAY
# Usage     :
# Function  :
# Example   :
# Keywords  : open_list_file_array
# Options   :
# Author    : jong@bio.cc,
# Category  :
# Returns   :
# Version   : 1.0
#------------------------------------------------------------------------------
sub open_list_file_ARRAY{
    my($list_file, @list);
    $list_file=${$_[0]} || $_[0];
    open(LIST_FILE, "$list_file") || die "\n Can not open $list_file \n";
    while(<LIST_FILE>){
       if(/(\S+)/){
          push(@list, $1);
       }
    }
    close(LIST_FILE);
    return(\@list);
}

#______________________________________________________________________________
# Title     : get_total_memory_size_in_linux
# Usage     : $mem=${&get_total_memory_size_in_linux};
# Function  :
# Example   : The /proc/meminfo file looks like this:>>>>
#           total:    used:    free:  shared: buffers:  cached:
#   Mem:  395735040 233975808 161759232 65953792 111476736 41345024
#   Swap:  7319552   147456  7172096
#   MemTotal:    386460 kB
#   MemFree:     157968 kB
#   MemShared:    64408 kB
#   Buffers:     108864 kB
#   Cached:       40376 kB
#   SwapTotal:     7148 kB
#   SwapFree:      7004 kB
#
# Keywords  : get_memory_size_in_linux, get_mem_size, check_memory_size
#             get_system_memory_size
# Options   :
# Author    : jong@bio.cc
# Category  :
# Version   : 1.1
#------------------------------------------------------------------------------
sub get_total_memory_size_in_linux{
    my($total_memory_size);
    open(MEMORY,"</proc/meminfo") || die "Something is very wrong, can't get memory size\n";
    while (<MEMORY>){
         if(/MemTotal:\s+(\S+)\s+kB/i){
                $total_memory_size=$1;
         }
    }
    close(MEM);
    return(\$total_memory_size);
}


