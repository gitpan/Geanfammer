#!/usr/bin/perl
# Last Update by /gn0/jong/Perl/update_subroutines.pl: Mon Sep  1 18:37:22 BST 1997
#____________________________________________________________________________
# Title     : diviclus.pl
# Usage     : diviclus.pl xxxxx.msp [s=100][t=30][e=30][f=2][v]
#             diviclus.pl *.msp  <-- For multiple MSP file processing
#                           -while xxx.msp is a file with MSP file format
# Function  : divides complex single linkage cluster into smaller duplication
#               module level sub clusters.
#             1) merges similar sequences in MSP file fomrat making
#                hash output.
#             2) processes the merged msp contents to sort small things
#             3) connects the sequences when there is any common
#                region between preliminary clusters.
#             4) makes various files(xx.clu, xx.sat, xx.mrg) and
#                also shows in STDOUT.
#
#             To controll the division of clusters, play with
#              the below parameters(if you do not specify, defaults:
#               t=30, s=100, e=0.01, f=7
#              if you give them mulitiple files, it will process them
#              all together.
# Example   : diviclus.pl xxxxxx.msp s=90 t=40 e=0.001
#             Above is for score 90, seqlet leng 40, evalue 0.001.
#             However, usually you dont need options. Just put xxx.msp
#
# Keywords  : divide_clusters, diviclust, find_linker, subcluster, divicl,
# Options   :
#   f=<digit>   for determing the factor in filtering out non-homologous
#                  regions, 7 = 70% now!!
#   l=<digit>   for seqlet(duplication module) length threshold
#   t=<digit>   for seqlet(duplication module) length threshold
#                  (same as l opt, confusing, huh? )
#   s=<digit>   for score threshold
#   e=<digit>   for evalue threshold
#   z           for activating remove_similar_sequences, rather than remove_dup....
#   o           for overwriting
#   v           for verbose printout (for debugging information for Sarah
#                              and Jong, NOT for John, Jason, Tom, Jessy, Pat,,)
#   d           for dynamic factor determination(if the original cluster size
#                    is very big(say, over 100 seqs), it increases the factor by
#                    10% etc.(above default or given factor)
#
#   S  $short_region=  S by S -S  # taking shorter region overlap in removing similar reg
#   L  $large_region=  L by L -L  # taking larger  region overlap in removing similar reg
#   A  $average_region=A by A -A  # taking average region overlap in removing similar reg
#
#  $short_region=  S by S -S # taking shorter region overlapped
#  $large_region=  L by L -L # taking larger  region overlapped
#  $average_region=A by A -A # taking average region overlapped
#  $verbose = v by v -v
#  $range = r by r -r
#  $merge = m by m -m
#  $sat_file = s by s -s
#  $dindom = d by d -d
#  $indup  = i by i -i
#  $over_write = w by w -w
#  $optimize   = o by o -o
#  $dynamic_factor=D by D -D
#  $score  = by s=   # Ssearch score cutoff, default 100
#  $factor = by f=   # factor is for the merge proces
#                      (misoverlap tolerance factor 3=33%, 2=50%)
#                      factor works within msp chunk for one sequence
#                      to filter a good mergable seqlets
#  $thresh = by t=   # seqlet length cutoff, default 30
#  $evalue = by e=   # maximum evalue cutoff default 30
#
# Author    : J. Park, Sarah Teichmann, sat@mrc-lmb.cam.ac.uk, jong@salt2.med.harvard.edu
# Version   : 2.4
#------------------------------------------------------------------------------

#============= Default parameters used ==========================
   $optimize='o'; ## default is set to optimized (for merging-> removing similar seqs in merging )
   $factor=7;     ## default is 7 (70% overlap is thought to be good), now 7 means 70%
   $range='r';    ## default is to write ranges in the sequences
   $merge='m';    ## default is to speed up clustering by merging
   $score=75;     ## default ssearch(or any smith waterman) score cutoff is 110 (rather low)
   $evalue=0.001; ## default ssearch evalue cutoff is 0.001 (rather stringent)
   $thresh=30;    ## default minimum seqlet length as a possible domain(cf. 40aa for DALI)
   $dynamic_factor=''; # default dynamic factor is not set.
   #$indup ='i';   # now defunt, do not use please
   #$dindom='d';   # now defunt, do not use please
   #$sat_file='s'; # now defunt, do not use please
#============= Default parameters used ==========================

my $good_bad;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (1) Getting the input MSP (matching seq pair) format file
#____________________________________________________________
my @file=@{&parse_arguments(1)};


#============= BELOW: Actual Run of divide_clusters which is the main sub of divclus ====================
&divide_clusters(
      \@file,
	  "s=$score",
      "f=$factor", ## this is a very impo. parameter in the behaviour of divclus, Sarah!
	  "t=$thresh",
	  "e=$evalue", ## this is a very impo. parameter in the behaviour of divclus, Sarah!
      $dynamic_factor,
      $verbose,
      $range,
      $merge,
      $sat_file,
	  $dindom,
	  $indup,
	  $over_write,
	  $optimize,
	  $short_region,
	  $large_region,
	  $average_region
);
#============= ABOVE: Actual Run of divide_clusters ====================

# Sub routine listing starts

#_______________________________________________________________________
# Title     : divide_clusters
# Usage     : &divide_clusters(\@file);
# Function  : This is the main funciton for divclus.pl
#               divides complex single linkage cluster into smaller duplication
#               module level sub clusters.
# Example   : &divide_clusters(\@file, $verbose, $range, $merge, $sat_file,
# 	                $dindom, $indup, "t=$thresh", "e=$evalue", $over_write,
#                   $optimize, "s=$score", "f=$factor");
#
# Keywords  : divicl, divclus, div_clus, divide clusters
# Options   : _  for debugging.
#   f=<digit>   for determing the factor in filtering out non-homologous
#                  regions, 7 = 70% now!!
#   l=<digit>   for seqlet(duplication module) length threshold
#   t=<digit>   for seqlet(duplication module) length threshold
#                  (same as l opt, confusing, huh? )
#   s=<digit>   for score threshold
#   e=<digit>   for evalue threshold
#   z           for activating remove_similar_sequences, rather than remove_dup....
#   o           for overwriting
#   v           for verbose printout (infor)
#   D           for dynamic factor
#   S  $short_region=  S by S -S  # taking shorter region overlap in removing similar reg
#   L  $large_region=  L by L -L  # taking larger  region overlap in removing similar reg
#   A  $average_region=A by A -A  # taking average region overlap in removing similar reg
#
# Version   : 2.4
#------------------------------------------------------------------------
sub divide_clusters{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my($merge, $verbose, $sat_file, $thresh, $factor, $indup, $indup_percent,
	   $score, @temp_show_sub, $optimize, $file, $evalue, $over_write, $din_dom,
	   $sum_seq_num, $base_1, $output_clu_file, $short_region, $large_region,
       $average_region, $dynamic_factor);

	$factor=7; # default factor is 7 for 70%

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dealing with options
    #_________________________________________
	if($char_opt=~/m/){	       $merge='m';
	}if($char_opt=~/v/){       $verbose='v'; # for showing debugging information
	}if($char_opt=~/i/){	   $indup='i';
	}if($char_opt=~/z/){	   $optimize='z';
	}if($char_opt=~/o/){       $over_write='o';
	}if($char_opt=~/d/){	   $din_dom='d';
    }if($char_opt=~/s/){       $sat_file='s';
    }if($char_opt=~/D/){       $dynamic_factor='D';
	}if($char_opt=~/S/){       $short_region  ='S';
	}if($char_opt=~/L/){	   $large_region  ='L';
    }if($char_opt=~/A/){       $average_region='A';
    }if($vars{'t'}=~/\d+/){    $thresh= $vars{'t'};
    }if($vars{'l'}=~/\d+/){    $thresh= $vars{'l'}; ## synonym of 't'
	}if($vars{'f'}=~/\d+/){    $factor= $vars{'f'};
	}if($vars{'s'}=~/\d+/){	   $score = $vars{'s'};
    }if($vars{'e'}=~/\d+/){    $evalue= $vars{'e'};
    }if($vars{'E'}=~/\d+/){    $evalue= $vars{'E'}; # synonym of e
    }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # When one file input was given (yes, divclus can handle multiple files, Sarah!)
   #________________________________________________________________________________
   if(@file == 1){  #<=== @file has xxxx.msp, yyyy.msp  zzzz.msp ....,
        print "\n# (1) divide_clusters: One single file was given=> \"@file\"\n" if $verbose;
		$file=$file[0];
		$base_1=${&get_base_names($file)};

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# (2) Define the output cluster file name:  eg, 3-232_cluster_F7.clu , F7 means factor used is 7
		#______________________________________________________________________________________________
		$output_clu_file="$base_1\_F${factor}\.clu";

		if( !$over_write and -s $output_clu_file){
			print "\n# $output_clu_file Already EXISTS, skipping. Use \'w\' opt to overwrite\n"; exit;
		}
		if($verbose){
		    print "# (2) divide_clusters: processing ONE single file \"@file\" with merge_sequence_in_msp_file\n";
        }

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # (3) merge_sequence_in_msp_file does not do much. Just filtering and producing
        #     sequences in ISPA_PBS_21-215 VPR_PBS_160-354 format from msp format
        #________________________________________________________________________________
        @out=@{&merge_sequence_in_msp_file(\@file, "s=$score", $optimize, $din_dom, $sat_file,
		    $optimize, "t=$thresh", "e=$evalue", "f=$factor", "$range", "$merge", $verbose,
		    $short_region, $large_region, $average_region, $over_write, $dynamic_factor)};

        if($verbose){
            print "\n\n# (3) divide_clusters: finished running \"merge_sequence_in_msp_file\"\n  ==>";
            for($i=0; $i< @out; $i++){
               print "\n-->> $out[$i]";
            }
        }

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # (4) This is critical seqlet merging step
        #________________________________________________________________________________
        @out=@{&cluster_merged_seqlet_sets(\@out, $dynamic_factor,
               "f=$factor",
		       $short_region, $large_region, $average_region, $optimize)};
		$percent_fac=$factor*10;


        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # (5) This is showing the result in clu file format
        #________________________________________________________________________________
        @temp_show_sub=&show_subclusterings(\@out, $file, $sat_file, $dindom, $indup,
						   "e=$evalue", "p=$percent_fac", "f=$factor" );
		$good_bad       = $temp_show_sub[0];
		$indup_c        = $temp_show_sub[1];
		$sum_seq_num   += $temp_show_sub[2];

		if($good_bad==1){      push(@good, $file);
		}else{                 push(@bad, $file);       }

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # (6) Final write up stage (unecessary)
        #_______________________________________________________________
		&write_good_bad_list_in_divide_clusters(\@good, \@bad);

   }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # when more than one single file input is given
   #____________________________________________________________
   elsif(@file >1 ){
		my (@good, @bad);
		if($indup =~/i/i){   open (INDUP, ">indup_stat\.txt");  } # this is not essential.

		for($i=0; $i< @file; $i++){
		     my (@out, @temp_show_sub);
			 my $indup_c=0;
			 $file=$file[$i];
			 $base_1=${&get_base_names($file)};
			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 # Define the output cluster file name:  eg, 3-232_cluster_F7.clu , F7 means factor used is 7
			 #______________________________________________________________________________________________
			 $output_clu_file="$base_1\_F${factor}\.clu";

			 if( !$over_write and -s $output_clu_file){
			     print "\n# $output_clu_file Already EXISTS, skipping. Use \'w\' opt to overwrite\n";
				 next;  }

			 #print "\n# (1)  divide_clusters: processing file \"$file\" for $output_clu_file";

			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 #  If clu file(eg 2-1618_ss.clu ) is in pwd, tries to skip
			 #____________________________________________________________
			 if((-s $output_clu_file) > 100 and $over_write !~/o/){
				print "# $output_clu_file exists, skipping, use \"w\" option to overwrite\n";  next;
			 }

			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 #  Merging similar sequences
			 #____________________________________________________________
             #print "\n# divide_clusters : I am merging seq in msp files \n";
             @out=@{&merge_sequence_in_msp_file(\$file, "s=$score", $din_dom, $sat_file, $optimize,
				 "t=$thresh", "e=$evalue", "f=$factor", "$range", "$merge", $verbose, $over_write,
				  $short_region, $large_region, $average_region, $dynamic_factor )};
             #print "\n# \@file has more than one input file\n # The result of \"merge_sequence_in_msp_file\"\n";
             #print "@out";

			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 #  Clustering the sets of merged seqlets => CORE algorithm
			 #____________________________________________________________

			 @out=@{&cluster_merged_seqlet_sets(\@out, "f=$factor", $optimize,
			        $short_region, $large_region, $average_region)};
			 $percent_fac=int(100-(1/$factor)*100);
			 @temp_show_sub=&show_subclusterings(\@out, $file, $sat_file, $dindom, $indup,
			                             "e=$evalue", "p=$percent_fac", "f=$factor");
			 $good_bad       = $temp_show_sub[0];
			 $indup_c        = $temp_show_sub[1];
			 $sum_seq_num   += $temp_show_sub[2];

			 if($good_bad==1){          push(@good, $file);
			 }else{         push(@bad, $file);       }

		}
		######################################### Writing stuff ##
		&write_good_bad_list_in_divide_clusters(\@good, \@bad);
		sub write_good_bad_list_in_divide_clusters{
			my  (@good, @bad);
			@good=@{$_[0]}; @bad=@{$_[1]};
			open(GOODBAD, ">good_bad.list");
			print GOODBAD "GOOD: all link    : 000\n";
			for($i=0; $i< @good; $i++){
			   print GOODBAD "$good[$i]\n";
			}
			print GOODBAD "BAD : Not all link: 000\n";
			for($i=0; $i< @bad; $i++){
			   print GOODBAD "$bad[$i]\n";
			}
			close(GOODBAD);
		}
		##########################################################

   }
}

#______________________________________________________________
# Title     : get_overlapping_range
# Usage     : @n1=@{&get_overlapping_range(\@ranges1, \@ranges2)};
# Function  :
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  : get_overlapping_range_in_msp, get_overlapping_range_in_msp_file,
#             get_overlapping_seq_match_range, get_overlap_seq_match_range
# Options   : _  for debugging.
#             #  for debugging.
# Returns   :
# Argument  :
# Version   : 1.1
#--------------------------------------------------------------
sub get_overlapping_range{
   my (@new_range, $R_start1, $R_start2);
   ($R_start1, $R_end1)=@{$_[0]}[0..1];
   ($R_start2, $R_end2)=@{$_[1]}[0..1];

   if(($R_start1 <= $R_start2)&&        # ------------
	 ( $R_end1 >= $R_end2) ){           #   -------
	   @new_range= ($R_start2, $R_end2);
   }elsif(($R_start1 <= $R_start2)&&    # -----------
	 ( $R_end1 <= $R_end2) &&           #    -----------
	 ( $R_end1 >  $R_start2) ){
	   @new_range= ($R_start2, $R_end1);
   }elsif(($R_start1 >= $R_start2)&&    #    -----------
	 ( $R_end1 >= $R_end2  ) &&         # -----------
	 ( $R_end2 >  $R_start1) ){
	   @new_range= ($R_start1, $R_end2);
   }elsif(($R_start1 >= $R_start2)&&    #   ------
	 ( $R_end1 <= $R_end2) ){           # -----------
	   @new_range= ($R_start1, $R_end1);
   }else{                                #  ----
	  @new_range=(0,0);                  #        --------
   }
   return(\@new_range);
}


#______________________________________________________________
# Title     : get_msp_enquiry_sequence
# Usage     :
# Function  : gets the name of sequence used as enquiry(target)
# Example   :
# Warning   :
# Class     :
# Keywords  : get_msp_target_sequence, get_msp_enquiry_sequence_name
# Options   : _  for debugging.
#             #  for debugging.
# Package   : Bio
# Reference : http://sonja.acad.cai.cam.ac.uk/perl_for_bio.html
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Scientist
# Version   : 1.0
#--------------------------------------------------------------
sub get_msp_enquiry_sequence{
   my $lines1=${$_[0]} || $_[0];
   my ($SEQ, $matched_SEQ);
   if($lines1 =~/^ *\d+ +\d+\.?[e\-\d]* +(\d+) +(\d+) +(\S+) +(\d+) +(\d+) +(\S+)/){
	  $SEQ        =$3;
	  $matched_SEQ=$6;
   }
   return \$SEQ;
}



#______________________________________________________________
# Title     : get_msp_matched_sequence
# Usage     :
# Function  : gets the name of sequence used as enquiry(target)
# Example   :
# Warning   :
# Class     :
# Keywords  : get_msp_matched_sequence_name
# Options   : _  for debugging.
#             #  for debugging.
# Package   : Bio
# Reference : http://sonja.acad.cai.cam.ac.uk/perl_for_bio.html
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Scientist
# Version   : 1.0
# Used in   :
# Enclosed  :
#--------------------------------------------------------------
sub get_msp_matched_sequence{
   my $lines1=${$_[0]} || $_[0];
   my ($SEQ, $matched_SEQ);
   if($lines1 =~/^ *\d+ +\d+\.?[e\-\d]* +(\d+) +(\d+) +(\S+) +(\d+) +(\d+) +(\S+)/){
	  $SEQ        =$3;
	  $matched_SEQ=$6;
   }
   return \$matched_SEQ;
}
#______________________________________________________________
# Title     : get_msp_range
# Usage     : @range=@{&get_msp_range($seqlet)};
#             @temp=&get_msp_range($seqlet);
#
# Function  :
# Example   :
# Warning   :
# Class     :
# Keywords  : get_msp_file_ranges
# Options   : _  for debugging.
#             #  for debugging.
# Package   : Bio
# Reference : http://sonja.acad.cai.cam.ac.uk/perl_for_bio.html
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Scientist
# Version   : 1.5
# Used in   :
# Enclosed  :
#--------------------------------------------------------------
sub get_msp_range{
   my $lines1=${$_[0]} || $_[0];
   my ($SEQ, $num_seq, $matched_SEQ, @Ranges);
   if($lines1 =~/^ *\d+ +\d+\.?[e\-\d]* +(\d+) +(\d+) +(\S+) +(\d+) +(\d+) +(\S+)/){
	  $SEQ        =$3;
	  $matched_SEQ=$6;
	  if($SEQ eq $matched_SEQ){ ## skipping self match
		  $num_seq++;
	  }else{
		  @Ranges=($1, $2, $4, $5);  ## <-- example. (10-20, 30-45)
	  }
   }
   return wantarray ? (\@Ranges, \$SEQ, \$matched_SEQ): \@Ranges;
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
  my($i, $perl_dir, $arg_num_limit, $head ,$arg_num_limit,
	  @entries, @entries_I_want_write );
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
 $_        $entries{$_}
.
}


#________________________________________________________________________________________
# Title     : merge_sequence_in_msp_file
# Usage     :
# Function  :
# Example   : INPUT: (MSP file) ===>
#  59     2.6        47    64     d2pia_3                    10    30     d1erd___10-30
#  161    1.1e-07    24    91     d2pia_3                    16    85     d1frd___16-85
#
#  722    0          1     106    d1put__                    1     106    d1put___1-106
#  66     4.9        2     68     d1put__                    43    106    d2lbp___43-106
#  69     1.3        12    49     d1put__                    81    120    d1cgo___81-120
#
#  60     3.3        13    38     d1frd__                    32    57     d1orda1_32-57
#  65     1.7        21    58     d1frd__                    40    69     d2mtac__40-69
#
#   ==== OUTPUT ===>
#    d1frd___1-98 d1frd___1-98_1-98 d1frd___16-85 d2pia_3_24-91_24-91
#    d1frd___16-85_16-85 d2pia_3_24-91
#    d1put___1-106 d1put___1-106_1-106
#    d2pia_3_1-98 d2pia_3_1-98_1-98
#
# Keywords  : mergr_seq_in_msp_file, merge_sequence_in_msp, merge_sequences_in_msp_file
# Options   :
#  $dynamic_factor = d by d -d   # adjusting factor value dynamically(more seq higher factor)
#  $short_region   =  S by S -S  # taking shorter region overlapped in removing similar regions
#  $large_region   =  L by L -L  # taking larger  region overlapped in removing similar regions
#  $average_region =  A by A -A # taking average region overlapped in removing similar regions
#
# Version   : 2.4
#----------------------------------------------------------------------------------------
sub merge_sequence_in_msp_file{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
	my ($msp_value, @all_seqlets, %temp_hash, @msp_chunks, $clu_out, $size_of_all_seqlets,
	    $base, $optimize, $mrg_out, @arr, $sat_out, %final_hash_out, @final_pre_hash,
        $thresh, $merge, $factor, $evalue, $score, $dynamic_factor, $given_factor,
		$short_region, $large_region, $average_region, $original_clu_size);
    $factor=$default_factor=7; #~~~~ default connection factor U, 7 means 70% now!
	$thresh=30;
	$evalue=10;
	$score =75;

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Following changes the defaults with given parameters
	#_____________________________________________________________
    if($char_opt=~/z/i){       $optimize='z';    ## This will cause using remove_similar_seqlets than remove_dup_in_array !
	}if($char_opt=~/m/){	   $merge='m';
    }if($char_opt=~/D/){       $dynamic_factor='D';
    }if($char_opt=~/v/){       $verbose='v';
	}if($char_opt=~/S/){       $short_region='S';
	}if($char_opt=~/L/){	   $large_region='L';
	}if($char_opt=~/A/){	   $average_region='A';
	}if($vars{'t'}=~/\d+/){	   $thresh=$vars{'t'};
	}if($vars{'f'}=~/\d+/){	   $factor=$vars{'f'};  ## Here I give a generous $factor !
	}if($vars{'s'}=~/\d+/){	   $score = $vars{'s'};
	}if($vars{'e'}=~/\d+/){	   $evalue= $vars{'e'};	}

    $given_factor=$factor;

   	for($c=0; $c< @file; $c++){
		open(MSP, "$file[$c]");
		$base=${&get_base_names($file[$c])};
		$clu_out="$base.clu"; # <-- This is the most important output. Sarah's program will process this
		$sat_out="$base.sat";
        print "# (2) merge_sequence_in_msp_file : processing $file[$c] for $clu_out\n" if $verbose;
		my @msp1=<MSP>;

		for($i=0; $i< @msp1; $i++){
            #~~~~~~~~~~ Include range or NOT in the seq name ~~~~~~~~~~~~~~~~~~~~~~~~~~`
			# %temp_hash is just to get the chunk of MSP block. As msp file uses empty line as a delimiter
			#____________________________________________________________________________
			if($char_opt=~/r/){
                if($msp1[$i]=~/^ *(\d+) +(\S+) *\S* +(\d+) +(\d+) +(\S+)[_\d+\-\d+]?[ \t]+(\d+)[ \t]+(\d+)[ \t]+(\S+)[_\d+\-\d+]?/){
                     $new_seq1="$5\_$3\-$4";
					 $new_seq2="$8\_$6\-$7";
                     if($new_seq1 eq $new_seq2){
                         print "# (3.0) merge_sequence_in_msp_file: Skipped = $msp1[$i]"  if $verbose;
                         next };

                     if($1 < $score or $2 > $evalue){
                         print "# (3.1) merge_sequence_in_msp_file: Skipped = $msp1[$i]"  if $verbose;
                         next };
                     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					 # Modifying $msp1[$i] line !!!
					 #______________________________
                     $msp1[$i]=sprintf("%s %-3s %s %s %s %s %s %s",
										$1, $2, $3, $4, $new_seq1, $6, $7, $new_seq2);
					 $temp_hash{$5}.="$msp1[$i]\n";
				}
			}else{
                if($msp1[$i]=~/^ *(\d+)[ \t]+(\S+)[ \t]*\S*[ \t]+\d+[ \t]+\d+[ \t]+(\S+)[_\d+\-\d+]? +\d+[\t ]+\d+[ \t]+\S+/){
                     if($1 < $score or $2 > $evalue){
                     print "# (3.0) merge_sequence_in_msp_file: Skipped = $msp1[$i]"  if $verbose;
                     next };
					 $temp_hash{$3}.="$msp1[$i]\n";
				}
			}#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		}
		close(MSP);
	}
    $original_clu_size=@msp_chunks= values(%temp_hash); ## Using temp hash is more than 2 times faster than push

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Determining the dynamic factor here (when 'd' opt is set)
    #____________________________________________________________
    if($dynamic_factor){
        #--> 100 => 10, 1000 => 15, 10000 => 20
        print "\n# ### Original \$factor: $factor";
        $factor += (log($original_clu_size)*5)/10 - 1; ## This is a simplistic.
        if($factor > 9.5){ $factor=9.5 } # this is the very upper limit for any factor.
        print "\n# ### After adjusted, \$factor: $factor\n";
    }

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#  Just to inform what parameters have been chosen
	#_____________________________________________________________
    print "\n# (4) merge_sequence_in_msp_file : default \$score  : $score";
	print "\n#                                : default \$evalue : $evalue";
	print "\n#                                : used    \$thresh : $thresh";
    print "\n#                                : given   \$factor : $given_factor";
    print "\n#                                : used    \$factor : $factor";
    print "\n#                                : \$dynamic_factor : $dynamic_factor\n";

    if($verbose){
        print "# (5) merge_sequence_in_msp_file: The msp chunks used are:\n";
        for($i=0; $i< @msp_chunks; $i++){
            print "$msp_chunks[$i]\n";
        }
    }

	for($i=0; $i< @msp_chunks; $i++){
       print "\n# (6) merge_sequence_in_msp_file: Processing eash chunk with merge_sequence_in_msp_chunk\n"  if $verbose;
       @arr=@{&merge_sequence_in_msp_chunk($msp_chunks[$i], $verbose, $optimize,
	         "$merge", "e=$evalue", "s=$score", "f=$factor", "t=$thresh",
	         $short_region, $large_region, $average_region)};
       push(@all_seqlets,  @arr);
	}

	#~~~~~~~~~ sorting inner sequences in strings ~~~~~~~~~
	#______________________________________________________
    @all_seqlets=@{&sort_words_in_string(@all_seqlets)}; ## This speeds up about 2 times !!!

	#~~~~~~~ Sort by the _digit-  in seqlet names ~~~~~~~~~
	@all_seqlets= map{$_->[0]} sort{$a->[1] cmp $b->[1] or $a->[2] <=> $b->[2]  }
			             map {/^ *((\S+)_(\d+)\-(\d+).*)/ && [$1, $2, $3, $4]} @all_seqlets;
    if( $verbose){
        print "\n# (7) merge_sequence_in_msp_chunk: Showing the very final result before merging\n";
        for($i=0; $i< @all_seqlets; $i++){
            print "\n$all_seqlets[$i]";
        }
    }

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# merge sequences in a simple way until there is no change in the array
	#  This is an incomplete merge(merges first seqlets of string ...
	#______________________________________________________________________
	for($i=0; $i< @msp_chunks; $i ++){
        ITERATION_RETURN_POINT:
        $size_of_all_seqlets=@all_seqlets;
        @all_seqlets = @{&merge_similar_seqlets(\@all_seqlets, $optimize,
		                 $short_region, $large_region, $average_region, "f=$factor")};
		if($size_of_all_seqlets > @all_seqlets){
		   @all_seqlets = @{&merge_similar_seqlets(\@all_seqlets, $optimize,
		                 $short_region, $large_region, $average_region, "f=$factor")};
           goto ITERATION_RETURN_POINT;
		}else{
		   last;
		}
	}

    if($optimize){
		@all_seqlets=@{&remove_similar_seqlets(\@all_seqlets,
		                $short_region, $large_region, $average_region)};
	}else{
		@all_seqlets=@{&remove_dup_in_array(\@all_seqlets)};
	}
	return(\@all_seqlets);
}



#_______________________________________________________________________________
# Title     : add_ranges_in_msp_line
# Usage     :
# Function  : this adds ranges to the seqnames of msp files
#             mmp line is msp line with additional sequences at the end
# Example   :
# Keywords  : convert_msp_to_mmp, convert_msp, convert_msp_2_mmp
#             change_msp_to_mmp, add_range_in_msp, convert_msp_line_to_mmp_line
# Options   : _  for debugging.
#             #  for debugging.
# Returns   :
# Argument  :
# Version   : 1.5
#-------------------------------------------------------------------------------
sub add_ranges_in_msp_line{
   my $input_msp=${$_[0]} || $_[0];
   my($score, $evalue, $long_1, $new_seq1, $new_seq2, $middle,
	  $start1, $end1, $start2, $end2, $seq1, $seq2, $new);

   if($input_msp=~/^ *(\d+) +(\S+) *\S*[ \t]+(\d+)[ \t]+(\d+)[ \t]+(\S+)[ \t]+(\d+)[ \t]+(\d+)[ \t]+(\S+)/){
	  ($score, $evalue, $start1, $end1, $start2, $end2)=($1, $2, $3, $4, $6, $7);
	  ($seq1, $seq2)=($5, $8);
	  if($seq1=~/(\S+)\_\d+\-\d+/){
		 $new_seq1="$1\_$start1\-$end1";
	  }else{
		 $new_seq1="$seq1\_$start1\-$end1";
	  }
	  if($seq2=~/(\S+)\_\d+\-\d+/){
		 $new_seq2="$1\_$start2\-$end2";
	  }else{
		 $new_seq2="$seq2\_$start2\-$end2";
	  }
	  $new=sprintf("%-6s %-8s %-5s %-5s %-32s %-5s %-5s %-32s",
					$score, $evalue, $start1, $end1, $new_seq1, $start2, $end2, $new_seq2);
   }
   return(\$new);
}




#______________________________________________________________
# Title     : sort_by_digits_in_string
# Usage     :
# Function  : sorts arrays of strings like
#
#   MJ0228_314-573 MJ1197_348-601
#   MJ0228_451-576 sll0078_502-594 sll1425_489-611
#   MJ0228_479-572 sll0078_502-594
#
#   According to the digits after seq names _314-, _451-, _479-
#    in the above
#   This only looks at the very first sequence in the string
#
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  :
# Options   : _  for debugging.
#             #  for debugging.
# Returns   :
# Argument  :
# Version   : 1.4
#--------------------------------------------------------------
sub sort_by_digits_in_string{
   my (@out, $i,  @temp1, @temp2, $old, @T);
   my @array_of_string=sort @{$_[0]};

   for($i=0; $i<= @array_of_string; $i++){
	  if($array_of_string[$i]=~/^((\S+)_(\d+)\-(\d+) *.*)$/){
		 unless(defined($old)){
			$old=$2;
			push(@temp1, $1);
			push(@temp2, $3);
		    next;
		 }elsif($2 eq $old){
			push(@temp1, $1);
			push(@temp2, $3);
			next;
		 }elsif( ($2 ne $old)||($i==$#array_of_string) ){
			&sort_and_put_strings_to_out;
		    push(@temp1, $1);
		    push(@temp2, $3);
			$old  =$2;

			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			sub sort_and_put_strings_to_out{
			   my ($j, $k, $num);
			   @temp2=sort { $a<=>$b } @temp2; ## sort numerically
			   F1: for($j=0; $j< @temp2; $j++){
				  $num=$temp2[$j];
				  for($k=0; $k< @temp1; $k++){
					 if($temp1[$k]=~/^(\S+)_$num\-/){
						push(@out, $temp1[$k]);
						splice(@temp1, $k, 1);
						$k--;
						splice(@temp2, $j, 1);
						$j--;
						next F1;
					 }
				  }
			   }
			   @temp1=@temp2=();

			 }#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	      }
	  }elsif($i > 0){ ## for the very last sort
		  &sort_and_put_strings_to_out;
	  }
   }
   return(\@out);
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
# Version   : 1.8
#--------------------------------------------------------------------
sub parse_arguments{
  my( $c, $d, $f, $arg_num, $option_table_seen, $n, $option_filtered,
		$option_table_example, $input_line, @input_files,
		$extension);
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
	my($table_found, @option_tb, $head);
	 open(SELF, "${$_[0]}");
	 while(<SELF>){
		if( (/^ *#+/) && ( $table_found== 1) ){
		  push (@option_tb, "$_");
		}elsif( ($table_found != 1)&&(/^ *\#+ *[Oo]ption *[Tt]able */) ){
			$table_found=1; $head="############## Option Table for $logname\'s \"$0\"\n"; ##
			push(@option_tb, $head);
		}
		if( ($table_found==1)&&(/^ *###################+ *$/)){  ### to find the end point of reading
			$table_found =0; last; }
	 }
	 return(\@option_tb);
}
#______________________________________________________________
# Title     : get_internal_dup_in_a_cluster
# Usage     :
# Function  :
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  :
# Options   : _  for debugging.
#             #  for debugging.
#  $short_region=  S by S -S  # taking shorter region overlapped in removing similar regions
#  $large_region=  L by L -L  # taking larger  region overlapped in removing similar regions
#  $average_region=A by A -A # taking average region overlapped in removing similar regions
#
# Version   : 1.2
#--------------------------------------------------------------
sub get_internal_dup_in_a_cluster{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

	$cluster_line=$_[0] || ${$_[0]};
	my(@seq, %out, $seq_name, $short_region, $large_region, $average_region);
	my $overlap_factor=30;
	my $min_inside_dom_size=30;
	@seq=split(/ +/, $cluster_line);  ## These sequence are single seq with different regions
	@seq= map{$_->[0]} sort{$a->[1] cmp $b->[1] or $a->[2] <=> $b->[2] }
			             map {/^((\S+)_(\d+)\-(\d+) *.*)$/ && [$1, $2, $3, $4]} @seq;
	if($char_opt=~/S/){       $short_region='S'; }
	if($char_opt=~/L/){    $large_region='L';   }
	if($char_opt=~/A/){    $average_region='A'; }

	F1:for($i=0; $i< @seq; $i++){
	   $seq1=$seq[$i];
	   if($seq1=~/^(\S+)_(\d+)\-(\d+)/){
		  $seq_name=$1;
		  $start1=$2;
		  $end1=$3;
	   }
	   F:for($j=1; $j< @seq; $j++){
		  $seq2=$seq[$j];
		  if($seq1 eq $seq2){ next } ### Skip IDENTICAL ones (xxxx_1-10, xxxx_1-10)
		  if($seq2=~/^(\S+)_(\d+)\-(\d+)/){
			 $start2=$2;
			 $end2=$3;
		  }
		  $leng2=$end2-$start2;
		  $margin=$leng2/12;   ## 8% overlap is regarded as not overlapping

		  if(( ($start1+$margin) > $end2)||
		    ( ($start2+$margin) > $end1)){ # skips non overlapping seqlets

			$out{"$start1\-$end1"}.="$start2\-$end2 ";

			splice(@seq, $j, 1);
			$j--;
		  }
	   }
	}
	#@out=sort (@out);
	#@out=@{&remove_dup_in_array(\@out)};
	#@out=@{&remove_similar_seqlets(\@temp, "f=2", $short_region, $large_region, $average_region)};
	return(\%out);
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
	  $title_entry_exist, $entry_value, $temp_W, $Warning_part
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
	return(\%Final_out);
}               ##<<--- ENd of the sub read_head_box

#______________________________________________________________
# Title     : check_linkage_of_2_similar_seqlet_sets
# Usage     :
# Function  : connects two clusters of seqlets if they share
#              identical or near identical seqlets
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  :
# Options   : _  for debugging.
#  $factor = by f=  # eg)  "f=$factor" in the higher level sub
#
# Returns   :
# Argument  :
# Version   : 1.7
#--------------------------------------------------------------
sub check_linkage_of_2_similar_seqlet_sets{
   my ($seq1, $name1, $start1, $end1, $seq2,
	   $leng1, $leng2, $name2, $start2, $end2, $diff_start,
	   $diff_end);
   my @splited1=@{$_[0]};
   my @splited2=@{$_[1]};

   my $link_or_not=0;
   my $factor=7;  # this means 70% sequence region overlap of the intermediate is chosen

   if($_[2]=~/f=(\S+)/i){
	  $factor=$1;
   }

   F1: for($s=0; $s<@splited1; $s++){
	  if($splited1[$s]=~/^ *((\S+)_(\d+)\-(\d+))/){
		  $seq1=$1;
		  $name1=$2;
		  $start1=$3;
		  $end1=$4;
	  }
	  F2: for($t=0; $t< @splited2; $t++){
		 if($splited2[$t]=~/^ *((\S+)_(\d+)\-(\d+))/){
			 $seq2=$1;
			 $name2=$2;
			 $start2=$3;
			 $end2=$4;
		 }
		 if($seq1 eq $seq2){ $link_or_not=1; return(\$link_or_not) }
		 if($name1 ne $name2){
			 next F2;
		 }elsif($name1 eq $name2){ ## ~~~~~~~~~~~~~ THIS is the MOST IMP CORE PART ~~~~~~~~~~~~~
			 $leng1=$end1-$start1;
		     $leng2=$end2-$start2;
			 if($leng1 >= $leng2){ $smaller_leng=$leng2; }else{ $smaller_leng=$leng1; }
			 $diff_start=abs($start1-$start2);
			 $diff_end  =abs($end1  -$end2  );
			 if((($diff_start+$diff_end)/2) <= ($smaller_leng/$factor) ){
			 	$link_or_not=1;
				return(\$link_or_not);
			 }
		 }## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  }
   }
   return(\$link_or_not);
}


#__________________________________________________________________________
# Title     : show_subclusterings
# Usage     : &show_subclusterings(\@out);
# Function  : This is the very final sub of divclus.pl
# Example   : @temp_show_sub=&show_subclusterings(\@out, $file, $sat_file, $dindom, $indup);
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  : print_subclusterings, sum_subclusterings, write_subclustering
#             show_clusterings, display_subclusterings
# Options   :
#             f  for file output, eg: xxxxxxx.sat
#
# Version   : 2.6
#-------------------------------------------------------------------------
sub show_subclusterings{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    my ($max_size, $sat_file_name, $clu_file_name,
		$ori_cluster_size, $ori_cluster_num, $good_bad, @keys, $percentage_fac,
		$indup, @sizes, $sum_seq_num, $indup_percent, $indup_count, %tem4,
		@sub_clustering_out_files);  # clusall_1e-5_clu_14-324_ss.sat
	my @out=@{$array[0]};
	$indup_count=0;

	if($char_opt=~/d/){	    $dindom=1;	}
	if($char_opt=~/i/){		$indup=1;	}
	if($vars{'f'}=~/\d+/){     $factor= $vars{'f'}; }
	if($vars{'p'}=~/\d+/){ $percentage_fac= int($vars{'p'}); }
	if($vars{'s'}=~/\d+/){	   $score = $vars{'s'};	}
	if($vars{'e'}=~/\d+/){	   $evalue= $vars{'e'};	}

	#print "\n# (1) show_subclusterings : \@file has : @file\n";
	if( $file[0]=~/([\S+_]*?(\d+)\-(\d+)[_\w]*)\.msp/  or
		$file[0]=~/([\S+_]*?(\d+)\-(\d+)[_\w]*)\.sat/   ){
		 $ori_cluster_size=$2;
		 $ori_cluster_num =$3;
		 $base=$1;
		 $sat_file_name="$base\.sat";
		 $clu_file_name="$base\.clu";
	}else{
         print "\n# (2) LINE:",__LINE__," The \@file input to show_subclusterings is not the right format, dying\n";
         print "\n     Sarah!, right format looks like: 13-234.msp or 8-420_cluster.msp \n";  exit;
	}

    open(CLU, ">$clu_file_name") or die "\n# (E) show_subclusterings failed miserably to open \"$clu_file_name\" \n";
    push(@sub_clustering_out_files, $clu_file_name);


	@out=@{&sort_string_by_length(\@out)};

	for($i=0; $i< @out; $i++){ # @out has ( 'YAL054C_98-695 YBR041W_90-617', 'YBR115C_230-842 YBR222C_16-537 YER015W_121-686', etc)
	   my $count+=$i+1;
	   my ( $int_dup_number, $sub_clu_size, $seq_with_range, @sp, $new_clus_NAME,
	        %tem, %tem2, %tem3, $j, @keys, $num_seq);
	   if($out[$i]=~/^ *$/){ next }
	   @sp=sort split(/ +/, $out[$i]);

	   for($j=0; $j < @sp; $j++){
		  $seq_with_range=$sp[$j];
		  if($seq_with_range=~/^((\S+)_((\d+)\-(\d+)))/){
			 $tem{$2}++;
			 $tem2{$2}.=sprintf("%-15s ", $1);
			 $tem3{$2} =$3;
			 $tem4{$2} =$5-$4;
		  }
	   }

	   @keys=sort keys %tem;
	   $num_seq=$sub_clu_size=@keys;

	   if($max_size < $sub_clu_size){
		  $max_size=$sub_clu_size; ## This is to collect the sizes of clusters to see if it is good.
	   }
	   $indup_count= &print_summary_for_divclus(
		         $count, \%tem2, \%tem,
		         $ori_cluster_num,
		         $ori_cluster_size,
		         $dindom,
		         $clu_file_name,
                 \%tem3, \%tem4,
                 $indup, );

       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # Local subroutine
       #_______________________________________________________________
	   sub print_summary_for_divclus{ #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			my $count=$_[0]; # count of cluster
			my %tem2=%{$_[1]};	my $num_seq=@keys=sort keys %tem2;
			my %tem=%{$_[2]};	my $ori_cluster_num=$_[3];
			my $new_clus_NAME=$ori_cluster_num.'0'.$count.'0'.$num_seq;
			my $ori_cluster_size=$_[4];
			my $dindom=$_[5];	my %tem3=%{$_[7]};
			my $indup=$_[9];	my (%internal_dup);
            my %tem4=%{$_[8]};

			#~~~~~~~~~~ Domain Inside Domain ~~~~~~~~~~~~~~~~~
			if($dindom){
				for($x=0; $x <@keys; $x++){
					@domain_inside_domain=@{&get_domain_inside_domain($tem2{$keys[$x]})};
					@domain_inside_domain=@{&remove_dup_in_array(\@domain_inside_domain)};
					for($m=0; $m< @domain_inside_domain; $m++){ print "  # Dindom: $m : $domain_inside_domain[$m]\n";   }
					print "\n";
				}
			}
			#==========================================================================================

			#~~~~~~~~~~ Internal duplication  ~~~~~~~~~~~~~~
			if($indup==1){
			   # @keys is the same as sub cluster size,
			   for($x=0; $x < @keys; $x++){
				  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				  # Checks each sequence for duplication
				  #___________________________________________________
				  my %internal_dup=%{&get_internal_dup_in_a_cluster( $tem2{$keys[$x]} )};
				  my @dup_keys=keys %internal_dup;
				  if(@dup_keys > 0){
					  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
					  #  This calculates the actual duplicated number rather than jus tthe sequences
					  #______________________________________________________________________________
					  $indup_count++;
					  printf ("%-14s %-12s %-4s", $keys[$x], $new_clus_NAME, $num_seq);
					  for($m=0; $m< @dup_keys; $m++){
						  printf ("%-19s=> %s\n", $dup_keys[$m], $internal_dup{ $dup_keys[$m] } );
					  }
				  }
			   }
			}

			#~~~~~~~~~~ Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~
			print  CLU  "Cluster size $num_seq\n";
            printf CLU ("Cluster number %-12s # E:%-5s Factor:%-2s P:%-2s, Ori size:%-4s Sub:%-4s From:%-12s\n",
						  $new_clus_NAME, $evalue, $factor, $percentage_fac,
						  $ori_cluster_size, $num_seq, $ori_cluster_num);
			print       "Cluster size $num_seq\n";
			printf     ("Cluster number %-12s # E:%-5s Factor:%-2s P:%-2s, Ori size:%-4s Sub:%-4s From:%-12s\n",
			              $new_clus_NAME, $evalue, $factor, $percentage_fac,
			              $ori_cluster_size, $num_seq, $ori_cluster_num);
			for($x=0; $x <@keys; $x++){
               printf CLU ("   %-4s %-5s %-17s %-10s %-3s leng: %-s\n",
                  $num_seq, $ori_cluster_num, $keys[$x], $tem3{$keys[$x]}, $tem{$keys[$x]}, $tem4{$keys[$x]});
			   printf     ("   %-4s %-5s %-17s %-10s %-3s leng: %-s\n",
				  $num_seq, $ori_cluster_num, $keys[$x], $tem3{$keys[$x]}, $tem{$keys[$x]}, $tem4{$keys[$x]});
			}
			return($indup_count);
	   }
	}
    close(CLU); ## this is a bug fix

	if($max_size == $ori_cluster_size){   $good_bad=1;
	}else{	                              $good_bad=0;	}

    print "\n# Sarah, Do you think the subclusterings are O.K.?" if $verbose;
    print "\n#   Tell me, if you feel suspicious, jong\@salts.med.harvard.edu\n\n" if $verbose;
    return($good_bad, $indup_count, $ori_cluster_size, \@sub_clustering_out_files);
}




#______________________________________________________________________
# Title     : sort_string_by_length (synonym of sort_str_by_length  )
# Usage     : @output = @{&sort_string_by_length(@any_input_strings, [-r], @more)};
# Function  : sorts strings in array according to their sizes
#             bigger comes first.
# Example   :
# Warning   :
# Keywords  : sort_array_by_length, sort_str_by_length, sort_array_string_by
#             sort_string_by_leng, sort_by_length, sort_by_leng,
# Options   : -r  reverse the order
# Version   : 1.2
#-------------------------------------------------------------------
sub sort_string_by_length{
	my(@input, $i, $small_first, @output);
	for($i=0; $i<@_; $i++){
		if( $_[$i]=~/^\-?r$/i){
			$small_first =1;
			splice(@_, $i, 1);
		}elsif(ref($_[$i]) eq 'ARRAY'){
		    push(@input, @{$_[$i]});
		}elsif(ref($_[$i]) eq 'SCALAR'){
			if(${$_[$i]}=~/^\-?r$/i){
			   $small_first=1;
			   splice(@_, $i, 1);
			}else{
			   push(@input, ${$_[$i]});
			}
		}elsif( !ref($_[$i]) ){
		    push(@input, $_[$i]);
		}
	}
	if($small_first ==1){
	    @output = sort {length($a) <=> length($a) || ($b cmp $a)} @input;
	}else{
	    @output = sort {length($b) <=> length($a) || ($a cmp $b)} @input;
	}
	return (\@output);
}


#______________________________________________________________
# Title     : cluster_merged_seqlet_sets
# Usage     : @out=@{&cluster_merged_seqlet_sets(\@lines)};
# Function  :
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  :
# Options   : _  for debugging.
#             #  for debugging.
#  $short_region=  S by S -S  # taking shorter region overlapped in removing similar regions
#  $large_region=  L by L -L  # taking larger  region overlapped in removing similar regions
#  $average_region=A by A -A # taking average region overlapped in removing similar regions
#
# Version   : 1.5
#--------------------------------------------------------------
sub cluster_merged_seqlet_sets{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   my ($optimize, @splited1, @splited2, $verbose, $link_or_not);
   my @lines=@{$array[0]};
   $link_or_not=0;
   my $factor=7; # 7 means 70% now

   if($vars{'f'}=~/(\d+)$/){ $factor=$1 }
   if($char_opt=~/o/){ $optimize=1 }
   if($char_opt=~/S/){ $short_region='S'; }
   if($char_opt=~/L/){ $large_region='L';   }
   if($char_opt=~/A/){ $average_region='A'; }
   if($char_opt=~/v/){ $verbose=1 }

   if($verbose){
       print "\n# (1) cluster_merged_seqlet_sets: Checking linkage and merging <<<<<>>>>>\n@lines\n";
   }

   F1: for($i=0; $i< @lines; $i++){
      @splited1=split(/ +/, $lines[$i]);
	  for($j=0; $j< @lines; $j++){
		 if($lines[$i] eq $lines[$j]){ next  }
		 @splited2=split(/ +/, $lines[$j]);

         $link_or_not=${&check_linkage_of_2_similar_seqlet_sets(\@splited1, \@splited2, "f=$factor")};
         print "\n +++++ \$link_or_not is  $link_or_not +++" if $verbose;
		 if($link_or_not==1){
             if($verbose){
                 print "\n# (2) cluster_merged_seqlet_sets: \n $lines[$i]  \n and $lines[$j] \n are linked \n";
             }

             if($optimize){ ##---- This will also remove similar seqlets, not only identical ones
                $lines[$i]=join(' ', sort @{&remove_similar_seqlets( [@splited1, @splited2],
                                            $short_region, $large_region, $average_region)} );
             }else{
                $lines[$i]=join(' ', sort @{&remove_dup_in_array( [@splited1, @splited2])} );
             }
             splice(@lines, $j,1);
             $j--; $i--;
             next F1;
		 }
	  }
   }
   return(\@lines);
}

#______________________________________________________________
# Title     : sort_words_in_string
# Usage     :
# Function  :
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Class     :
# Keywords  : sort_words_in_sequences, sort_sequences_in_string,
#             sort_strings_in_string,
# Options   : _  for debugging.
#             #  for debugging.
# Package   : Bio
# Reference : http://sonja.acad.cai.cam.ac.uk/perl_for_bio.html
# Returns   :
# Tips      :
# Argument  :
# Todo      :
# Author    : A Scientist
# Version   : 1.1
#--------------------------------------------------------------
sub sort_words_in_string{
   my @in=@{$_[0]} || @_;
   my (@OUT);
   for (@_){
	  push(@OUT, join(' ', sort split(/ +/) ));
   }
   return(\@OUT);
}

#________________________________________________________________________
# Title     : show_array
# Usage     : &show_array(\@input_array);
# Function  : for debugging purpose. Shows any array elem line by line.
# Example   : Output:      item1
#             Output:      item2
#             Output:      item3
# Warning   : can handle scalar ref, too.
# Keywords  :
# Options   : -h  for horizontal display of elements
#             c   for compact (do not put new line between array chunk)
#             s   for putting new line between arrays
# Returns   :
# Argument  :
# Version   : 2.4
#--------------------------------------------------------------------
sub show_array{
  my($k, $i, $t,  @in2, $in, $space, $show_horizontally, $compact);
  my(@in)=@_;

  ## This is to get the option of 'horizontal' to make horizontal output.
  for($t=0; $t < @in ; $t++){
	 if($in[$t] =~/\-?[hH][orizontal]*$/){   ### No ref.
		 $show_horizontally = "h";
		 splice(@in, $t, 1);  $t--;
	 }elsif(${in[$t]} =~/-?[hH][orizontal]*$/){  ### ref.
		 $show_horizontally = "h";
		 splice(@in, $t, 1);  $t--;
	 }elsif(${in[$t]} =~/^s$/i){  ### ref.
		 $space = "s";
		 $compact='';
		 splice(@in, $t, 1);  $t--;
	 }elsif(${in[$t]} =~/^c$/i){  ### ref.
		 $compact = "c";
		 $space='';
		 splice(@in, $t, 1);  $t--;
	 }
  }

  for($k=0; $k < @in; $k++){
	 if(ref($in[$k]) eq 'ARRAY'){
		 &show_array(@{$in[$k]}, "$show_horizontally", "$compact", "$space" );
	 }elsif(ref($in[$k]) eq 'SCALAR'){
		 if($show_horizontally eq "h"){
			 print ${$in[$k]}, ",  ";
		 }elsif(  $show_horizontally ne "h"){
			 print ${$in[$k]}, "\n";
		 }
	 }elsif( !ref($in[$k]) ){
		 if($show_horizontally eq 'h'){
			 print  $in[$k] , ",  ";
		 }elsif(  $show_horizontally ne "h"){
			 print  $in[$k] , "\n";
		 }
	 }
  }
  if($compact !~/^c$/i){
	print "\n"; #### This is necessary to distinguish different arrays.
  }
}


#_________________________________________________________________________
# Title     : merge_similar_seqlets
# Usage     : @all_seqlets = @{&merge_similar_seqlets(@all_seqlets)};
# Function  : merges seqlet sets which have identical
#             sequences and share similar regions by connection factor of 30%
#             This means, if any two seqlets from the same sequences which
#             share more than 70% seqlet regions overlapping are merged
#             This only sees the very first sequence in the seqlets line!!!
#             (so, PARTIAL MERGE !!)
# Example   : INPUT:
#
#   @input=( 'seq1_1-30 seq2_1-40 seq3_1-50',
#            'seq1_2-49 seq3_4-40 seq4_2-99'....)
#
#   @output=('seq1_1-30 seq2_1-45 seq3_2-45 seq4_2-99');
#
# Keywords  : merge_similar_sequences, merge_sequence_names, merge_sequences,
#              merge_sequence_ranges, merge_similar_sequences_with_ranges,
#              merge_seqlets, merge_duplication_modules
# Options   :
#
#   f=<digit>   for determing the factor in filtering out non-homologous
#                  regions, 7 = 70% now!!
#   l=<digit>   for seqlet(duplication module) length threshold
#   z           for activating remove_similar_sequences, rather than remove_dup....
#   S  $short_region=  S by S -S  # taking shorter region overlap in removing similar reg
#   L  $large_region=  L by L -L  # taking larger  region overlap in removing similar reg
#   A  $average_region=A by A -A  # taking average region overlap in removing similar reg
#
# Version   : 2.0
#-------------------------------------------------------------------------------
sub merge_similar_seqlets{
   my (@all_seqlets, @result_all_seqlets, $i, $j, $k, $seq1, $start1, $end1, $seq2,
	   $smaller_leng, $start2, $end2, @split, @split1, @split2,
       $short_region, $large_region, $average_region, $overlapping_seq_match_size);
   my $factor=7;     #  30% sequence mismatch region is allowed(3)
   my $leng_thresh=30;
   my $optimize=1;
   $average_region='A'; # default

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   # Sorting (parsing) input to get options and input array
   #_________________________________________________________
   for($i=0; $i< @_; $i++){
	   if(ref($_[$i]) eq 'ARRAY'){
           @all_seqlets=@{$_[$i]};
       }elsif($_[$i]=~/f=(\S+)/){  $factor=$1;
       }elsif($_[$i]=~/z/i){       $optimize=1;
       }elsif($_[$i]=~/l=(\d+)/i){ $leng_thresh=$1;
       }elsif($_[$i]=~/^S/){       $short_region='S';   $large_region=$average_region='';
       }elsif($_[$i]=~/^L/){       $large_region='L';   $short_region=$average_region='';
       }elsif($_[$i]=~/^A/){       $average_region='A'; $short_region=$large_region  =''; }
   }
   if(@all_seqlets==1){

       print "\n# (1) merge_similar_seqlets: \@all_seqlets == 1, returning\n" if $verbose;
       print "\n#  \@all_seqlets = @all_seqlets\n" if $verbose;
       return(\@all_seqlets);
   }else{
       print "\n# (1) merge_similar_seqlets: \@all_seqlets > 1, Processing\n" if $verbose;
       if( $verbose){
           for($i=0; $i< @all_seqlets; $i++){
               print "\n# $all_seqlets[$i]" ;
           }
       }
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # This is to remove which are identical in @all_seqlets;
   #_________________________________________________________
   F1: for($i=0; $i< @all_seqlets; $i++){
      my $merged_two_seqlet_lines;

      if($all_seqlets[$i] eq $all_seqlets[$i+1]){
	      print "\n# (2) merge_similar_seqlets: \$all_seqlets\[$i\] equals \$all_seqlets\[$i+1\]\n" if $verbose;
          splice(@all_seqlets, $i+1, 1);
          $i-- if $i >0;
		  next F1;
	  }else{
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # @split1 and 2 are arrays from different string entry in @all_seqlets
          #______________________________________________________________________
          @split1=sort split(/ +/, $all_seqlets[$i]);
          @split2=sort split(/ +/, $all_seqlets[$i+1]);
      }

	  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
	  #  (3) If the first elements of @split1 and 2 are identical, lets merge the two arrays
	  #________________________________________________________________________________
      if($split1[0] eq $split2[0] or $split1[0] eq $split2[1] or $split1[0] eq $split2[2]){
	      print "\n# (3) \$split1[0] and \$split2[0] are identical, good!\n" if $verbose;
		  @split=(@split1, @split2);
		  if(1){ #~~~~~ optimize option removes similar seqlets
             $all_seqlets[$i]= join(' ', sort @{&remove_similar_seqlets(\@split,
			                              $short_region, $large_region, $average_region)} );
		  }else{
			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 # Only removes exactly identical ones
			 #__________________________________________________________
             $all_seqlets[$i]=  join(' ', @{&remove_dup_in_array(\@split, 's')} );
		  }
          print "\n# (3) New \$all_seqlets\[\$i\] is \n$all_seqlets[$i]\n" if $verbose;
		  splice(@all_seqlets, $i+1, 1);
		  $i-- if $i >0;
		  next F1;
	  }

	  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # (4) If the first elements of @split1 and 2 are NOT identical, lets check the sequence ranges
	  #_____________________________________________________________________________________________
      if($verbose){
          print "\n# (4) merge_similar_seqlets: the first elements of \@split1 and \@split2 are NOT identical\n";
          print "# (4) >>>>>>>>> Checking the regions of them, \$factor= $factor\n";
      }

      F2: for($j=0; $j < @split1; $j++){
          if($split1[$j] =~/^ *(\S+)_(\d+)\-(\d+)/){
               my ($seq1, $start1, $end1)=($1, $2, $3);

               F3: for($k=0; $k<@split2; $k++){
                   if($split2[$k] =~/^ *(\S+)_(\d+)\-(\d+)/){
                       my($seq2, $start2, $end2)=($1, $2, $3);

                       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````
                       # Check if the seqs are identicl (from the two arrays), no point to merge which are not identical from the first
                       #__________________________________________________________________________________________
                       if($seq1 eq $seq2){
                            if($verbose){
                                print "# (5) <<<<< The sequence names are identical with different regions, See if we can merge!\n";
                            }
                            $diff_start=abs($start1-$start2);
                            $diff_end  =abs($end1  -$end2  );
                            $leng1=$end1-$start1;
                            $leng2=$end2-$start2;
                            if($leng1 >= $leng2){
                                $smaller_leng=$leng2;
                                $larger_leng =$leng1
                            }else{
                                $smaller_leng=$leng1;
                                $larger_leng =$leng2
                            }

                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            # Checking the minimal seq region leng here
                            #______________________________________________________
                            if($smaller_leng < $leng_thresh){ next }

                            $overlapping_seq_match_size=${&get_overlapping_seq_match_size($start1, $end1, $start2, $end2)};
                            $averge_seq_leng_of_2_seqs=($leng1+$leng2)/2;

                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            # This is the critically important part
                            #_______________________________________________________________
                            if($average_region){
                                $finally_adjusted_seq_leng=$averge_seq_leng_of_2_seqs*($factor/10);
                            }elsif($short_region){
                                $finally_adjusted_seq_leng=$smaller_leng*($factor/10);
                            }elsif($large_region){
                                $finally_adjusted_seq_leng=$larger_leng*($factor/10);
                            }

                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
                            # Now let's check if we regard them homologous or not\
                            #_______________________________________________________
                            if( $overlapping_seq_match_size >=  $finally_adjusted_seq_leng){

                                @split= (@split1, @split2);
                                if ($verbose){
                                    print "# (5) merge_similar_seqlets: \$averge_seq_leng_of_2_seqs=$averge_seq_leng_of_2_seqs";
                                    print " $$$$ Merging occurs with \n@split  \n";
                                }
                                if($optimize){ #~~~~~ $optimize option removes similar seqlets
                                    $all_seqlets[$i]= join(' ', sort @{&remove_similar_seqlets(\@split,
			                                  $short_region, $large_region, $average_region)} );
                                }else{
                                    $all_seqlets[$i]= join(' ', @{&remove_dup_in_array(\@split, 's')} );
                                }
                                $merged_two_seqlet_lines=1;
                                splice(@all_seqlets, $i+1, 1);
                                $i-- if $i >0;
                                next F1;

                            }
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
                            # We believe they are not homologous
                            #____________________________________________
                            else{
                                if ($verbose){
                                    print "\n# (5) merge_similar_seqlets !! MERGING DID NOT occur->
                                    $split1[$j] $split2[$k] \n";
                                    print "\n#   \$averge_seq_leng_of_2_seqs= $averge_seq_leng_of_2_seqs, \$overlapping_seq_match_size= $overlapping_seq_match_size \$finally_adjusted_seq_leng= $finally_adjusted_seq_leng\n";
                                    print "\n# (5) merge_similar_seqlets, \$all_seqlets[$i]\$all_seqlets[$i+1]\n$all_seqlets[$i]\n$all_seqlets[$i+1]\n";
                                }
                                next F3;
                            }
                       }
                   }
                   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   # If there is no range (region) in seq naem, let's skip, as there is no way to check
                   #__________________________________________________________________________________
                   else{ # when split2 does not match xxx_10-20 format
                       next;
                   }
               }
          }else{ # when split1 does not match xxx_10-20 format
               next;
          }
      }
      unless($merged_two_seqlet_lines){

      }
   }
   print "\n# \@all_seqlets is  @all_seqlets =======\n" if $verbose;
   return(\@all_seqlets);
}


#________________________________________________________________________________
# Title     : get_overlapping_seq_match_size
# Usage     : $ovlapsize=${&get_overlapping_seq_match_size($st1, $en1, $st2, $en2)
# Function  :
# Example   :
# Keywords  :
# Options   :
# Version   : 1.1
#--------------------------------------------------------------------------------
sub get_overlapping_seq_match_size{
    my($start1, $end1, $start2, $end2, $overlapping_region_matched);
    if(@_ == 4){
       $start1=$_[0]; $end1 =$_[1];  $start2=$_[2]; $end2  =$_[3];
    }elsif(@_==2){
       if( $_[0]=~/(\d+)\-(\d+)/ ){
           $start1=$1;      $end1  =$2;
       }elsif($_[1]=~/(\d+)\-(\d+)/ ){
           $start2=$1;      $end2  =$2;
       }else{
           print "\n# (E) get_overlapping_seq_match_size: I need 2 or 4 arguments for regions\n";
           print "   They look like ($start1, $end1, $start2, $end2) or ('10-100', '20-211')\n";
           print "   You got it, Sarah?? Try again my dear!\n";
       }
    }else{
           print "\n# (E) get_overlapping_seq_match_size: I need 2 or 4 arguments for regions\n";
           print "   They look like ($start1, $end1, $start2, $end2) or ('10-100', '20-211')\n";
           print "   You got it, Sarah?? Try again my dear!\n";
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #     ---------
    #  ------
    #___________________________________
    if($start1 >= $start2 and $end1 >= $end2){
        $overlapping_region_matched=$end2-$start1;
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ---------
    #     ----------
    #___________________________________
    elsif($start1 <= $start2 and $end1 <= $end2){
        $overlapping_region_matched=$end1-$start2;
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #      -----
    #    ----------
    #___________________________________
    elsif($start1 >= $start2 and $end1 <= $end2){
        $overlapping_region_matched=$end1-$start1;
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  ---------
    #    ----
    #___________________________________
    elsif($start1 <= $start2 and $end1 >= $end2){
        $overlapping_region_matched=$end2-$start2;
    }
    return(\$overlapping_region_matched);
}


#_____________________________________________________________________________
# Title     : remove_similar_seqlets
# Usage     : @seqlets=@{&remove_similar_seqlets(\@split)};
# Function  : merges(gets average starts and ends ) of similar
#             seqlets to reduce them into smaller numbers. This can also handle
#              names like XLBGLO2R_8-119_d1hlm__.
#
# Example   : @seqlets=@{&remove_similar_seqlets(\@mrg1, $mrg2, \@mrg3)};
#               while @mrg1=qw(M_2-100 M_2-110 M_8-105 M_4-108 N_10-110 N_12-115);
#                     $mrg2='Z_3-400 Z_2-420';
#                     @mrg3=('X_2-300 X_3-300', 'X_2-300', 'X_5-300 X_2-301' );
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  : merge_sequence_names, merge_seq_names, merge_sequence_ranges
#             merge_seq_ranges
# Options   : _  for debugging.
#             #  for debugging.
#             f= for factor
#             S  for shorter region matched is used
#             A  for average region matched is used
#             L  for larger region matched is used
#
# Version   : 2.0
#-------------------------------------------------------------------------------
sub remove_similar_seqlets{
   my ($i, $j, $seq1, $smaller_leng, $leng1, $leng2, $start1, $end1, $seq2, $start2,
	   $av_diff, $num_of_seq, $av_end, $av_start, $end2, @seqlets,
	   @array_input, @seqlet, $tail1, $tail2, $shorter_region, $larger_region,
	   $average_region);
   my $factor=7;  ## !!! This var makes big difference in the final clustering
   $average_region = 'A'; ## default is to get the average of comparing regions

   for($i=0; $i< @_; $i++){
	   if(ref($_[$i]) eq 'ARRAY'){
		   @array_input=@{$_[$i]};
		   for($j=0; $j<@array_input; $j++){
			   @seqlet=split(/ +/, $array_input[$j]);
               push(@seqlets, @seqlet);
		   }
		   if($verbose){
               print "\n# remove_similar_seqlets: ARRAY ref is given as input\n";
               print "#  They are: @seqlets, Sarah\n";
           }
	   }elsif($_[$i]=~/f=(\S+)/){ $factor=$1
	   }elsif($_[$i]=~/^(S) *$/){     $shorter_region=$1 ; $average_region=0;
	   }elsif($_[$i]=~/^(L) *$/){     $larger_region =$1 ; $average_region=0;
	   }elsif($_[$i]=~/^(A) *$/){     $average_region=$1 ; $shorter_region=$larger_region=0;
	   }elsif($_[$i]=~/\S+\_\d+\-\d+/){
		   push(@seqlets, split(/ +/, $_[$i]) );
	   }elsif(ref($_[$i]) eq 'SCALAR' and ${$_[$i]}=~/\S+\_\d+\-\d+/){
	       push(@seqlets, split(/ +/, ${$_[$i]}) );
	   }
   }
   print "\n# remove_similar_seqlets : I am using \$factor : $factor\n" if $verbose;

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Sorting is necessary as I am not doing the real thorough comparison
   #______________________________________________________________________
   $num_of_seq=@seqlets=sort @seqlets;
   if($verbose){
       print "\n# (1) remove_similar_seqlets: Num of seq to merge: $num_of_seq (from \@seqlets)";
       print "\n# (2) remove_similar_seqlets: \@seqlets are  @seqlets\n";
   }
   my ($short_start, $large_start, $short_end, $large_end);

   for($i=0; $i< @seqlets; $i++){
	  if($seqlets[$i]=~/^ *(\S+)_(\d+)\-(\d+)(\S*)/){  ## last (\S*) is necessary for XLBGLO2R_8-119_d1hlm__
		 my($seq1, $start1, $end1, $tail1)=($1, $2, $3, $4);
	     if($seqlets[$i+1]=~/^(\S+)_(\d+)\-(\d+)(\S*)/){
			($seq2, $start2, $end2, $tail2)=($1, $2, $3, $4);
			if($seq1 eq $seq2){
			   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			   $diff_start=abs($start1 - $start2);
			   $diff_end  =abs($end1   - $end2  );
			   $leng1=$end1-$start1;
		       $leng2=$end2-$start2;

			   if($leng1 >= $leng2){ $smaller_leng=$leng2; }else{ $smaller_leng=$leng1; }
			   if( ($diff_start+$diff_end)/2 <= $smaller_leng*($factor/10) ){

				   if($average_region){
					   $av_start=int(($start1+$start2) / 2);
					   $av_end  =int(($end1 + $end2) / 2);
                       $seqlets[$i]="$seq1\_$av_start\-${av_end}$tail1";  # $tail1 is for names like XLBGLO2R_8-119_d1hlm__
                       print "\n# new seqlet : $seqlets[$i]\n" if $verbose;
					   splice(@seqlets, $i+1, 1);
					   $i--;
				   }else{
					   if($start1 < $start2){
							$short_start=$start2; $large_start=$start1;  ## note that short start should be $start2 if $start2 is bigger
					   }else{
							$short_start=$start1; $large_start=$start2;
					   }
					   if($end1 < $end2){
							$short_end=$end1;  $large_end=$end2;
					   }else{
							$short_end=$end2;  $large_end=$end1;
					   }
					   if($shorter_region){
						   $seqlets[$i]="$seq1\_$short_start\-${short_end}$tail1";
					   }elsif($larger_region){
						   $seqlets[$i]="$seq1\_$large_start\-${large_end}$tail1";
                       }

					   splice(@seqlets, $i+1, 1);
					   $i--;
			       }
			   }
			}
		 }
	  }
   }
   print "\n# (3) remove_similar_seqlets: The final out are: @seqlets\n" if $verbose;
   return(\@seqlets);
}


#______________________________________________________________
# Title     : get_domain_inside_domain
# Usage     :
# Function  :
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  : find_dindoms, domain_inside_domain, domain_in_domain
# Options   : _  for debugging.
#             #  for debugging.
# Returns   :
# Argument  :
# Version   : 1.0
#--------------------------------------------------------------
sub get_domain_inside_domain{
	$cluster_line=$_[0] || ${$_[0]};
	my($i, $j, @seq, @out);
	my $overlap_factor=40;
	my $min_inside_dom_size=40;
	@seq=split(/ +/, $cluster_line);
	F1:for($i=0; $i< @seq; $i++){
	   $seq1=$seq[$i];
	   if($seq1=~/^(\S+)_(\d+)\-(\d+)/){
		  $seq_name=$1;
		  $start1=$2;
		  $end1=$3;
	   }
	   F:for($j=0; $j< @seq; $j++){
		  $seq2=$seq[$j];
		  if($seq1 eq $seq2){ next } ### Skip IDENTICAL ones (xxxx_1-10, xxxx_1-10)
		  if($seq2=~/^(\S+)_(\d+)\-(\d+)/){
			 $start2=$2;
			 $end2=$3;
		  }
		  if(($start1 > $end2)||($start2 > $end1)){ # skips non overlapping seqlets
			 next;
		  }
		  if(($start1 > $start2)&&($end1 < $end2)){  #   -----
			 $leng_seq1=$end1-$start1;               # ----------
			 $leng_seq2=$end2-$start2;
			 if(( ($leng_seq2/2) >= $leng_seq1 )&&
			    ($leng_seq1 > $min_inside_dom_size) ){   # if seq1 is less than 60% of seq2, it is a hidden domain
				push(@out, "$seq2\($seq1\)");
			 }
		  }elsif(($start1 < $start2)&&($end1 > $end2)){  # -----------
			 $leng_seq1=$end1-$start1;                   #   ------
			 $leng_seq2=$end2-$start2;
			 if(( ($leng_seq1/2) >= $leng_seq2)&&
			    ($leng_seq2 > $min_inside_dom_size) ){   # if seq1 is less than 60% of seq2, it is a hidden domain
				push(@out, "$seq1\($seq2\)");
			 }
		  }
	   }
	}
	return(\@out);
}
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
# Title     : handle_arguments
# Usage     : Just put the whole box delimited by the two '###..' lines below
#             to inside of your subroutines. It will call 'handle_arguments'
#             subroutine and parse all the given input arguments.
#             To use, claim the arguments, just use the variable in the box.
#             For example, if you had passed 2 file names for files existing
#             in your PWD(or if the string looks like this: xxxx.ext),
#             you can claim them by $file[0], $file[1] in
#             your subroutine.
# Function  : Sorts input arguments going into subroutines and returns default
#             arrays of references for various types (file, dir, hash, array,,,,)
#             If you give (\@out, @file), it will put @out into @array as a ref
#             and also the contents of @out will be dereferenced and put to
#             raw_string regardless what is in it).
#
# Example   : 'handle_arguments(\@array, $string, \%hash, 8, 'any_string')
# Warning   :
# Keywords  : handling arguments, parsing arguments,
# Options   :
# Returns   : Following GLOBAL variables
#
#             $num_opt,    @num_opt     @file          @dir
#             $char_opt,   @char_opt    %vars          @array,
#             @hash        @string,     @raw_string    @range,
#
#             $num_opt has 10,20
#             @num_opt has (10, 20)
#             @file has  xxxx.ext
#             @dir has  dir  or /my/dir
#             $char_opt has 'A,B'
#             @char_opt has (A, B)
#             @array has  (\@ar1, \@ar2)
#             @hash has (\%hash1, \%hash2)
#             @string  ('sdfasf', 'dfsf')
#             @raw_string (file.ext, dir_name, 'strings',,)
#             @range has values like  10-20
#             %vars deals with x=2, y=3 stuff.
#
# Argument  : any type, any amount
# Version   : 4.8
#--------------------------------------------------------------------
sub handle_arguments{
	my($c, $d, $e, $f, $i, $j, $k, $l, $s, $t, $x, $y, $z, $char_opt, $dir, @hash,
		$file, $in_dir, $num_opt, @char_opt, @dir, @file, @string, @file_dir, @k,
		@num_opt, @raw_string,@string, @array, %vars, @range, @temp, $temp,
		@char_options);

  &set_debug_option;
  if(@_<1){ print chr(7),"\n This is handle_arguments. No args Passed, Error?\n"}
  elsif( (@_ ==1)&& (ref($_[0]) eq 'ARRAY') ){ # when there is only 1 argument
	  push(@array, $_[0]);
	  push(@k, $_[0]);
  }elsif( (@_==1)&&( !ref($_[0]) ) ){
	  if(-f $_[0]){ push(@file, $_[0]);   push(@string, $_[0]) }
	  elsif(-d $_[0]){ push(@dir, $_[0]); push(@string, $_[0]) }
	  elsif($_[0]=~/^\d+$/){ push(@num_opt, $_[0]); $num_opt.=$_[0] }
	  elsif($_[0]=~/^\w+$/){ push(@string, $_[0]); }
  }elsif(@_ >=1){ @k = @_ }

  #####______Start of  general argument handling______######
  for($k=0; $k < @k ;$k++){
	  if( !ref($k[$k]) ){
		  if($k[$k]=~ /^[\-]?([a-zA-Z]\d*) {0,5}$/){  push(@char_opt, $1); $char_opt .= "$1\,";
		  }elsif($k[$k]=~ /^\-([a-zA-Z]+)$/){          ## When multiple option is given,
			  @char_options = split(/\,|/, $1);  push(@char_opt, @char_options);
			  $char_opt .= join("\,", @char_options); ## '-' should be used. eg. '-HEGI'
		  }elsif($k[$k]=~ /^(\w+)\=(\S* *)$/){  $vars{$1}=$2;  $vars .= "$1\,";
		  }elsif($k[$k]=~ /^(\-?\d+)$/){ push(@num_opt, $1);  $num_opt .= "$1\,";
		  }elsif($k[$k]=~ /^\d+\.?\d*\-\d+\.?\d*$/){  push(@range,  $k[$k] );
		  }elsif(-f $k[$k]){                          push(@file,   $k[$k] );
		  }elsif(-d $k[$k]){                          push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /\/[\w\d\.\-]+[\/].+[\/]$/){push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /^\/[\w\d\.\-]+[\/]*$/){    push(@dir,    $k[$k] );
		  }elsif($k[$k]=~ /^[\/\w\d\-\.]+\.\w+$/){    push(@file,   $k[$k] );
		  }elsif($k[$k]=~ /\S\/[\/\w\d\-\.]+\.\w+$/){ push(@file,   $k[$k] );
		  }elsif($k[$k]=~/^\w+[\/\\\w\d\.\-]+$/){     push(@string, $k[$k] );
		        # string does not have space, but includes '\', '/', '.'
		  }else{                                      push(@raw_string, $k[$k] );  }

	  }elsif( ref($k[$k]) ){
		  if( ref($k[$k]) eq "SCALAR"){
			 if(${$k[$k]} =~ /^[\-]?([a-zA-Z]\d*) {0,5}$/){ push(@char_opt, $1); $char_opt  .= "$1\,";
				}elsif(${$k[$k]}=~ /^\-([a-zA-Z]+)$/){ push(@char_opt, @char_options);
					$char_opt  .= join("\,", @char_options);  ## as an option string.
				}elsif(${$k[$k]}=~ /^(\w+)\=(\S* *)$/){  $vars{$1}=$2;  $vars .= "$1\,";
				}elsif(${$k[$k]}=~ /^(\-?\d+)$/){ $num_opt .= "$1\,";  push(@num_opt, $1);
			    }elsif(${$k[$k]}=~ /^\d+\.?\d*\-\d+\.?\d*$/){    push(@range,  $k[$k] );
				}elsif(-f ${$k[$k]}){                            push(@file,   ${$k[$k]} );
				}elsif(-d ${$k[$k]}){                            push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~ /\/[\/\w\d\.\-]+[\/].+[\/]/){ push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~/^\/[\/\w\d\.\-]+[\/]*$/){     push(@dir,    ${$k[$k]} );
				}elsif(${$k[$k]}=~ /^[\/\w\d\-\.]+\.\w+$/){      push(@file,   ${$k[$k]} );
				}elsif(${$k[$k]}=~/^\w+[\w\d\.\-]+$/){           push(@string, ${$k[$k]} );
				}else{                                           push(@raw_string, ${$k[$k]}); }
		  }elsif(ref($k[$k]) eq "ARRAY"){ my @temp_arr = @{$k[$k]}; push(@array, $k[$k]);
			for ($i=0; $i<@temp_arr; $i++){
			   if(-f $temp_arr[$i]){                            push(@file, $temp_arr[$i]);
			   }elsif($temp_arr[$i]=~/^\d+\.?\d*\-\d+\.?\d*$/){ push(@range,$temp_arr[$i] );
			   }elsif(-d $temp_arr[$i]){                        push(@dir , $temp_arr[$i]);
			   }elsif($temp_arr[$i]=~/\/[\/\w\d\.\-]+[\/].+[\/]/){ push(@dir, $temp_arr[$i] );
			   }elsif($temp_arr[$i]=~/^\/[\/\w\d\.\-]+[\/]*$/){ push(@dir, $temp_arr[$i] );
			   }elsif($temp_arr[$i]=~/^[\/\w\d\-\.]+\.\w+$/){   push(@file,$temp_arr[$i] );
																push(@string,$temp_arr[$i] );
			   }elsif($temp_arr[$i]=~/^\w+[\w\d\.\-]+$/){       push(@string,$temp_arr[$i]);
			   }else{                                           push(@raw_string, $temp_arr[$i]); }
			 }
		  }elsif(ref($k[$k]) eq "HASH"){                             push(@hash,   $k[$k] ); }
	  }
  }
  @raw_string=(@raw_string, @string);
  @file = @{&remove_dup_in_arrayH(\@file)};
  #-----------------------------------------------------
	 sub remove_dup_in_arrayH{  my($i, @nondup, @out_ref, %duplicate, @orig, @out_ref);
		for($i=0; $i<@_; $i++){  undef(%duplicate);
	       if(ref($_[$i]) eq 'ARRAY'){    @orig = @{$_[$i]};    }
		   @nondup = grep { ! $duplicate{$_}++ } @orig; push(@out_ref, \@nondup);  }
		if(@out_ref ==1){ return($out_ref[0]);}
		elsif(@out_ref >1){  return(@out_ref);}
	 }
  #-----------------------------------------------------
  return(\@hash, \@array, \@string, \@dir, \@file, \@num_opt,
			\@char_opt, \$num_opt, \$char_opt, \@raw_string, \%vars, \@range );
}

#______________________________________________________________
# Title     : convert_msp_line_to_mmp_line
# Usage     :
# Function  : this adds ranges to the seqnames of msp files
#             mmp line is msp line with additional sequences at the end
# Example   :
# Keywords  : convert_msp_to_mmp, convert_msp, convert_msp_2_mmp
#             change_msp_to_mmp, add_range_in_msp
# Options   : _  for debugging.
#             #  for debugging.
# Returns   :
# Argument  :
# Version   : 1.5
#--------------------------------------------------------------
sub convert_msp_line_to_mmp_line{
   my $input_msp=${$_[0]} || $_[0];
   my($score, $evalue, $long_1, $new_seq1, $new_seq2, $middle,
	  $start1, $end1, $start2, $end2, $seq1, $seq2, $new);

   if($input_msp=~/^ *(\d+) +(\S+) *\S*[ \t]+(\d+)[ \t]+(\d+)[ \t]+(\S+)[ \t]+(\d+)[ \t]+(\d+)[ \t]+(\S+)/){
	  ($score, $evalue, $start1, $end1, $start2, $end2)=($1, $2, $3, $4, $6, $7);
	  ($seq1, $seq2)=($5, $8);
	  if($seq1=~/(\S+)\_\d+\-\d+/){
		 $new_seq1="$1\_$start1\-$end1";
	  }else{
		 $new_seq1="$seq1\_$start1\-$end1";
	  }
	  if($seq2=~/(\S+)\_\d+\-\d+/){
		 $new_seq2="$1\_$start2\-$end2";
	  }else{
		 $new_seq2="$seq2\_$start2\-$end2";
	  }
	  $new=sprintf("%-6s %-8s %-5s %-5s %-32s %-5s %-5s %-32s",
					$score, $evalue, $start1, $end1, $new_seq1, $start2, $end2, $new_seq2);
   }
   return(\$new);
}
#________________________________________________________________________
# Title     : remove_dup_in_array
# Usage     : @out2 = @{&remove_dup_in_array(\@input1, \@input2,,,,)};
#             @out1 = &remove_dup_in_array(\@input1 );
# Function  : removes duplicate entries in an array. You can sort the
#             result if you wish by 's' opt. Otherwise, result will keep
#             the original order
# Example   : (1,1,1,1,3,3,3,3,4,4,4,3,3,4,4);  --> (1,3,4);
# Warning   :
# Keywords  : merge array elements, remove_repeting_elements,
#             remove_same_array_elements
# Options   :
#   s for sorting the array output
# Returns   : one or more references.
# Argument  : one or more refs for arrays or one array.
# Version   : 1.4
#--------------------------------------------------------------------
sub remove_dup_in_array{
  my($i, $sort_opt, @out_ref, @nondup,%duplicate, @orig, @out_ref);
  my @in=@_;
  for($i=0; $i<@in; $i++){
	 if($in[$i] eq 's'){
		$sort_opt=1;  splice(@in, $i, 1); $i--;
	 }elsif( ref($in[$i]) eq 'SCALAR'  and  ${$in[$i]} eq 's' ){
		$sort_opt=1;  splice(@in, $i, 1); $i--;
	 }
  }
  for($i=0; $i<@in; $i++){
	  undef(%duplicate);
	  if(ref($in[$i]) eq 'ARRAY'){    @orig = @{$in[$i]};    }
	  else{ @orig=@in }
	  @nondup = grep { ! $duplicate{$_}++ } @orig;    ## NOTE -> $_
	  if($sort_opt==1){ @nondup= sort @nondup }
	  push(@out_ref, \@nondup);
  }
  if(@out_ref ==1){ return($out_ref[0]);}
  elsif(@out_ref >1){  return(@out_ref);}
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
# Version   : 2.5
#--------------------------------------------------------------------
sub assign_options_to_variables{
  my($i, $j, $op, $z, $n, $symb, $value, $var, %val, @val, $option_table_example, @input_options);

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


#__________________________________________________________________________
# Title     : merge_sequence_in_msp_chunk
# Usage     :
# Function  : merges sequences which are linked by common regions
#             This filters the sequences by evalue and ssearch score
#             This is the main algorithm of merging similar sequences.
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  : connect_sequence_in_msp, link_sequence_in_msp_chunk
#             connect_sequence_in_msp_chunk, link_sequence_in_msp
#             merge_sequence, link_sequence, connect_sequence
# Options   : _  for debugging.
#             #  for debugging.
#             m  for merge file output format (.mrg)
#             t= for threshold of seqlet length eg)  "t=30"
#             f= for overlap factor (usually between 2 to 7 )
#                 2 means, if the two regions are not overlapped
#                  by more than HALF of of the smaller region
#                  it will not regard as common seqlet block
#             s= for ssearch score minimum
#             e= for ssearch e value maximum
#             S  for S -S  # taking shorter region overlapped in removing similar regions
#             L  for L -L  # taking larger  region overlapped in removing similar regions
#             A  for A -A # taking average region overlapped in removing similar regions
#
# Returns   :
# Argument  :
# Version   : 2.3
#--------------------------------------------------------------
sub merge_sequence_in_msp_chunk{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   my ($ssearch_score2, $evalue_found2, $evalue_found1, $ssearch_score1, $optimize );
   my ($L, %out_hash, @out, $LL, @Final_out, $verbose, $final_factor, $R_diff,
		$short_region, $large_region, $average_region);
   my $factor =4; # default factor for around 30% sequence mis-overlap is the threshold for common block
	  #~~~~~~~~~~~~~~ The lower the factor the larger clustering will occur ~~~~~~~~~~~~
   my $score  =75; # default ssearch score. seq below this will be chucked out
   my $evalue =40; # default maximum e value used. Seq higher than this will be thrown out
   my $thresh =40; # sequence length threshold. overlap less than this will be ignored

   if($char_opt=~/v/){     $verbose = 'v'
   }if($char_opt=~/z/){    $optimize = 'z'
   }if($char_opt=~/S/){    $short_region='S';
   }if($char_opt=~/L/){	   $large_region='L';
   }if($char_opt=~/A/){	   $average_region='A'; }

   if($vars{'t'}=~/\d+/){
	  $thresh=$vars{'t'}; print "\n# merge_sequence_in_msp_chunk: Thresh is $thresh\n" if (defined $verbose);
   }if($vars{'f'}=~/\d+/){
	  $factor=$vars{'f'}; print "\n# merge_sequence_in_msp_chunk: Factor is $factor\n" if (defined $verbose);
   }if($vars{'s'}=~/\d+/){
	  $score = $vars{'s'}; print "\n# merge_sequence_in_msp_chunk: Score is $score\n" if (defined $verbose);
   }if($vars{'e'}=~/\d+/){
	  $evalue= $vars{'e'}; print "\n# merge_sequence_in_msp_chunk: Evalue is $evalue\n" if (defined $verbose);
   }
   my @seqlets=split(/\n+/, (${$_[0]} || $_[0]) );
   print "@seqlets" if $verbose;

   F1: for($i=0; $i < @seqlets; $i ++){
	  if($seqlets[$i]=~/^ *((\d+) +(\d+\.?[e\-\d]*) +(\d+) +(\d+) +(\S+) +(\d+) +(\d+)) +(\S+) *(.*)/){
		  if($6 eq $9){ splice(@seqlets, $i, 1); $i--; next };
		  ($long_match1, $enq_seq1, $mat_seq1, $R_start1, $R_end1 )=($1, $6, $9, $4, $5);
		  $R_leng1=$R_end1-$R_start1;  $ssearch_score1= $2;  $evalue_found1 = $3;
	  }
	  if( ($R_leng1 < $thresh) || ($ssearch_score1 < $score) ){ splice(@seqlets, $i, 1); $i--; next; }
	  if( $evalue_found1 > $evalue){ splice(@seqlets, $i, 1); $i--; next; }

	  F2: for($j=0; $j < @seqlets; $j ++){
		 if($seqlets[$i] eq $seqlets[$j]){ next };
		 if($seqlets[$j]=~/^ *((\d+) +(\d+\.?[e\-\d]*) +(\d+) +(\d+) +(\S+) +(\d+) +(\d+)) +(\S+) *(.*)/){
			($long_match2, $enq_seq2, $mat_seq2, $R_start2, $R_end2)=($1, $6, $9, $4, $5);
			$R_leng2=$R_end2-$R_start2;			$ssearch_score2=$2;			$evalue_found2= $3;
	     }
		 if( ($R_leng2 < $thresh)||($ssearch_score2 < $score) ){ splice(@seqlets, $j, 1); $j--; next; }
		 if( $evalue_found2 > $evalue){ splice(@seqlets, $j, 1); $j--; next; }

		 $R_diff=abs($R_leng1-$R_leng2)/2;   ## <<<---- Note it is div by 2

		 if($R_leng2 < $R_leng1){ $smaller_leng=$R_leng2; }else{ $smaller_leng=$R_leng1; }

		 $Start_diff=abs($R_start1-$R_start2)/2; ## <<<---- Note it is div by 2
		 $final_factor=$smaller_leng/$factor;


		 #~~~~~~~~~~ If average R_diff and average Start_diff are less then 1/7 of the smaller seqlet
		 #~~~~~~~~~~ we regard they are same selqets
		 if(( $R_diff < $final_factor ) &&       ### $Start_diff is essential!
			($Start_diff < $final_factor ) ){  ### if diff is less than around 30% of the smaller length
			if($verbose=~/v/){
			   print "\n\$R_diff:$R_diff \$Start_diff:$Start_diff $smaller_leng $final_factor $factor";
			}
			if($R_leng2 >= $R_leng1){
			       #~~~~~ $mat_seq1 or $mat_seq2 can increase to 'slr1453,sll0238', so you need ',' in the middle only
				   $extended_name="$mat_seq2,$mat_seq1";
				   $L=length($extended_name);
				   $LL=length($long_match2)+2;
				   $seqlets[$i]= sprintf("%-${LL}s %-${L}s", $long_match2, $extended_name);
				   splice(@seqlets, $j, 1);
				   $i-- unless($i==0);
				   $j--;
				   next F1;
			}elsif( $R_leng1 >= $R_leng2){  ## chooses the bigger range seq
				   $extended_name="$mat_seq1,$mat_seq2"; # must be ',' not ' '
				   $L=length($extended_name);
				   $LL=length($long_match1)+2;
				   $seqlets[$i]=sprintf("%-${LL}s %-${L}s", $long_match1, $extended_name);
				   splice(@seqlets, $j, 1);
                   $i-- unless($i <= 0);
				   $j--;
				   next F1;
			}
	     }else{
			next F2;
		 }
	  }
   }
   if($char_opt=~/m/){
	  for($i=0; $i< @seqlets; $i++){
		if($seqlets[$i]=~/^ *\d+ +\d+\.?[e\-\d]* +\d+ +\d+ +(\S+) +\d+ +\d+ +(\S+) *$/){
		   if($1 eq $2){ next }
		   $leading_seq=$1; $long=$2; $long=~s/\,/ /g;
		   push(@Final_out, "$leading_seq $long" );
		}
	  }
   }
   @Final_out=sort @Final_out;
   print "\n\n\n# \@Final_out\n@Final_out \n=================\n " if $verbose;
   return(\@Final_out);
}

#______________________________________________________________
# Title     : sort_words_in_string
# Usage     :
# Function  : sort words in strings sperated by ' ' or "\n"
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  : sort_words_in_sequences, sort_sequences_in_string,
#             sort_strings_in_string,
# Options   : _  for debugging.
#             #  for debugging.
# Returns   :
# Argument  :
# Version   : 1.0
#--------------------------------------------------------------
sub sort_words_in_string{
   my @in=@{$_[0]} || @_;
   my @OUT;
   for (@_){
	  push(@OUT, join(' ', sort split(/ +|\n/) ));
   }
   return(\@OUT);
}



#________________________________________________________________________
# Title     : get_base_names
# Usage     : $base =${&get_base_names(\$file_name)};
#             :   or @bases = &get_base_names(\@files);  # <-- uses `pwd` for abs directory
# Function  : produces the file base name(eg, "evalign"  out of "evalign.pl" ).
# Example   : $base => 'test'  with 'test.txt' or '/home/dir/of/mine/text.txt'
# Warning   :
# Keywords  : get_base_name{, base_name, file_base_name ,  get_file_base_name
#             get_basename, basename, get_root_name
# Options   :
# Returns   :
# Argument  : handles both ref and non-ref.
# Version   : 1.3
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
			$pos = rindex($file_only, ".");
	        $base= substr($file_only, 0, $pos);
		}else{
			$file = $file[$x];
			$pos1=rindex($file, "/");
	        $file_only=substr($file, ($pos1+1));
			$pos = rindex($file_only, ".");
	        $base= substr($file_only, 0, $pos);
		}
		push(@base, $base);
	}
	if(@base == 1 ){ \$base[0] }else{ \@base }
}

__END__

#______________________________________________________________
# Title     : merge_similar_seqlets
# Usage     : @all_seqlets = @{&merge_similar_seqlets(@all_seqlets)};
# Function  : merges seqlet sets which have identical
#             sequences and share similar regions by connection factor of 30%
#             This means, if any two seqlets from the same sequences which
#             share more than 70% seqlet regions overlapping are merged
#             This only sees the very first sequence in the seqlets line!!!
#             (so, PARTIAL MERGE !!)
# Example   : INPUT:
#
#   @input=( 'seq1_1-30 seq2_1-40 seq3_1-50',
#            'seq1_2-49 seq4_4-40 seq8_2-99'....)
#
# Keywords  : merge_similar_sequences, merge_sequence_names,
#              merge_sequence_ranges, merge_similar_sequences_with_ranges
# Options   : _  for debugging.
#             #  for debugging.
#  $short_region=  S by S -S  # taking shorter region overlapped in removing similar regions
#  $large_region=  L by L -L  # taking larger  region overlapped in removing similar regions
#  $average_region=A by A -A # taking average region overlapped in removing similar regions
#
# Version   : 1.7
#--------------------------------------------------------------
sub merge_similar_seqlets{
   my (@all_seqlets, @result_all_seqlets, $i, $seq1, $start1, $end1, $seq2,
	   $smaller_leng, $start2, $end2, @split, @split1, @split2,
	   $short_region, $large_region, $average_region);
   my $factor=6.5;     #  33% sequence mismatch region is allowed(3)
   my $leng_thresh=30;
   my $optimize=0;
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   # Sorting (parsing) input to get options and input array
   #_________________________________________________________
   for($i=0; $i< @_; $i++){
	   if(ref($_[$i]) eq 'ARRAY'){
		   @all_seqlets=@{$_[$i]};
	   }elsif($_[$i]=~/f=(\S+)/){ $factor=$1
	   }elsif($_[$i]=~/o/i){      $optimize=1
	   }elsif($_[$i]=~/^S/){      $short_region='S';
	   }elsif($_[$i]=~/^L/){      $large_region='L';
	   }elsif($_[$i]=~/^A/){      $average_region='A'; }
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # This is to remove which are identical in @all_seqlets;
   #_________________________________________________________
   for($i=0; $i< @all_seqlets; $i++){
	  if($all_seqlets[$i] eq $all_seqlets[$i+1]){
		  push(@result_all_seqlets, $all_seqlets[$i]);
		  $i++;
		  next;
	  }
	  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  # @split1 and 2 are arrays from different string entry in @all_seqlets
	  #_________________________________________________________
	  @split1=sort split(/ +/, $all_seqlets[$i]);
	  @split2=sort split(/ +/, $all_seqlets[$i+1]);

	  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
	  #  (1) If the first elements of @split1 and 2 are identical, lets merge the two arrays
	  #________________________________________________________________________________
	  if($split1[0] eq $split2[0]){
		  @split=(@split1, @split2);
		  if($optimize==1){ #~~~~~ optimize option removes similar seqlets
			 push(@result_all_seqlets, join(' ', sort @{&remove_similar_seqlets(\@split,
			                              $short_region, $large_region, $average_region)} ));
		  }else{
			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 # Only removes exactly identical ones
			 #__________________________________________________________
			 push(@result_all_seqlets, join(' ', @{&remove_dup_in_array(\@split, 's')} ));
		  }
		  $i++;
		  next;
	  }
	  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
	  # (2) If the first elements of @split1 and 2 are NOT identical, lets check the sequence ranges
	  #________________________________________________________________________________
	  if($split1[0] =~/^(\S+)_(\d+)\-(\d+)/){
		   ($seq1, $start1, $end1)=($1, $2, $3);
		   if($split2[0] =~/^(\S+)_(\d+)\-(\d+)/){
			   ($seq2, $start2, $end2)=($1, $2, $3);

			   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~````
			   # Check if the seqs are identicl (from the two arrays), no point to merge which are not identical from the first
			   #__________________________________________________________________________________________
			   if($seq1 eq $seq2){
					$diff_start=abs($start1-$start2);
					$diff_end  =abs($end1  -$end2  );
					$leng1=$end1-$start1;
					$leng2=$end2-$start2;
					if($leng1 >= $leng2){ $smaller_leng=$leng2; }else{ $smaller_leng=$leng1; }

					#~~~~~~ If the sum of overhangs are smaller than a third of average length
					if( ( ($diff_start+$diff_end)/2 <= $smaller_leng/$factor ) &&
						($smaller_leng > $leng_thresh ) ){
						@split=(@split1, @split2);
						if($optimize==1){ #~~~~~ optimize option removes similar seqlets
						   push(@result_all_seqlets, join(' ', sort @{&remove_similar_seqlets(\@split,
						                            $short_region, $large_region, $average_region )} ));
						}else{
						   push(@result_all_seqlets, join(' ', @{&remove_dup_in_array(\@split, 's')} ));
						}
						$i++;
						next;
					}else{
						push(@result_all_seqlets, join(' ', @split1));
						next;
					}
			   }
			   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			   # As they are not teh same, lets just check the next one in @split2
			   #_____________________________________________________________________
			   else{
					push(@result_all_seqlets, join(' ', @split1));
					next;
			   }
		   }
		   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		   # If there is no range (region) in seq naem, let's skip, as there is no way to check
		   #__________________________________________________________________________________
		   else{
			   push(@result_all_seqlets, join(' ', @split1));
			   next;
		   }
	  }
   }
   return(\@result_all_seqlets);
}

#______________________________________________________________
# Title     : merge_sequence_in_msp_chunk
# Usage     :
# Function  : merges sequences which are linked by common regions
#             This filters the sequences by evalue and ssearch score
#             This is the main algorithm of merging similar sequences.
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  : connect_sequence_in_msp, link_sequence_in_msp_chunk
#             connect_sequence_in_msp_chunk, link_sequence_in_msp
#             merge_sequence, link_sequence, connect_sequence
# Options   : _  for debugging.
#             #  for debugging.
#             m  for merge file output format (.mrg)
#             t= for threshold of seqlet length eg)  "t=30"
#             f= for overlap factor (usually between 2 to 7 )
#                 2 means, if the two regions are not overlapped
#                  by more than HALF of of the smaller region
#                  it will not regard as common seqlet block
#             s= for ssearch score minimum
#             e= for ssearch e value maximum
#             S  for S -S  # taking shorter region overlapped in removing similar regions
#             L  for L -L  # taking larger  region overlapped in removing similar regions
#             A  for A -A # taking average region overlapped in removing similar regions
#
# Returns   :
# Argument  :
# Version   : 2.2
#--------------------------------------------------------------
sub merge_sequence_in_msp_chunk{
	#"""""""""""""""""< handle_arguments{ head Ver 4.1 >"""""""""""""""""""
	my(@A)=&handle_arguments(@_);my($num_opt)=${$A[7]};my($char_opt)=${$A[8]};
	my(@hash)=@{$A[0]};my(@file)=@{$A[4]};my(@dir)=@{$A[3]};my(@array)=@{$A[1]};
	my(@string)=@{$A[2]};my(@num_opt)=@{$A[5]};my(@char_opt)=@{$A[6]};
	my(@raw_string)=@{$A[9]};my(%vars)=%{$A[10]};my(@range)=@{$A[11]};
	my($i,$j,$c,$d,$e,$f,$g,$h,$k,$l,$m,$n,$o,$p,$q,$r,$s,$t,$u,$v,$w,$x,$y,$z);
	if($debug==1){print "\n\t\@hash=\"@hash\"
	\@raw_string=\"@raw_string\"\n\t\@array=\"@array\"\n\t\@num_opt=\"@num_opt\"
	\@char_opt=\"@char_opt\"\n\t\@file=\"@file\"\n\t\@string=\"@string\"\n" }
	#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   my ($ssearch_score2, $evalue_found2, $evalue_found1, $ssearch_score1, $optimize );
   my ($L, %out_hash, @out, $LL, @Final_out, $verbose, $final_factor, $R_diff,
		$short_region, $large_region, $average_region);
   my $factor =6; # default factor for around 30% sequence mis-overlap is the threshold for common block
	  #~~~~~~~~~~~~~~ The lower the factor the larger clustering will occur ~~~~~~~~~~~~
   my $score  =105; # default ssearch score. seq below this will be chucked out
   my $evalue =40; # default maximum e value used. Seq higher than this will be thrown out
   my $thresh =40; # sequence length threshold. overlap less than this will be ignored

   if($char_opt=~/v/){     $verbose = 'v'
   }if($char_opt=~/o/){    $optimize = 'o'
   }if($char_opt=~/S/){    $short_region='S';
   }if($char_opt=~/L/){	   $large_region='L';
   }if($char_opt=~/A/){	   $average_region='A'; }

   if($vars{'t'}=~/\d+/){
	  $thresh=$vars{'t'}; print "\n# merge_sequence_in_msp_chunk: Thresh is $thresh\n" if (defined $verbose);
   }if($vars{'f'}=~/\d+/){
	  $factor=$vars{'f'}; print "\n# merge_sequence_in_msp_chunk: Factor is $factor\n" if (defined $verbose);
   }if($vars{'s'}=~/\d+/){
	  $score = $vars{'s'}; print "\n# merge_sequence_in_msp_chunk: Score is $score\n" if (defined $verbose);
   }if($vars{'e'}=~/\d+/){
	  $evalue= $vars{'e'}; print "\n# merge_sequence_in_msp_chunk: Evalue is $evalue\n" if (defined $verbose);
   }
   my @seqlets=split(/\n+/, (${$_[0]} || $_[0]) );

   F1: for($i=0; $i < @seqlets; $i ++){
	  if($seqlets[$i]=~/^ *((\d+) +(\d+\.?[e\-\d]*) +(\d+) +(\d+) +(\S+) +(\d+) +(\d+)) +(\S+) *(.*)/){
		  if($6 eq $9){ splice(@seqlets, $i, 1); $i--; next };
		  ($long_match1, $enq_seq1, $mat_seq1, $R_start1, $R_end1 )=($1, $6, $9, $4, $5);
		  $R_leng1=$R_end1-$R_start1;
		  $ssearch_score1= $2;
		  $evalue_found1 = $3;
	  }
	  if( ($R_leng1 < $thresh) || ($ssearch_score1 < $score) ){ splice(@seqlets, $i, 1); $i--; next; }
	  if( $evalue_found1 > $evalue){ splice(@seqlets, $i, 1); $i--; next; }

	  F2: for($j=0; $j < @seqlets; $j ++){
		 if($seqlets[$i] eq $seqlets[$j]){ next };
		 if($seqlets[$j]=~/^ *((\d+) +(\d+\.?[e\-\d]*) +(\d+) +(\d+) +(\S+) +(\d+) +(\d+)) +(\S+) *(.*)/){
			($long_match2, $enq_seq2, $mat_seq2, $R_start2, $R_end2)=($1, $6, $9, $4, $5);
			$R_leng2=$R_end2-$R_start2;
			$ssearch_score2=$2;
			$evalue_found2= $3;
	     }
		 if( ($R_leng2 < $thresh)||($ssearch_score2 < $score) ){ splice(@seqlets, $j, 1); $j--; next; }
		 if( $evalue_found2 > $evalue){ splice(@seqlets, $j, 1); $j--; next; }

		 $R_diff=abs($R_leng1-$R_leng2)/2;   ## <<<---- Note it is div by 2

		 if($R_leng2 < $R_leng1){ $smaller_leng=$R_leng2; }else{ $smaller_leng=$R_leng1; }

		 $Start_diff=abs($R_start1-$R_start2)/2; ## <<<---- Note it is div by 2
		 $final_factor=$smaller_leng/$factor;


		 #~~~~~~~~~~ If average R_diff and average Start_diff are less then 1/7 of the smaller seqlet
		 #~~~~~~~~~~ we regard they are same selqets
		 if(( $R_diff < $final_factor ) &&       ### $Start_diff is essential!
			($Start_diff < $final_factor ) ){  ### if diff is less than around 30% of the smaller length
			if($verbose=~/v/){
			   print "\n\$R_diff:$R_diff \$Start_diff:$Start_diff $smaller_leng $final_factor $factor";
			}
			if($R_leng2 >= $R_leng1){
			       #~~~~~ $mat_seq1 or $mat_seq2 can increase to 'slr1453,sll0238', so you need ',' in the middle only
				   $extended_name="$mat_seq2,$mat_seq1";
				   $L=length($extended_name);
				   $LL=length($long_match2)+2;
				   $seqlets[$i]= sprintf("%-${LL}s %-${L}s", $long_match2, $extended_name);
				   splice(@seqlets, $j, 1);
				   $i-- unless($i==0);
				   $j--;
				   next F1;
			}elsif( $R_leng1 >= $R_leng2){  ## chooses the bigger range seq
				   $extended_name="$mat_seq1,$mat_seq2"; # must be ',' not ' '
				   $L=length($extended_name);
				   $LL=length($long_match1)+2;
				   $seqlets[$i]=sprintf("%-${LL}s %-${L}s", $long_match1, $extended_name);
				   splice(@seqlets, $j, 1);
				   $i-- unless($i <= 0);
				   $j--;
				   next F1;
			}
	     }else{
			next F2;
		 }
	  }
   }
   if($char_opt=~/m/){
	  for($i=0; $i< @seqlets; $i++){
		if($seqlets[$i]=~/^ *\d+ +\d+\.?[e\-\d]* +\d+ +\d+ +(\S+) +\d+ +\d+ +(\S+) *$/){
		   if($1 eq $2){ next }
		   $leading_seq=$1; $long=$2; $long=~s/\,/ /g;
		   push(@Final_out, "$leading_seq $long" );
		}
	  }
   }
   sort @Final_out;
   return(\@Final_out);
}



