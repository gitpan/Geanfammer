#!/usr/bin/perl
# Last Update by /usr/users/jpark/Perl/update_subroutines.pl: Sat Jan 24 20:29:50 EST 1998
#________________________________________________________________________
# Title     : geanfammer.pl
#
# Usage     : geanfammer.pl DATABASE(or GENOME) [e= ] [f=]
#                  * look at the Example section down below!
#
# Function  : Creates a domain level clustering file from a given
#              FASTA format sequence DB. It has been used for complete
#              genome sequence analysis.
#
#              ------------ USAGE INFORMATION -------------------
#             The parameters you put are important for the meaningful
#               protein family maker.
#             The most important one is the E and e options (Mostly,
#               they will have same value).
#             Large E is for setting the threshold for the single
#               linkage clustering.
#             This means, any sequence hit BELOW the threshold
#               (which is good ) will be linked.
#             For example, if Seq1 matched with Seq2 with E value
#              of FASTA search:
#              0.001, and you set the threshold 0.1, then YOU
#              ordered the geanfammer to regard them a family.
#
#             The second small e option is for the dividing a complex
#              and wrong cluster into correct more correct
#              duplication modules. This is necessary as a
#              lot of multidomain proteins can be clustered together
#              WRONGLY by single linkage.
#             At this stage, the e value is irrelevant to E value
#              and you can set a higher or lower one. Or you can set
#              the same as E(just set the 2 the same!)
#
#             Rough guide from our experience for E and e values:
#              We know that with 1000 sequence database, 0.01
#              produces around 1% error in grouping sequences
#              according to the evalue.
#              With 180,000, 0.081 gave us less than 1% error.
#             Evalue of FASTA and SSEARCH is DEPENDENT on DB size,
#              so you need to play a little bit to know the best
#              E value for your OWN database or genome.
#             The best approach is :
#               1) You run geanfammer.pl with any of your target DB
#                  with certain E value you like
#               2) Check sequence families which are clustered
#                  in the final resultant file xxxx.gclu and decide
#                  if the E value is low or high. Lower evalues will
#                  make sure you do not make wrong clusters while
#                  high evalue will include more probable sequence
#                  family members.
#               3) Put all the xxxx.msp files in subdirectory(s)
#                  created by geanfammer and run divclus.pl (which
#                  is accompanied in the package) with different
#                  Evalues. Divclus will not run any search algorithm
#                  etc, so it can be done fairly quickly.
#
#          * Most of the subroutines are found in Bio.pl or Bio.pm
#            Bioperl library for perl.
#
#
# Example   :geanfammer.pl E_gnme.fa             # simplest form
#            geanfammer.pl E_gnme.fa a=ssearch   # use SSEARCH
#            geanfammer.pl E_gnme.fa o           # for overwriting
#                                                   when you want a
#                                                   fresh run ovr old
#            geanfammer.pl E_gnme.fa c         # For keeping
#                                                 SSO files
#                                                 (fasta output)
#            geanfammer.pl E_gnme.fa k=2       # changing default
#                                                 k tuple for
#                                                 FASTA to 2
#            geanfammer.pl E_gnme.fa E=0.01     # set the E value
#                                                 for initial single
#                                                 linkage clustering
#            geanfammer.pl E_gnme.fa e=0.01    # set the E value
#                                                for domain level linkage
#       -->  geanfammer.pl E_gnme.fa e=0.01 E=0.01 # set the 2 E values
#                                                    separately (no need
#                                                    to do this)
#
# Keywords  : genome_analysis_and_protein_family_maker,
#             genome_ana_protein_fam_maker
# Options   :
#             o  for overwrite existing xxxx.fa files for search
#             c  for create SSO file (sequence search out file)
#             d  for very simple run and saving the result in
#                    xxxx.gz format in sub dir starting with one char
#             N
#             s
#             m
#             v  for debugging purpose. It says more to you while running
#             z
#             D  for making subdir like ./MG or /FA in PWD. For clean PWD
#             L  for Lean output(removes all the intermediate
#                                     outputs to save space)
#             u  for making separate summary file (redundant now)
#
#             DB=
#             File=
#             k= for k-tuple value. default is 1 (ori. FASTA prog.
#                                                   default is 2)
#             a= for choosing either fasta or ssearch algorithm
#                    You can set absolute path like (/usr/bin/fasta)
#             E= for Evalue cutoff for single linkage clustering
#                    $E_cut_main
#             e= for Evalue cutoff for divide_clusters subroutine.
#             u=
#             l=
#             d=
#             T= for minimal domain size (default is 30 aa residue)
#
#   !! Do not remove the following lines down to # Author line.
#                This program parses them!!
#
#  $factor=                 by f=     ## overlapping factor
#  $Lean_output=L           by L -L
#  $dynamic_factor=y        by y  Y -y -Y
#  $over_write=o            by o -o
#  $create_sso_file=c       by c -c
#  $k_tuple=                by k=
#  $upper_expect_limit=     by u=
#  $lower_expect_limit=     by l=
#  $algorithm=              by a=
#  $No_processing=N         by N -N
#  $single_msp=s            by s -s
#  $sequence_db_fasta=      by DB=
#  $query_file=             by File=
#  $machine_readable=M      by M -M
#  $make_subdir_out=D       by D
#  $make_subdir_gzipped=d   by d -d
#  $direct_MSP_conversion=m by m -m
#  $verbose=v               by v -v
#  $sub_dir_size=           by d=
#  $Evalue_cut_single_link= by E=
#  $Evalue_cut_divclus=     by e=
#  $optimize=z              by z -z
#  $make_separate_summary=u by u -u
#  $length_thresh=          by T=      # minimal sequence domain length threshold
#
# Author    : Sarah A Teichmann, Jong Park, sat@mrc-lmb.cam.ac.uk,
#                                      jong@salt2.med.harvard.edu
# Version   : 1.7
#------------------------------------------------------------------

print "\n\n\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "# (1) $0: The Geanfammer started ... BEFORE parsing ARGUments\n";
print "#   I will call subroutine geanfammer_main now.\n";
print "#   -Remember that default evalue used is 0.01 \n";
print "#_____________________________________________________________________\n";


@final_clu_files=@{&geanfammer_main};


#________________________________________________________________________________
# Title     : geanfammer_main
# Usage     : &geanfammer_main;
# Function  : The main sub of geanfammer
# Example   :
# Keywords  : main_geanfammer, geanfammer
# Options   :
# Version   : 2.0
#--------------------------------------------------------------------------------
sub geanfammer_main{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # All the defaults, Evalues are determined in geanfammer sub
    #__________________________________________________________________
    $|=1;
    unless($algorithm){
        $algorithm='fasta'; # default search algorithm setting(will be overridden
             # by 'a=xxx' prompt argument
    }
    $sub_dir_size=2; # this is the subdirectory name char size
    $machine_readable='M';  # this is to invoke FASTA's m=10 option
    $make_msp_in_sub_dir_opt='m';
    $make_subdir_gzipped='d';
    $make_subdir_out=$make_msp_in_sub_dir_opt='D';
    $Evalue_cut_single_link=0.01;
    $Evalue_cut_divclus    =0.01;
    $length_thresh=30;  # default

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
    # preprocessing the inputs,  parse_arguments reads in the options in the headerbox
    #__________________________________________________________________________________
    @your_genome_or_db_to_analyse_file=@{&parse_arguments(1)};

    print "\n\n\n# (2) geanfammer_main (with $algorithm): \@your_genome_or_db_to_analyse_file are(is)\n";
    print "\n => @your_genome_or_db_to_analyse_file with $algorithm. Min domain size is \"$length_thresh\"
    \$make_subdir_out is $make_subdir_out\n\n";

    if(@your_genome_or_db_to_analyse_file < 1){
        print "\n# (E) geanfammer_main: ERROR!\n";
        print "\n# Dear $ENV{'USER'}, $0: failed to find input file!\n
              Did you put FASTA format DB file as input?\n
              Or I guess your INPUT file DOES NOT exist in PWD.\n\n";
        print " As like:  $0 MG.fa \n\n\n";
        print chr(7);
        exit;
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
    # Checking if they have 'fasta' or 'ssearch' in the path
    #__________________________________________________________________________________
    if($algorithm=~/\/([^\/]+) *$/){
        print "\n# Checking \$algorithm ($algorithm) name\n";
        $algorithm_name=$1;
        if(-s $algorithm){
           print "\n# (i) $0 will use \"$algorithm\"\n";
        }else{
           $result_of_which_run=`which $algorithm_name`;
           if($result_of_which_run=~/^ *(\/\S+\/)[fastassearch]+\d* *$/){ ## after Lily Fu's suggestion
               print "\n# $0: Your $algorithm_name is in $1, good, I am running it\n";
           }else{
               print "\n# (E) \$algorithm value $algorithm is not found\n";
               print "\n# (E) $0 ran \'which\' Linux command and the result is:\n";
               print "  $result_of_which_run\n";
               print "\n# Please check your path for $algorithm\n\n"; print chr(7);
               exit;
           }
        }
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Running the actual big sub
    #______________________________
    @final_clu_files=@{&geanfammer(\@your_genome_or_db_to_analyse_file,
          "f=$factor",
          "a=$algorithm",
          "T=$length_thresh",
           $verbose,
           "d=$sub_dir_size",
           $over_write,
           $dynamic_factor,
           $create_sso_file,
           $reverse_sequence,
           $machine_readable,
           $make_msp_in_sub_dir_opt,
           $make_subdir_out,
           "E=$Evalue_cut_single_link",
           "e=$Evalue_cut_divclus",
           $make_subdir_gzipped,
           $Lean_output,
         )};
    if($verbose){
        print "\n\n\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
        print "\n#";
        print "\n# $0  is finished \n\n";
    }
    print "\n# $0 : the final output \'@final_clu_files\' is created " if $make_separate_summary;
    print "\n#__________________________________________________________________\n\n" if $verbose;
    return(\@final_clu_files);
}




#------------------------- SUB ROUTINES -----------------------------------


#______________________________________________________________________
# Title     : geanfammer
# Usage     : &geanfammer(\@your_genome_or_db_to_analyse_file,
#                          $verbose);
#
# Function  : Creates a domain level clustering file from a given
#              FASTA format sequence DB. It has been used for complete
#              genome sequence analysis. Can use PSI-blat, fasta, ssearch
#
#              ------------ USAGE INFORMATION -------------------
#             The parameters you put are important for the meaningful
#               protein family maker.
#             The most important one is the E and e options (Mostly,
#               they can have same value).
#             Large E is for setting the threshold for the single
#               linkage clustering.
#             This means, any sequence hit BELOW the threshold
#               (which is good ) will be linked.
#             For example, if Seq1 matched with Seq2 with E value
#              of FASTA search:
#              0.001, and you set the threshold 0.1, then YOU
#              ordered the geanfammer to regard them a family.
#
#             The second small e option is for the dividing a complex
#              and wrong cluster into correct more correct
#              duplication modules. This is necessary as a
#              lot of multidomain proteins can be clustered together
#              WRONGLY by single linkage.
#             At this stage, the e value is irrelevant to E value
#              and you can set a higher or lower one. Or you can set
#              the same as E.
#
#             Rough guide from our experience for E and e values:
#              We know that with 1000 sequence database, 0.01
#              produces around 1% error in grouping sequences
#              according to the evalue.
#              With 180,000, 0.081 gave us less than 1% error.
#             Evalue of FASTA and SSEARCH is DEPENDENT on DB size,
#              so you need to play a little bit to know the best
#              E value for your database or genome.
#             The best approach is :
#               1) You run geanfammer.pl with any of your target DB
#                  with certain E value you like
#               2) Check sequence families which are clustered
#                  in the final resultant file xxxx.gclu and decide
#                  if the E value is low or high. Lower evalues will
#                  make sure you do not make wrong clusters while
#                  high evalue will include more probable sequence
#                  family members.
#               3) Put all the xxxx.msp files in subdirectory(s)
#                  created by geanfammer and run divclus.pl (which
#                  is accompanied in the package) with different
#                  Evalues. Divclus will not run any search algorithm
#                  etc, so it can be done fairly quickly.
#
# Example   :geanfammer.pl E_gnme.fa             # simplest form
#            geanfammer.pl E_gnme.fa a=ssearch   # use SSEARCH
#            geanfammer.pl E_gnme.fa o           # for overwriting
#                                                   when you want a
#                                                   fresh run ovr old
#            geanfammer.pl E_gnme.fa c         # For keeping
#                                                 SSO files
#                                                 (fasta output)
#            geanfammer.pl E_gnme.fa k=2       # changing default
#                                                 k tuple for
#                                                 FASTA to 2
#            geanfammer.pl E_gnme.fa E=0.01     # set the E value
#                                                 for initial single
#                                                 linkage clustering
#            geanfammer.pl E_gnme.fa e=0.01    # set the E value
#                                                for domain level linkage
#       -->  geanfammer.pl E_gnme.fa e=0.01 E=0.01 # set the 2 E values
#                                                    separately (no need
#                                                    to do this)
#
# Keywords  : genome_analysis_and_protein_family_maker,
#             genome_ana_protein_fam_maker
# Options   :
#             o  for overwrite existing xxxx.fa files for search
#             c  for create SSO file (sequence search out file)
#             d  for very simple run and saving the result in
#                    xxxx.gz format in sub dir starting with one char
#             N
#             s
#             m
#             v
#             z
#             D
#             y  for dynamic factor
#             L  for Lean output(removes all the intermediate
#                                     outputs to save space)
#             u  for making separate summary file (redundant now)
#
#             DB=
#             File=
#             k= for k-tuple value. default is 1 (ori. FASTA prog.
#                                                   default is 2)
#             a= for choosing either fasta or ssearch algorithm
#             E= for Evalue cutoff for single linkage clustering
#                    $E_cut_main
#             e= for Evalue cutoff for divide_clusters subroutine.
#             u=
#             l=
#             d=
#             T= for $length_thresh (minimal domain size thresh)
#
#   !! Do not remove the following lines down to # Author line.
#                This program parses them
#
#  $factor=                 by f=       ## overlapping factor
#  $Lean_output=L           by L -L
#  $dynamic_factor=y        by y  Y -y -Y
#  $over_write=o            by o -o
#  $create_sso_file=c       by c -c
#  $k_tuple=                by k=
#  $upper_expect_limit=     by u=
#  $lower_expect_limit=     by l=
#  $algorithm=              by a=
#  $No_processing=N         by N -N
#  $single_msp=s            by s -s
#  $sequence_db_fasta=      by DB=
#  $query_file=             by File=
#  $machine_readable=M      by M -M
#  $make_subdir_out=D       by D
#  $make_subdir_gzipped=d   by d -d
#  $direct_MSP_conversion=m by m -m
#  $verbose=v               by v -v
#  $sub_dir_size=           by d=
#  $Evalue_cut_single_link= by E=
#  $Evalue_cut_divclus=     by e=
#  $optimize=z              by z -z
#  $make_separate_summary=u by u -u
#  $length_thresh=       by T=
#
# Author    : Sarah A Teichmann, Jong Park, sat@mrc-lmb.cam.ac.uk,
#                                      jong@salt2.med.harvard.edu
# Version   : 2.8
#--------------------------------------------------------------------------------
sub geanfammer{
    my ($i, $algorithm, $upper_expect_limit, $k_tuple, $sub_dir_size,$factor,
        @msp_files_main, $msp_directly_opt, %hash_main, $single_linkage_file,
        @written_msp_files, $num_of_seq_in_fa_file, $Evalue_cut_divclus,
        $Evalue_cut_single_link, %fasta_seqs, $final_summary_file,@file,
        @sub_clustering_clu_files, $base, $final_gclu_output, $Lean_output,
        $dynamic_factor, $length_thresh, $make_msp_in_sub_dir_opt );
    $algorithm         ='fasta'; # default algorithm set(you can change of course by using a= opt
    $upper_expect_limit= 10;     # default
    $k_tuple=1;
    $sub_dir_size=2;          # default
    $msp_directly_opt='m';
    $Evalue_cut_single_link=0.01;
    $Evalue_cut_divclus    =0.01;
    $length_thresh=30;

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # STATIC set of default E value:
    #   0.001 Evalue is very conserved for a database size roughly over 2000 seq
    #_________________________________________________________________________
    # $Evalue_cut_single_link=$Evalue_cut_divclus=0.001;

    for($i=0; $i< @_; $i++){
        if(ref($_[$i]) eq 'ARRAY'){ @file=@{$_[$i]}; next }
        elsif($_[$i]=~/a=(\S+)/){ $algorithm          = $1            }
        elsif($_[$i]=~/f=(\S+)/){ $factor = $1            }
        elsif($_[$i]=~/l=(\S+)/){ $lower_expect_limit = $1            }
        elsif($_[$i]=~/k=(\d+)/){ $k_tuple            = $1            }
        elsif($_[$i]=~/t=(\d+)/){ $Score_thresh       = $1            }
        elsif($_[$i]=~/m=(\d+)/){ $margin             = $1            }
        elsif($_[$i]=~/d=(\d+)/){ $sub_dir_size       = $1            }
        elsif($_[$i]=~/s=(\S+)/){ $single_big_msp     = $1            }
        elsif($_[$i]=~/DB=(\S+)/){ $sequence_DB       = $1            }
        elsif($_[$i]=~/File=(\S+)/){ $input_file_name = $1;           }
        elsif($_[$i]=~/Query_seqs=(\S+)/){ %seq_input = %{$1}         }
        elsif($_[$i]=~/u=(\S+)/){ $upper_expect_limit = $1            }
        elsif($_[$i]=~/E=(\S+)/){ $Evalue_cut_single_link = $1        }
        elsif($_[$i]=~/e=(\S+)/){ $Evalue_cut_divclus = $Evalue_cut_single_link = $1            }
        elsif($_[$i]=~/T= *(\S+)/){ $length_thresh = $1  }  # Length threshold Suggested by morozov@lion-ag.de
        elsif($_[$i]=~/R *$/){    $add_range          = 'r' }
        elsif($_[$i]=~/y *$/){    $dynamic_factor     = 'y' }
        elsif($_[$i]=~/z *$/){    $optimize           = 'z' }
        elsif($_[$i]=~/o *$/){    $over_write         = 'o' }
        elsif($_[$i]=~/c *$/){    $create_sso         = 'c' }
        elsif($_[$i]=~/s *$/){    $single_big_msp     = 's'; print "\n# Single file opt is set\n"; }
        elsif($_[$i]=~/m *$/){    $msp_directly_opt   = 'm' }
        elsif($_[$i]=~/M *$/){    $machine_readable   = 'M' }
        elsif($_[$i]=~/d *$/){$make_gz_in_sub_dir_opt = 'd' } # for simple search and storing in gz file (sso file will be zipped
        elsif($_[$i]=~/D *$/){$make_msp_in_sub_dir_opt= 'D' } # for simple search and storing msp file
        elsif($_[$i]=~/b *$/){    $do_in_batch        = 'b' } # for reading in all the
        elsif($_[$i]=~/n *$/){    $new_format         = 'n' }
        elsif($_[$i]=~/y *$/){    $answer_for_search  = 'y' }
        elsif($_[$i]=~/x *$/){    $answer_for_search  = 'x' }
        elsif($_[$i]=~/r *$/){    $reverse_sequence   = 'r'  }
        elsif($_[$i]=~/L *$/){    $Lean_output        = 'L'  };
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # creating xxx.msso.gz or pbla files in subdirs
    #____________________________________________
    for($i=0; $i< @file; $i++){
        $input_db_or_genome=$file[$i];
        print "\n# (i) geanfammer : input file : $input_db_or_genome, Evalue for single link:$Evalue_cut_single_link\n";
        print "# (i) geanfammer : Min. Domain size is $length_thresh\n";
        print "# (i) geanfammer : Factor is $factor\n";
        print "# (i) geanfammer : \$Evalue_cut_divclus is $Evalue_cut_divclus\n";

        $base=${&get_base_names($file[$i])};
        $final_gclu_output="$base\.gclu";

        %fasta_seqs=%{&open_fasta_files(\$input_db_or_genome)};
        $num_of_seq_in_fa_file=keys %fasta_seqs;

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # When you want to reverse your protein sequence for whatever crazy reason
        #____________________________________________________________________________
        if($reverse_sequence){ ## reverse the query seqs.
            %fasta_seqs=%{&reverse_sequences(\%fasta_seqs)};
        }

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Ask if she wants to search the db first, as it is the first time
        #_____________________________________________________________________
        unless($answer_for_search and !$over_write){
            print "\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
            print "\n# (3) You should have working FASTA (or SSEARCH) in the path";
            print "\n#     If you are doing the very first search now, Say, \"y\"es.     ";
            print "\n#_______________________________________________________________";
            print "\n# You did not set \$answer_for_search, shall I do the self search or not", "\n  (y/n)\n\n \t-->";
            $answer_for_search=getc;
        }

        if($answer_for_search=~/y/i){
            $over_write='o';

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #  By default, the search level E value is always higher by 2 !!!
            #__________________________________________________________________________
            $Evalue_cut_single_link_AT_SEARCH=$Evalue_cut_single_link+2; ## adding 2 to the given E value

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #  Searching the db self to self using default algorithm of FASTA ktup=1
            #__________________________________________________________________________
            $num_of_seq_in_fa_file=@{&do_sequence_search(
                      $Lean_output,
                      "a=$algorithm",
                      \$file[$i],
                      $over_write,
                      $msp_directly_opt,
                      "u=$upper_expect_limit",
                      $do_in_batch,
                      $create_sso,
                      $reverse_query,
                      $single_big_msp,
                      $machine_readable,
                      "DB=$input_db_or_genome",
                      "File=$input_db_or_genome",
                      $make_msp_in_sub_dir_opt,
                      $make_gz_in_sub_dir_opt,
                      "d=$sub_dir_size",
                      $new_format,
                      $add_range,
                      "E=$Evalue_cut_single_link_AT_SEARCH",
                    ) };
        }else{
            print "\n# (i) $0; You chose NO search. I will cluster and divide now\n";
            print "# (i) The min domain size var \$length_thresh is $length_thresh\n";
            $over_write='o';
            sleep(1);
        }

        @msp_files_main=@{&get_all_msp_files};

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # When users didn't define E value and no STATIC Eval is set,
        #    estimate according to the db size
        #______________________________________________________________________
        if( (!$Evalue_cut_single_link and !$Evalue_cut_divclus) or $dynamic_factor){
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # following is a very rough guide for a reasonable E value thresh for different DB size
            #______________________________________________________________________________________________
            if(@msp_files_main > 280,000){     $Evalue_cut_single_link=$Evalue_cut_divclus=0.2;
            }elsif(@msp_files_main > 180000){  $Evalue_cut_single_link=$Evalue_cut_divclus=0.1;
            }elsif(@msp_files_main > 50000){   $Evalue_cut_single_link=$Evalue_cut_divclus=0.08;
            }elsif(@msp_files_main > 10000){   $Evalue_cut_single_link=$Evalue_cut_divclus=0.07;
            }elsif(@msp_files_main > 1000){    $Evalue_cut_single_link=$Evalue_cut_divclus=0.05;
            }elsif(@msp_files_main > 500){     $Evalue_cut_single_link=$Evalue_cut_divclus=0.03;
            }elsif(@msp_files_main > 100){     $Evalue_cut_single_link=$Evalue_cut_divclus=0.02;
            }elsif(@msp_files_main > 50){      $Evalue_cut_single_link=$Evalue_cut_divclus=0.01;
            }elsif(@msp_files_main > 20 ){     $Evalue_cut_single_link=$Evalue_cut_divclus=0.001;
            }
        }elsif(!$Evalue_cut_single_link ){ $Evalue_cut_single_link=$Evalue_cut_divclus
        }elsif(!$Evalue_cut_divclus ){     $Evalue_cut_divclus    =$Evalue_cut_single_link }


        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Making single linkage cluster file xxxxxxxx.sclu
        #___________________________________________________________________________________________
        print "\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`#\n";
        print "# BEfore running Sinlge linkage sub                    #\n";
        print "# Evalue single for clus & Divclus  : $Evalue_cut_single_link, $Evalue_cut_divclus     #\n";
        print "#______________________________________________________#\n";
        %hash_main=%{&msp_single_link_hash(\@msp_files_main, $Evalue_cut_single_link)};


        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # From Single LINKAGE file, I compile each wrong cluster MSP files resulting in bigger MSPs
        #___________________________________________________________________________________________
        $single_linkage_file=${&print_clusfile_from_hash(\%hash_main)};
        if(-s $single_linkage_file){
            print "\n# (i) Good! \$single_linkage_file $single_linkage_file is made\n";
            @written_msp_files=@{&convert_clu_to_msp(\$single_linkage_file, $Lean_output)};
        }else{
            print "\n# (E) BAD! \$single_linkage_file $single_linkage_file is NOT here\n";
            exit;
        }
        print "\n\n# (i) \@written_msp_files by convert_clu_to_msp:\n @written_msp_files\n\n\n";

        $over_write='o';
        $average_region='A';
        $range='r';
        $merge='m';
        $sat_file=0;
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # >>>>>> !! The very core of geanfammer. DIVCLUS
        #___________________________________________________________________________
        @sub_clustering_clu_files=@{&divide_clusters(
                $Lean_output,
                \@written_msp_files,
                "s=$score",
                "f=$factor", ## this is a very impo. parameter in the behaviour of divclus, Sarah!
                "T=$length_thresh",
                "e=$Evalue_cut_divclus", ## this is a very impo. parameter in the behaviour of divclus, Sarah!
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
               )};

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  This is the final result now.
        #_______________________________________________________________

        print "\n\n# create_sorted_cluster: Final gclu file name is $final_gclu_output\n  with @sub_clustering_clu_files\n" if $verbose;
        $good_cluster_file=${&create_sorted_cluster(\@sub_clustering_clu_files, \$final_gclu_output)};
        print "\n\n\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
        print "\n# ! Finished create_sorted_cluster(geanfammer), result is \"$good_cluster_file\"";
        print "\n#_________________________________________________________________\n\n\n";

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  make_clustering_summary  became redundant after peer review
        #_______________________________________________________________
        $final_summary_file=${&make_clustering_summary(\$final_gclu_output, $make_separate_summary)};

        push(@files_created, $final_summary_file);

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Make the whole output Lean
        #_______________________________________
        if($Lean_output){
            for($k=0; $k < @sub_clustering_clu_files; $k++){
                $base_name=${&get_base_names($sub_clustering_clu_files[$k])};
                $sub_clustering_msp_file="$base_name\.msp";
                unlink($sub_clustering_clu_files[$k], $sub_clustering_msp_file);
            }
        }
   }# end of for loop for @file
   return(\@files_created);
}



#_______________________________________________________________________
# Title     : divide_clusters
# Usage     : &divide_clusters(\@file);
# Function  : This is the main funciton for divclus.pl
#               divides complex single linkage cluster into smaller duplication
#               module level sub clusters.
# Example   : &divide_clusters(\@file, $verbose, $range, $merge, $sat_file,
# 	                $dindom, $indup, "T=$length_thresh", "e=$evalue", $over_write,
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
#   o  for $over_write
#
# Version   : 2.9
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
    my($merge, $verbose, $sat_file, $length_thresh, $factor, $indup, $indup_percent,
       $score, @temp_show_sub, $optimize, $file, $evalue, $over_write, $din_dom,
       $sum_seq_num, $base_1, $output_clu_file, $short_region, $large_region,
       $average_region, $dynamic_factor, @sub_clustering_clu_files);

    $factor=7; # default factor is 7 for 70%
    $length_thresh=30;

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dealing with options
    #_________________________________________
    if($char_opt=~/m/){        $merge='m';
    }if($char_opt=~/v/){       $verbose='v'; # for showing debugging information
    }if($char_opt=~/i/){       $indup='i';
    }if($char_opt=~/z/){       $optimize='z';
    }if($char_opt=~/o/){       $over_write='o';
    }if($char_opt=~/d/){       $din_dom='d';
    }if($char_opt=~/s/){       $sat_file='s';
    }if($char_opt=~/y/){       $dynamic_factor='y';
    }if($char_opt=~/S/){       $short_region  ='S';
    }if($char_opt=~/L/){       $large_region  ='L';
    }if($char_opt=~/A/){       $average_region='A';
    }if($vars{'T'}=~/\d+/){    $length_thresh= $vars{'T'};
    }if($vars{'l'}=~/\d+/){    $length_thresh= $vars{'l'}; ## synonym of 't'
    }if($vars{'f'}=~/\S+/){    $factor= $vars{'f'};
    }if($vars{'s'}=~/\d+/){    $score = $vars{'s'};
    }if($vars{'e'}=~/\d+/){    $evalue= $vars{'e'};
    }if($vars{'E'}=~/\d+/){    $evalue= $vars{'E'}; # synonym of e
    }

   $percent_fac=$factor*10;

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # (0) When one file input was given (yes, divclus can handle multiple files, Sarah!)
   #________________________________________________________________________________
   if(@file == 1){  #<=== @file has xxxx.msp, yyyy.msp  zzzz.msp ....,
        print "\n# (1) divide_clusters: One single file was given=> \"@file\"\n";
        $file=$file[0];
        $base_1=${&get_base_names($file)};

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # (2) Define the output cluster file name:  eg, 3-232_cluster_F7.clu , F7 means factor used is 7
        #______________________________________________________________________________________________
        $output_clu_file="$base_1\_F${factor}\.clu";

        if( !$over_write and -s $output_clu_file){
            print "\n# $output_clu_file Already EXISTS, skipping. Use \'o\' opt to overwrite\n"; exit;
        }
        print "# (2) divide_clusters: processing ONE single file \"@file\" with merge_sequence_in_msp_file\n" if $verbose;

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # (3) merge_sequence_in_msp_file does not do much. Just filtering and producing
        #     sequences in ISPA_PBS_21-215 VPR_PBS_160-354 format from msp format
        #________________________________________________________________________________
        @grouped_seq_names=@{&merge_sequence_in_msp_file(\@file, "s=$score", $optimize, $din_dom, $sat_file,
                $optimize, "T=$length_thresh", "e=$evalue", "f=$factor", "$range", "$merge", $verbose,
                $short_region, $large_region, $average_region, $over_write, $dynamic_factor)};

        if($verbose){
            print "\n\n# (3) divide_clusters: finished running \"merge_sequence_in_msp_file\"\n  ==>";
            for($i=0; $i< @grouped_seq_names; $i++){
               print "\n-->> $grouped_seq_names[$i]";
            }
        }

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # (4) This is critical seqlet merging step
        #________________________________________________________________________________
        @out=@{&cluster_merged_seqlet_sets(\@grouped_seq_names, $dynamic_factor,
               "f=$factor", $short_region, $large_region, $average_region, $optimize)};

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # (5) This is showing the result in clu file format
        #________________________________________________________________________________
        @temp_show_sub=&show_subclusterings(\@out, $file, $sat_file, $dindom, $indup,
						   "e=$evalue", "p=$percent_fac", "f=$factor" );
        $good_bad       = $temp_show_sub[0];
        $indup_c        = $temp_show_sub[1];
        $sum_seq_num   += $temp_show_sub[2];
        push(@sub_clustering_out_files, @{$temp_show_sub[3]});

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
            my (@grouped_seq_names, @temp_show_sub);
            my $indup_c=0;
            my $big_msp_file=$file[$i];
            unless(-s $big_msp_file){ print "\n# (E) \$big_msp_file does not exist\n"; exit }

            $base_1=${&get_base_names($big_msp_file)};
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # (1) Define the output cluster file name:  eg, 3-232_cluster_F7.clu , F7 means factor used is 7
            #______________________________________________________________________________________________
            $output_clu_file="$base_1\_F${factor}\.clu";

            print "\n# DIVCLUS: just before merge_sequence_in_msp_file, \$output_clu_file is $output_clu_file from input file $big_msp_file" if $verbose;
            if( !$over_write and -s $output_clu_file){
                print "\n# $output_clu_file Already EXISTS, skipping. Use \'w\' opt to overwrite\n";
                next;  }

            print "\n# (1)  divide_clusters: processing file \"$big_msp_file\" for $output_clu_file" if $verbose;

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #  (2) If clu file(eg 2-1618_ss.clu ) is in pwd, tries to skip
            #____________________________________________________________
            if((-s $output_clu_file) > 10 and $over_write !~/o/){
                print "# $output_clu_file exists, skipping, use \"o\" option to overwrite\n";  next;
            }

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # (3) merge_sequence_in_msp_file does not do much. Just filtering and producing
            #     sequences in ISPA_PBS_21-215 VPR_PBS_160-354 format of STRING from msp format
            #     $big_msp_file is an MSP file
            #________________________________________________________________________________
            print "\n# (i) divide_clusters : I am merging seq in $big_msp_file file for $output_clu_file\n" if $verbose;
            @grouped_seq_names=@{&merge_sequence_in_msp_file(\$big_msp_file, "s=$score", $din_dom, $sat_file, $optimize,
                                "T=$length_thresh", "e=$evalue", "f=$factor", "$range", "$merge", $verbose, $over_write,
                                 $short_region, $large_region, $average_region, $dynamic_factor )};

            if($verbose){
                print "\n# \@file has more than one input file\n # The result of \"merge_sequence_in_msp_file\"\n";
                print "@grouped_seq_names";
            }

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #  (4) Clustering the sets of merged seqlets => CORE algorithm
            #____________________________________________________________

            @out=@{&cluster_merged_seqlet_sets(\@grouped_seq_names, "f=$factor", $optimize, $dynamic_factor,
                   $short_region, $large_region, $average_region)};

            @temp_show_sub=&show_subclusterings(\@out, $big_msp_file, $sat_file, $dindom, $indup,
                                                    "e=$evalue", "p=$percent_fac", "f=$factor");
                        $good_bad       = $temp_show_sub[0];
                        $indup_c        = $temp_show_sub[1];
                        $sum_seq_num   += $temp_show_sub[2];
            push(@sub_clustering_out_files, @{$temp_show_sub[3]});

            if($good_bad==1){          push(@good, $big_msp_file);
            }else{         push(@bad, $big_msp_file);       }

          }
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
          &write_good_bad_list_in_divide_clusters(\@good, \@bad);
          sub write_good_bad_list_in_divide_clusters{
               my  (@good, @bad); @good=@{$_[0]}; @bad=@{$_[1]};
               open(GOODBAD, ">good_bad.list");
               print GOODBAD "GOOD: all link    : 000\n";
               for($i=0; $i< @good; $i++){  print GOODBAD "$good[$i]\n";  }
               print GOODBAD "BAD : Not all link: 000\n";
               for($i=0; $i< @bad; $i++){   print GOODBAD "$bad[$i]\n";   }
               close(GOODBAD);
          }
          #_______________________________________________________________

   }
   return(\@sub_clustering_out_files); # contains (xxxx.clu, yyy.clu,, )
}



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
        mkdir ($dir,  0777) unless -d $dir;  ## Jong changed
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

#___________________________________________________________________________________________
# Title     : sw_search
# Usage     : &sw_search($query_database, $target_database, $ssearch_version_to_use
# Function  : to search one database against the other using ssearch (default is simply "ssearch")
# Example   :
# Keywords  : ssearch, database_ssearch
# Options   :
# Author    : Sarah A. Teichmann
# Date      : 19th September 1997
# Version   : 1.0
#--------------------------------------------------------------------------------------------
sub sw_search{
my ($qdb, $tdb, @qdbcont, $fastaver, $gene, $seq, @genes, %genes, $genes,
    $out, $tmp, $sw_score, $e_val, @tmpcontent, $i, $dir, @dir);
    $qdb=$_[0];
    open (QDB, "$qdb");
    @qdbcont=<QDB>;
    close QDB;
    $tdb=$_[1];
    open(TDB, "$tdb");
    close TDB;
    if ($_[2]){$fastaver=$_[2];}
    else{$fastaver="ssearch";}

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
        @dir=splice(@dir,0,2);
        $dir=join('',@dir);
        mkdir ("$dir", 0777) unless -d $dir;
        $out="$dir"."/"."$genes".".sso";
        if (-e $out){next;}
        $tmp=&tempname;
        open (TMP, ">$tmp");
        print TMP ">$genes\n", $genes{"$genes"}, "\n";
        close TMP;
        $sw_score=0;
        $e_val=10;
        @tmpcontent=`$fastaver -E 0.1 -H -m 10 $tmp $tdb`;
        open (OUT, ">$out");
        print OUT "@tmpcontent\n";
        close OUT;
        unlink ("$tmp");
        next;
    }

}


#________________________________________________________________________________
# Title     : get_all_msp_files
# Usage     : @msp_files=@{&get_all_msp_files};
# Function  : puts the names of all msp or msp.gz files in the directory and its subdirectories into an array
# Example   :
# Keywords  : read_msp_files, make_msp_file_array
# Date      : 19th September 1997
# Options   :
# Version   : 1.2
# Author    : Sarah A. Teichmann, jong@salt2.med.harvard.edu
#--------------------------------------------------------------------------------
sub get_all_msp_files {
   my (@msp_files_main_final);
   my @msp_files_main=@{&read_file_names_only('.msp','.msp.gz')};
   for(@msp_files_main){
      if($_=~/cluster\.msp/){ next
      }else{ push(@msp_files_main_final, $_) }
   }
   my @dirs=@{&read_dir_names_only('n', '.')};
   my ($i, $j, @msp_files);
   for ($i=0; $i<@dirs; $i++){
       my $dir=$dirs[$i];
       unless( -d $dir){
          next;
       }
       if( -d $dir){
          chdir($dir);
       }

       my @msp_files=@{&read_file_names_only('.','msp','msp.gz')};
       for ($j=0; $j<@msp_files; $j++){
          if($msp_files[$j]=~/cluster\.msp/){ next } ## skipping already existing  20-1_cluster.msp like files
          my $file_in_dir=$msp_files[$j];
          my $dir_file="$dir"."/"."$file_in_dir";
          push(@msp_files_dirs, $dir_file);
          next;
       }
       chdir('..');
       next;
   }

   @msp_files=(@msp_files_main_final, @msp_files_dirs);
   sort @msp_files;
   return (\@msp_files);
}

#______________________________________________________________
# Title     : get_largest_file
# Usage     : $largest_file_name=${&get_largest_file(@ARGV)};
# Function  : checks the size of files and returns the largest
#             one's name. If a file is not present in pwd or
#             specified absolute path, it ignores it.
# Example   :
# Keywords  : choose_largest_file, largest_file, find_largest_file
#             get_the_largest_file, choose_the_largest_file, get_biggest_file
#             fetch_largest_file, take_largest_file, get_bigger_file, get_larger_file
# Options   : _  for debugging.
#             #  for debugging.
#             e  for extract the largest from the input array
#                       leaving it one element less, in this case
#                       there will be two returning refs.
# Category  :
# Version   : 1.4
#--------------------------------------------------------------
sub get_largest_file{
    my @in;
    if(ref $_[0] eq 'ARRAY'){
	 @in = @{$_[0]};
    }else{
	 @in = @_;
    }

    my ($largest_file, $largest, $i, $extract_opt);
    for($i=0; $i< @in; $i++){
	if(($in[$i]=~/^\-?e$/i)&&(!(-f $in[$i])) ){
	   $extract_opt=1;
	   splice(@in, $i, 1);
	   $i--;
	}
    }
    for($i=0; $i< @in; $i++){
	my $size=(-s $in[$i]);
	if($size > $largest){
             $largest=$size;
             if($extract_opt ==1){
                    print "\$extract_opt is $extract_opt \n";
                    push(@in, $largest_file) if defined($largest_file);
                    $largest_file = splice(@in, $i, 1);
                    print "\n $largest_file \n";
                    $i--;
             }else{
                    $largest_file=$in[$i];
             }
        }
    }
    if($extract_opt==1){
	 return(\$largest_file, \@in);
    }else{ return(\$largest_file); }
}


#______________________________________________________________
# Title     : get_smallest_file
# Usage     : $smallest_file_name=${&get_largest_file(@ARGV)};
# Function  : checks the size of files and returns the smallest
#             one's name. If a file is not present in pwd or
#             specified absolute path, it ignores it.
# Example   :
# Keywords  : choose_smallest_file, smallest_file, find_smallest_file
#             get_the_smallest_file, choose_the_smallest_file,
#             fetch_smallest_file, take_smallest_file, get_smaller_file,
# Options   : _  for debugging.
#             #  for debugging.
#             e  for extract the smallest from the input array
#                       leaving it one element less, in this case
#                       there will be two returning refs.
# Category  :
# Version   : 1.3
#--------------------------------------------------------------
sub get_smallest_file{
    my @in;
    if(ref $_[0] eq 'ARRAY'){
	 @in = @{$_[0]};
    }else{
	 @in = @_;
    }
    my $smallest=10000000000;
    my ($smallest_file, $i, $extract_opt);

    for($i=0; $i< @in; $i++){
	if(($in[$i]=~/^\-?e$/i)&&(!(-f $in[$i])) ){
	   $extract_opt=1;
	   splice(@in, $i, 1);
	   $i--;
	}
    }
    for($i=0; $i< @in; $i++){
	my $size=(-s $in[$i]);
	if($size < $smallest){
	   $smallest=$size;
	   if($extract_opt ==1){
		  print "\$extract_opt is $extract_opt \n";
		  push(@in, $smallest_file) if defined($smallest_file);
		  $smallest_file = splice(@in, $i, 1);
		  print "\n $smallest_file \n";
		  $i--;
	   }else{
		  $smallest_file=$in[$i];
	   }
	}
    }
    if($extract_opt==1){
         return(\$smallest_file, \@in);
    }else{ return(\$smallest_file); }
}
#______________________________________________________________________________
# Title     : find_program_in_path
# Usage     :
# Function  :
# Example   :
# Keywords  : which, whence
# Options   :
# Author    : Jong, jong@salt2.med.harvard.edu
# Category  :
# Version   : 1.3
#------------------------------------------------------------------------------
sub find_program_in_path{
    my(@found_prog_and_path, $file, $verbose);
    print "\n\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`";
    print "\n# (INFO) find_program_in_path: Running inside of $0  ";
    print "\n#________________________________________________\n\n";

    if(@_==0){
       print "\n# (ERROR) $0 requires at least one program in path to detect\n\n";
       print chr(7);
       exit;
    }
    foreach $file (@_) {
       if($file=~/^v/i and !-s $file){
          $verbose='v'; next;
       }
       for $dir (split(/:/,$ENV{'PATH'})) {
          if(-f ($path="$dir/$file") && -f _){
             $size= -s "$path";
             $path=~s/\/\//\//g;  ## removing /bin//gzip like entry
             push(@found_prog_and_path, sprintf("%-30s %12s byte", $path, $size) );
          }
       }
    }
    @found_prog_and_path=@{&remove_dup_in_array(\@found_prog_and_path)};

    if(@found_prog_and_path < 1){
       print "\n# (INFO) $0 could not find any match for @ARGV, Sorry. \n\n";
    }else{
       print "\n# (INFO) find_program_in_path: There were ", scalar(@found_prog_and_path)," matches found for $file\n";
    }
    if($verbose){
       &show_array(\@found_prog_and_path);
    }
    return(\@found_prog_and_path);
}



#____________________________________________________________________________________
# Title     : do_sequence_search
# Usage     : &do_sequence_search("Query_seqs=\%pdb_seq", "DB=$sequence_db_fasta",
#  		         "File=$ARGV[0]", $single_msp, $over_write,
# 	        	 "u=$upper_expect_limit", "l=$lower_expect_limit",
#       		 "k=$k_tuple", $No_processing );
# Function  : do FASTA, SSEARCH or BLASTPGP(psi-blast) search
# Example   : &do_sequence_search(\%pdb_seq, $owl_db_fasta, $ARGV[0], $single_msp, $over_write,
#                    "u=$upper_expect_limit", "l=$lower_expect_limit", "k=$k_tuple" );
#
# Keywords  : sequence_search
# Options   :
#             Query_seqs=  for enquiry sequences eg)  "Query_seqs=$ref_of_hash"
#             DB=   for target DB  "DB=$DB_used"
#             File= to get file base(root) name.  "File=$file[0]"
#             m  for MSP format directly from FASTA or Ssearch result than through sso_to_msp to save mem
#             s  for the big single output (msp file output I mean)
#             s= for the single big msp file name
#             o  for overwrite existing xxxx.fa files for search
#             c  for create SSO file (sequence search out file)
#             d  for very simple run and saving the result in xxxx.gz format in sub dir starting with one char
#             r  for reverse the query sequence
#             R  for attaching ranges of sequences
#             k= for k-tuple value. default is 1 (ori. FASTA prog. default is 2)
#             u= for $upper_expect_limit
#             l= for $lower_expect_limit
#             a= for choosing either fasta or ssearch algorithm
#             d= for defining the size of subdir made. 2 means it creates
#                    eg, DE while 1 makes D
#             d  for $make_gz_in_sub_dir_opt, putting resultant sso files in gz format and in single char subdir
#             D  for $make_msp_in_sub_dir_opt, convert sso to msp and put in sub dir like /D/, /S/
#             n  for new format to create new msp file format with sso_to_msp routine
#          PVM=  for PVM run of FASTA (FASTA only)
#             M  for machine readable format -m 10 option
#             M= for machine readable format -m 10 option
#             N  for 'NO' do not do any processing but, do the searches only.
#       FILE_AGE for defining the age of file in days to be overwritten.
#             L  for Lean output(removes xxxx.fa query seq file)
#
# Returns   : the names of files created (xxxxx.msp, yyy.msp,,)
# Version   : 5.4
#----------------------------------------------------------------------------------------
sub do_sequence_search{
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
    my (@final_out, $add_range, $single_big_msp, $base_name, $create_sso, @nondup,
       $Single_msp_out_file, %duplicate, $Evalue_thresh, $Score_thresh, @SSO, $sequence_DB,
       @sso, @temp, $algorithm, $margin, $out_msp_file, @MSP, @final_msp_file_names_out,
       $upper_expect_limit, $lower_expect_limit, $k_tuple, %seq_input, %MSP, $No_processing,
       $new_format, $PVM_FASTA_run, $over_write, $sub_dir_size, $age_in_days_of_out_file,
       $over_write_by_age, $Lean_output, $gzipped_msp_file, $gzipped_sso_file,
       $defined_all_ok, $make_msp_in_sub_dir_opt );
    my ($E_val) = 5;  ## default 5 <<<<<<<<<<<<<<<<<<<<<

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # DEFAULTS
    #________________________________________
    $k_tuple           =1;  # 1 or 2, 1 is more sensitive
    $algorithm         ='fasta';
    $sub_dir_size      =2;  # the default char number taken from seq name to make sub dirs
    $upper_expect_limit=1;
    $lower_expect_limit=0;
    $Score_thresh      =75; # FASTA or SSSEARCH score
    $margin            =0;  # sequence region margin. If it is 2, 2 more edged residues will be added
    $add_range         ='';
    $pwd               =`pwd`; chomp($pwd);
    $age_in_days_of_out_file=1000; ## If the files being tested is older than this, let's make anew

    if($vars{'a'}=~/\S+/){ $algorithm          = $vars{'a'}            };
    if($vars{'u'}=~/\d+/){ $upper_expect_limit = $vars{'u'}            };
    if($vars{'l'}=~/\d+/){ $lower_expect_limit = $vars{'l'}            };
    if($vars{'k'}=~/\d+/){ $k_tuple            = $vars{'k'}            };
    if($vars{'t'}=~/\d+/){ $Score_thresh       = $vars{'t'}            };
    if($vars{'m'}=~/\d+/){ $margin             = $vars{'m'}            };
    if($vars{'d'}=~/\d+/){ $sub_dir_size       = $vars{'d'}            };
    if($vars{'r'}=~/\S+/){ $add_range          = 'r'                   };
    if($vars{'s'}=~/\S+/){ $single_big_msp     = 's'                   };
    if($vars{'DB'}=~/\S+/){            $sequence_DB=$vars{'DB'} ;
        if(-s $sequence_DB){
        }elsif(-s "../$sequence_DB"){  $sequence_DB= "../$sequence_DB"
        }elsif(-s "../../$sequence_DB"){  $sequence_DB= "../../$sequence_DB";
        }else{
            print "\n# (E) do_sequence_search: You set DB param, but I can\'t find $sequence_DB\n";
            exit;
        }
    }else{  print "\n# (E) do_sequence_search: I need DB param defined, sorry, aborting\n"; }

    if($vars{'FILE'}=~/\S+/){ $input_file_name = $vars{'FILE'}; push(@file,$input_file_name) };
    if($vars{'File'}=~/\S+/){ $input_file_name = $vars{'File'}; push(@file,$input_file_name) };
    if($vars{'FILE_AGE'}=~/\S+/){ $age_in_days_of_out_file= $vars{'FILE_AGE'};  };
    if($vars{'Query_seqs'}=~/\S+/){ %seq_input = %{$vars{'Query_seqs'}}};
    if($vars{'Query'}=~/\S+/){      %seq_input = %{$vars{'Query'}}};
    if($vars{'u'}    =~/\S+/){ $E_val          = $vars{'u'}            };
    if($vars{'PVM'}  =~/\S+/){ $PVM_FASTA_run  = $vars{'PVM'}; print "\n# PVM opt is set\n";     };
    if($vars{'M'}  =~/\S+/){ $machine_readable = $vars{'M'};           };

    if($char_opt=~/r/){    $add_range          = 'r' }
    if($char_opt=~/L/){    $Lean_output        = 'L' }
    if($char_opt=~/o/){    $over_write         = 'o' }
    if($char_opt=~/c/){    $create_sso         = 'c' }
    if($char_opt=~/s/){    $single_big_msp     = 's'; print "\n# Single file opt is set\n"; }
    if($char_opt=~/m/){    $msp_directly_opt   = 'm' }
    if($char_opt=~/M/){    $machine_readable   = 'M' }
    if($char_opt=~/d/){    $save_in_gz_in_sub_dir  = 'd' }
    if($char_opt=~/D/){$make_msp_in_sub_dir_opt= 'D' } # for simple search and storing msp file
    if($char_opt=~/N/){    $No_processing      = 'N'; $create_sso='c'; }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
    # When no %seq is given, but files
    #___________________________________________
    if(@hash==0 and @file > 0){
       print "\n# (i) do_sequence_search: You did not put sequences as in \%seq, but raw sequence file @file!\n";
       print "        I will run \'open_fasta_files\' sub to fetch sequences to store in \%seq_input\n";
       %seq_input=%{&open_fasta_files(\@file)};
    }else{
       print "\n# (i) do_sequence_search: I will use given seqs in \%seq_input from \%\{\$hash\[0\]\}\n";
       %seq_input=%{$hash[0]};
    }
    my (@seq_names) = keys %seq_input;
    $base_name = ${&get_base_names($input_file_name)};
    print "\n# (i) line:",__LINE__, ", do_sequence_search, \$algorithm => $algorithm, \$base_name:$base_name
               $input_file_name <--> $sequence_DB\n";

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```
   # If one of the files is missing, ask if it is for self self search and
   #  make query=target, else die
   #_______________________________________________________________________
   $defined_all_ok=&check_if_defined($input_file_name, $sequence_DB);
   unless($defined_all_ok){
       print "\n  Did you want to do self self search? ->(y/n) ";
       $answer_for_self_self=getc;
       if($answer_for_self_self =~/y/i){
           if($input_file_name){  $sequence_DB=$input_file_name }
           else{ $input_file_name=$sequence_DB }
       }else{
           print "\n# You seemed made a mistake, O.K., I will kill myself!\n\n";
           print chr(7);  exit;
       }
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # (0) If blast is chosen run Blast
   #_________________________________________________________
   if($algorithm=~/[psi\-]*[pb][last]*/i){
      print "\n# (i) Doing PSI search with @file\n";
      @final_out=@{&do_psi_blast_search(\@file, "d=$source_DB_file",
                           "i=$input_seq_file",  $over_write,
                           $make_msp_in_sub_dir_opt, $Lean_output)};
      return(\@final_out); #<<<<<<----------- F I N I S H
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # (1) Controlling which kind of search it should do. Do save_in_gz_in_sub_dir first if d is set
   #______________________________________________________________________________________________
   if( $make_msp_in_sub_dir_opt ){  ## convert sso to msp and put in sub dir like /D/, /S/

       for($x=0; $x < @seq_names; $x++){
          my ($over_write_sso_by_age, $over_write_msp_by_age,  %single_seq,
              $out_file_sso_gz_name, $out_file_msp_name, $out_file_gz_name, $existing_sso);
          my ($seq_name, $seq)= ($seq_names[$x], $seq_input{$seq_names[$x]});
          my $first_char= substr("\U$seq_name", 0, $sub_dir_size);
          mkdir ("$first_char", 0777) unless -d $first_char;
          chdir("$first_char");
          print "\n# (i) do_sequence_search: You set \'d\' or \'D\' opt\n";
          print "# (i) making subDIRs ($first_char) with $seq_name $sequence_DB to store MSP files\n";

          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # Let's make each fasta file for each seq to be used in searching
          #_____________________________________________________________________
          my $temp_file_name="$seq_name.fa";
          %single_seq=($seq_name, $seq_input{$seq_name});
          &write_fasta(\%single_seq, $temp_file_name ); ## e for writing each file

          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # Making output file name according to the option given
          #_______________________________________________________________________
          if($machine_readable and $algorithm=~/[fastassearch]+/){
                 $out_file_sso_name="$seq_name\.msso";
          }else{ $out_file_sso_name="$seq_name\.sso";      }
          $out_file_sso_gz_name="$out_file_sso_name\.gz";
          $out_file_msp_name="$seq_name\.msp";
          $out_file_gz_name="$seq_name\.msp\.gz";

          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # Check if SSO file already there
          #_______________________________________________________________________
          if(-s $out_sso_file){ $existing_sso=$out_file_sso_name }
          elsif(-s $out_sso_gz_name){ $existing_sso=$out_file_sso_gz_name }
          if(-s $out_msp_name){ $existing_msp=$out_file_msp_name }
          elsif(-s $out_gz_name){ $existing_msp=$out_file_gz_name }

          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # If the dates of files created are long ago, overwrite to refresh
          #____________________________________________________________________
          if(  (localtime(time- (stat($existing_sso))[9]))[3] > $age_in_days_of_out_file ){
               $over_write_sso_by_age='o';
          }
          if(  (localtime(time- (stat($existing_msp))[9]))[3] > $age_in_days_of_out_file ){
               $over_write_msp_by_age='o';
          }

          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          #  To check if the target seq DB is in ../
          #________________________________________________
          if(-s $sequence_DB){
              print "\n# (i) Good, target \$sequence_DB $sequence_DB is in this working dir\n";
          }elsif( -s "../$sequence_DB"){ $sequence_DB="../$sequence_DB"; }

          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # (2) Searching: Making MSP files directly,  MSP file format is the major format used in geanfammer!
          #_________________________________________________________________________________________
          if($char_opt =~/D/){ #### To make MSP file
               if( !(-s $out_file_gz_name or -s $out_file_msp_name) or $over_write or $over_write_msp_by_age){
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    # (2.1) Running  run_fasta_sequence_search !!
                    #_______________________________________________________
                    print "\n# (i) Running  run_fasta_sequence_search !!\n";
                    $gzipped_msp_file=${&run_fasta_sequence_search(
                                       "a=$algorithm",
                                       "O=$out_file_msp_name",
                                       "File=$temp_file_name", "e=$E_val",
                                       "DB=$sequence_DB", "k=$k_tuple", "$machine_readable")};
                    $gzipped_sso_file=${&compress_files_by_gzip($out_file_sso_name)};
               }else{
                   print "\n#  Line No. ", __LINE__,", $out_file_gz_name already exists or
                         \$over_write is set or NOT older than $age_in_days_of_out_file\n";
               }
          }
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # To make gzipped SSO files and MSP files
          #_______________________________________________
          elsif($create_sso or $char_opt=~/m/){ ### To make gzipped SSO files
               if( !(-s $out_file_sso_name or -s $out_file_sso_gz_name ) or $over_write or $over_write_sso_by_age){
                   print "\n# (i) Running  run_fasta_sequence_search with \"\$create_sso option\"!!\n\n";
                   $gzipped_msp_file=${&run_fasta_sequence_search(
                                       "a=$algorithm",
                                       "O=$out_file_msp_name", "$create_sso",
                                       "File=$temp_file_name", "e=$E_val",
                                       "DB=$sequence_DB", "k=$k_tuple", "$machine_readable")};

                   $gzipped_sso_file=${&compress_files_by_gzip($out_file_sso_name)};
               }else{
                   print "\n#  Line No. ", __LINE__,", $out_file_gz_name already exists or
                         \$over_write is set or NOT older than $age_in_days_of_out_file\n";
               }
          }else{
               if( !(-s $out_file_sso_name or -s $out_file_sso_gz_name ) or $over_write or $over_write_sso_by_age){
                   system(" $algorithm -m 10 -H  -E $E_val $temp_file_name $sequence_DB $k_tuple > $out_file_sso_name");
                   system("gzip $out_file_sso_name");
               }else{
                   print "\n#  Line No. ", __LINE__,", $out_file_gz_name already exists or
                         \$over_write is set or NOT older than $age_in_days_of_out_file\n";
               }
          }
          unlink("$seq_name.fa") if -s "$seq_name.fa";
          unlink("$first_char/$seq_name\.fa") if -s "$first_char/$seq_name\.fa";
          print "\n# Sub dir $first_char and $seq_name\.msp has been made, finishing do_sequence_search\n";
          chdir ('..');
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # F I N I S H
      #________________________________________
      goto EXIT;
   } # if ($char_opt =~/[dD]/){  is finished


   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # (2) Writing on PWD. This is the big single MSP output
   #____________________________________________________________
   $Single_msp_out_file="$base_name\.msp" if($single_big_msp eq 's');
   if(-s $Single_msp_out_file and !$over_write ){
       print "\n# (i) $Single_msp_out_file exists, and no \$over_write is set, skipping \n";
       push(@final_out, $Single_msp_out_file);
   }else{  $over_write  ='o';  }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Check if it is necessary to write each sequences.fa files
   #______________________________________________________
   if( $over_write ){  &write_fasta_seq_by_seq(\%seq_input, 'e'); } ## e for writing each seq file

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   #  When, you didn't use "DB=$XXX" and "File=$FXXX" format, first file input is DB etc
   #_______________________________________________________________________________________
   $defined_all_ok=&check_if_defined($input_file_name, $sequence_DB);
   if(!$defined_all_ok){ print "\n# (E) FATAL: do_sequence_search: You did not use \"DB=\$XXX\" format\n"; exit   };

   print "\n# Finished writing the enquiry fasta files from \%seq_input by write_fasta";
   print "\n# I am in do_sequence_search sub, Target database used :  $sequence_DB with seqs of \'@seq_names\'\n";


   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  Main search with given @seq_names
   #______________________________________________________________
   for($j=0; $j< @seq_names; $j++){  # @seq_names has sequence names coming from  (@seq_names) = keys %seq_input;
       my ($over_write_sso_by_age, @temp, $existing_sso, $out_gz_name,
           $over_write_msp_by_age, $existing_msp, $out_msp_file, $seq_name);
       $seq_name=$seq_names[$j];
       $each_seq_fasta="$seq_name\.fa";
       $out_msp_file="$seq_name\.msp";
       $out_gz_name="$seq_name\.msp\.gz";
       $out_msso_file="$seq_name\.msso";

       &die_if_file_not_present($each_seq_fasta);

       print "\n# (i) :-) Found $each_seq_fasta is searched against $sequence_DB\n";
       if($algorithm=~/fasta/){       $out_sso_file="$seq_name\.fsso";
       }elsif($algorithm=~/ssearch/){ $out_sso_file="$seq_name\.ssso"; }
       $out_sso_gz_name="$out_sso_name\.gz";

       if(-s $out_sso_file){ $existing_sso=$out_sso_file }
       elsif(-s $out_sso_gz_name){ $existing_sso=$out_sso_gz_name }
       if(-s $out_msp_file){ $existing_msp=$out_msp_file }
       elsif(-s $out_gz_name){ $existing_msp=$out_gz_name }
       if(  (localtime(time- (stat($existing_sso))[9]))[3] > $age_in_days_of_out_file ){
            $over_write_sso_by_age='o';
       }
       if(  (localtime(time- (stat($existing_msp))[9]))[3] > $age_in_days_of_out_file ){
            $over_write_msp_by_age='o';
       }

       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       #  To check if the target seq DB is in ../
       #________________________________________________
       if(-s $sequence_DB){ print "\n# (i) \$sequence_DB $sequence_DB exists, Good\n";
       }elsif( -s "../$sequence_DB"){ $sequence_DB="../$sequence_DB";
       }elsif( -s "../../$sequence_DB"){ $sequence_DB="../../$sequence_DB"; }

       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # If MSP file already exist
       #_____________________________________________________________
       if( -s $out_msp_file and !$over_write_msp_by_age and !$over_write ){
            print "\n# (i) File: $out_msp_file exists, skipping, to overwrite use \'o\' opt or set days";
            push(@final_out, $out_msp_file);
       }else{  ## -E is for e value cutoff. -b is for num of seq fetched
           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~~~~~~~~~~~~~~
           #  K-tuple is 1 by default. If xxxx.sso exsts, skip running fasta or ssearch
           #________________________________________________________________________________
           if(-s $out_sso_file and !$over_write ){ ## If SSO is already present, JUST READ IT!
                print "\n# (i) Just opening existing $out_sso_file $out_sso_file $out_msp_file $over_write_msp_by_age $over_write\n";
                open(SSO_ALREADY, "$out_sso_file");
                @temp=<SSO_ALREADY>;
                print "\n# (i) \@temp has ", scalar(@temp), " lines\n";
                close(SSO_ALREADY);
                &compress_files_by_gzip($out_sso_file);
           }else{ ## Run FASTA HERE
                print "\n# (i) Running \"run_fasta_sequence_search\" ";
                $gzipped_msp_file=${&run_fasta_sequence_search(
                                   "a=$algorithm",
                                   "O=$out_msp_file", "$create_sso",
                                   "File=$each_seq_fasta", "e=$E_val",
                                   "DB=$sequence_DB", "k=$k_tuple", "$machine_readable")};
                push(@final_out, $gzipped_msp_file) if -s $gzipped_msp_file ;
                unlink($each_seq_fasta) if $Lean_output;
           }
       }
       if($machine_readable and $create_sso and -s $out_sso_file){ &cp($out_sso_file, $out_msso_file); }
   } # end of for($j=0; $j< @seq_names; $j++){
   return(\@final_out);
   EXIT:

} # do_sequence_search

#______________________________________________________________________________
# Title     : compress_files_by_gzip
# Usage     : &compress_files_by_gzip('file_name_to_be_compressed');
# Function  :
# Example   :
# Keywords  :
# Options   :
# Author    : jong@salt2.med.harvard.edu
# Category  :
# Version   : 1.3
#------------------------------------------------------------------------------
sub compress_files_by_gzip{
    my($file, $input_file_name, @final_output_file_names, $out_file_gz_name);
    for $file(@_){
       if(ref($file) eq 'SCALAR'){
           $file=${$file};
       }
       $out_file_gz_name="$file\.gz";
       if(-s $out_file_gz_name){
           unlink ($out_file_gz_name); # remove existing $out_file_gz_name
       }
       if(@{&find_program_in_path('gzip')} > 0){
           system("gzip $file"); ## gzipping it
       }else{
           print "\n# (ERROR) gzip is not found here. Please put gzip compressing program in exe path\n";
           print "\n# (INFO) No compression will be made by gzip \n\n";
       }
       push( @final_output_file_names, $out_file_gz_name) if -s $out_file_gz_name;
       if(@final_output_file_names < 1){ print "\n# (ERROR) No gz files were written, Error??\n"; }
    }
    if(@final_output_file_names> 1){ return(\@final_output_file_names);
    }else{  return(\$final_output_file_names[0]);   }
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
#
#             'e=xxxx' for filtering input files by extension xxxx
#
# Returns   : Filenames in a reference of array
#             and input files in an array (file1, file2)=@{&parse_arguments};
# Argument  : uses @ARGV
# Version   : 2.0
#--------------------------------------------------------------------
sub parse_arguments{
     my( $c, $d, $f, $arg_num, $option_table_seen, $n, $option_filtered,
           $option_table_example, $input_line, @input_files,
           $extension);

     &import_ENV_vars; # this enables  $PDB40D automatically assigns the shell var value of $ENV{'PDB40D'}

     #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
     #   Checks if there were arguments
     #_______________________________________________________
     if( @ARGV < 1 ){ #<-- If Argument is not given at prompt
          for(@_){
                 if($_ eq '0'){ ## this means, parsearguments do not need any argus. So no need to parse anything
                         last;
                 }else{
                         print "\n \"$0\" requires at least one Argument, suiciding.\n\n";
                         print chr(7); #<-- This is beeping
                         print "  To get help type \"$0  h\"\n\n\n ";
                         exit;
                 }
          }
     }
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     #   Some DEFAULT $debug variables for debugging purposes
     #__________________________________________________________
     &set_debug_option;

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     #  If there is only one prompt arg. and is [-]*[hH][elp]*, it calls
     #   &show_default_help and exits
     #___________________________________________________________________
     if( ( @ARGV == 1 ) && ($ARGV[0] =~ /^[\-]*[hH\?][elp ]*$/) ){
           &show_default_help;
           exit;
     }

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     #  Checking some input options like 'e=txt' for extension filtering
     #_____________________________________________________________________
     for($i=0; $i< @_; $i++){
         if($_[$i]=~/e=(\S+)/){
             push(@extension, $1);
             }elsif($_[$i]=~/^[\-]?r$/){ ## reversing file order in return stage
             $reverse_out_file_name_order='r';
         }
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
     if($HELP == 1){ &show_default_help }
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # reverse order if 'r' opt is set
     #_____________________________________________
     if($reverse_out_file_name_order){
         @input_files=reverse(@input_files);
     }
     return(\@input_files);
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
# Version   : 2.7
#--------------------------------------------------------------------
sub assign_options_to_variables{
  my($i, %vars, $j, $op, $z, $n, $symb, $value, $var, %val, @val, $ARG_REG,
	 $option_table_example, @input_options, $first_border_and_title, $sym, @arg);

  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #      Defining small variables for option table reading
  #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  my($g)='gets'; my($if)='if';
  my($is)='is';  my(@input_files);
  my($o)='or';   my(@arguments) = sort @ARGV;

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #  Assigning global arguments(@num_opt, %vars) variables
  #_______________________________________________________________
  for($i=0; $i< @arguments; $i++){
	 if(($arguments[$i]=~/^(\-?\d+[\.\d+]?)$/)&&   ### it mustn't be a file
              ( !(-f $arguments[$i]) ) ){                ### getting NUM opt
	      push(@num_opt, $1);
	 }elsif( $arguments[$i]=~/^(\S+)=(\S+) *$/){
                $vars{$1}=$2;
		# print "\n# $1 $2\n";
	 }
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Some DEFAULT $debug variables for debugging purposes
  #_________________________________________________________
  &set_debug_option;

  #"""""""""""""""""""""""""""""""""""""""""""""""""""""""
  #   The main processing of self
  #_______________________________________________________
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
                  if($op=~/^([\$\@\%])([\w\-]+) {0,30}[=|$g|$is] *[\$\@\%]*([\- \w\.\d]+) *[bB]y +$ARG_REG/){
                       ## $sym   $var        becomes          a [$a...]       by       a -a -A
                       my $sym = $1;  #### The symbols like ($, @, %), '$' in the above.
                       my $var = $2;  #### Actual variable name 'var' from $var, 'av_sc_segment' in the above.
                       my $val = $3;  #### The becoming value  first 'a' in the above.
                       my @arg = ($4, $5, $6, $7, $8);  ## The alternative prompt arguments, second 'a' in the above..
                       print "\n $sym $var $val \n" if $debug==1;
                       print "\n \@arg are @arg \n" if $debug==1;

                       #""""""""""""""""""""""""""""""""""""""""""""""""""""
                       #  Going through the PROMPT args.
                       #____________________________________________________
                       #print "\n# \@arguments are @arguments\n";
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
                                      ( $arguments[$z] =~ /^${1}= *([\w\.*\-*\/]+) *$/) &&
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
# Title     : get_base_names
# Usage     : $base =${&get_base_names(\$file_name)};
#             :   or @bases = &get_base_names(\@files);  # <-- uses `pwd` for abs directory
# Function  : produces the file base name(eg, "evalign"  out of "evalign.pl" ).
#              when xxxx.xx.gz form file is given, it removes gz as well
#
# Example   : $base => 'test'  with 'test.txt' or '/home/dir/of/mine/text.txt'
# Warning   :
# Keywords  : get_base_name{, base_name, file_base_name ,  get_file_base_name
#             get_basename, basename, get_root_name, base , root, get_file_root
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

#______________________________________________________________
# Title     : write_msp_files
# Usage     : &write_msp_files(\%in1, \%in2, ['s'], [$filename],,)
# Function  : Writes input which is already in msp file format to
#              files either the name is given or generated
#              If more than one ref of hash is given, this will
#              concatenate all the hashes to one big one to
#              make one file.
#             When NO output xxx.msp file name is given, it creates
#              with the query sequence name.
# Example   :  &write_msp_files(@sso, 's', $out_file);
# Warning   : When NO output xxx.msp file name is given, it creates
#              with the query sequence name.
# Keywords  : write_msp,
# Options   : _  for debugging.
#             #  for debugging.
#             s  for each single file output for each hash input
#      filename  for putting output to the specified filename, should be xxx.msp
#
# Returns   : if 's' option is set, it will make say,
#               HI001.msp HI002.msp HI003.msp  rather than
#
#               HI001HI002HI003.msp
#  eg of one output(single file case)
#
#   1027     0.0     1     154   HI0004     1     154   HI0004
#   40       0.0     84    132   HI0004     63    108   HI0001
#   31       0.0     79    84    HI0004     98    103   HI0003
#
# Category  :
# Version   : 2.8
#--------------------------------------------------------------
sub write_msp_files{
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

    my ($out_msp_file, $add_range, @final_out, $msp_file_out,
         @keys, $N, $temp_1, %hash, $query_seq_name, $single_out_opt);

    if($char_opt=~/r/){ $add_range      ='r' };
    if($char_opt=~/s/){ $single_out_opt ='s' };
    if(@file == 1){ $out_msp_file=$file[0]; $single_out_opt='' } # s is for single file output

    if($single_out_opt eq 's'){ #~~~~~~~~~~~` single files output option WITHOUT a given outfilename
         $msp_file_out='default_single_out.msp';
         for($i=0; $i< @hash; $i++){
             my %hash=%{$hash[$i]};
             my @keys =sort keys %hash;

             #------------------ Writing the first line ---------------------------
             for($j=0; $j< @keys; $j++){
                  if($keys[$j]=~/(\S+)_\d+\-\d+/){ $N = $1 }else{ $N = $keys[$j] }
                  if($hash{$keys[$j]}=~/ +$N[\_\d+\-\d+]* +\d+ +\d+ +$N[\_\d+\-\d+]*/){
                       open(MSP_FILE, ">$msp_file_out") ||
                           die "# write_msp_files: I can not create $msp_file_out, check permission\n";
                       chomp( $hash{$keys[$j]} ); ## precaution
                       print MSP_FILE "# (H) $0 write_msp_files: $keys[$j]\n";
                       print MSP_FILE $hash{$keys[$j]}, "\n";
                       splice(@keys, $j, 1);
                       $j--; last;
                  }
             }
             #------------- Writing the rest of the lines ____________________
             for($j=0; $j< @keys; $j++){
                  chomp( $hash{$keys[$j]} );
                  print MSP_FILE $hash{$keys[$j]}, "\n";
             }
             print MSP_FILE "\n";
         }
         if(-s $msp_file_out){
              print "\n# write_msp_files: $msp_file_out is written \n";
         }else{
              print "\n# Error, write_msp_files\n"; exit
         }
         push(@final_out, $msp_file_out);
         close(MSP_FILE);
         return(\@final_out);
    }else{
         #~~~~~~~~~~~~~ DEfault ~~~~~~~~~~~~~~~~~~
         #  When output file name was given!
         #________________________________________
         if(@file==1){
             my($temp_1);
             open(MSP_FILE, ">$out_msp_file") ||  die "# write_msp_files: I can not create $out_msp_file, check permission\n";
             print MSP_FILE "# (H) $0 write_msp_files: @file\n";
             for($i=0; $i< @hash; $i++){
                 my %hash=%{$hash[$i]};
                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 # Sorting %hash values by the second column(Evalue)
                 #_______________________________________________________
                 @keys= map {$_->[0]} sort { $a->[1] <=> $b->[1] } map { $hash{$_}=~/^ *\S+[\t ]+(\S+)[\t ]+/ and [$_, $1] } keys %hash;

                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 # for Final output
                 #_____________________________
                 push(@final_out, $out_msp_file);

                 #--------- Writing the first line only --------------
                 for($j=0; $j< @keys; $j++){
                      if($keys[$j]=~/(\S+)_\d+\-\d+$/){ $N = $1 }else{ $N = $keys[$j] }

                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      # Following is to put the self match on top of the list
                      #________________________________________________________
                      if($hash{$keys[$j]}=~/ +$N[\_\d+\-\d+]* +\d+ +\d+ +$N[\_\d+\-\d+]*/){
                          $temp_1=$keys[0]; $keys[0]=$keys[$j]; $keys[$j]=$temp_1;
                      }
                 }
                 for($j=0; $j< @keys; $j++){
                      chomp($hash{$keys[$j]});
                      print MSP_FILE $hash{$keys[$j]}, "\n";
                 }
                 print MSP_FILE "\n";
             }
             print MSP_FILE "\n";
             close(MSP_FILE);
             if(-s $out_msp_file and $out_msp_file !~/^ *\.msp$/){
                 print "\n# write_msp_files: $out_msp_file is written\n" if(-s $out_msp_file);
             }else{
                 print "\n# write_msp_files: ERROR. Either $out_msp_file is empty or \".msp\" is written\n";
             }
         }else{
             for($i=0; $i< @hash; $i++){
                 my %hash=%{$hash[$i]};
                 my @keys =sort keys %hash;
                 ($query_seq_name)=$hash{$keys[0]}=~/\S+ +\d+ +\d+ +(\S+) +\d+ +\d+ +\S+/;
                 $msp_file_out="$query_seq_name\.msp";
                 open(MSP_FILE, ">$msp_file_out") or die "\n# write_msp_files: Failed to open $msp_file_out\n";
                 print MSP_FILE "# (H) $0 write_msp_files: $query_seq_name\n";

                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 # for Final output
                 #_____________________________
                 push(@final_out, $msp_file_out);
                 #~~~~~~~ Writing the first line only ~~~~~~~~~~~~~~~~~~
                 for($j=0; $j< @keys; $j++){
                     if($keys[$j]=~/(\S+)_\d+\-\d+$/){ $N = $1 }else{ $N = $keys[$j] }
                     if($hash{$keys[$j]}=~/ +$N[\_\d+\-\d+]* +\d+ +\d+ +$N[\_\d+\-\d+]*/){
                      $keys[0]=$temp_1; $keys[0]=$keys[$j]; $keys[$j]=$temp_1;
                     }
                 }
                 for($j=0; $j< @keys; $j++){
                     chomp($hash{$keys[$j]});
                     print MSP_FILE $hash{$keys[$j]}, "\n";
                 }
                 print MSP_FILE "\n";
             }
             print MSP_FILE "\n";
             if(-s $out_msp_file and $out_msp_file !~/^ *\.msp$/){
                  print "\n# write_msp_files: $out_msp_file is written\n" if(-s $out_msp_file);
             }else{
                  print "\n# write_msp_files: ERROR. Either $out_msp_file is empty or only \".msp\" is written\n";
             }
             close MSP_FILE;
         }
     }
     if(@final_out ==1){ return( \$final_out[0] ); }else{
         return(\@final_out);
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
# Title     : show_default_help
# Usage     : &show_default_help2;  usually with 'parse_arguments' sub.
# Function  : Prints usage information and others when invoked. You need to have
#             sections like this explanation box in your perl code. When invoked,
#             show_default_help routine reads the running perl code (SELF READING) and
#             displays what you have typed in this box.
#             After one entry names like # Function :, the following lines without
#             entry name (like this very line) are attached to the previous entry.
#             In this example, to # Function : entry.
# Example   : &show_default_help2; &show_default_help2(\$arg_num_limit);   &show_default_help2( '3' );
#             1 scalar digit for the minimum number of arg (optional),
#             or its ref. If this defined, it will produce exit the program
#             telling the minimum arguments.
# Warning   : this uses format and references
# Keywords  : default_help, show_headbox_help, show_help
# Options   :
# Returns   : formated information
# Argument  :
# Category  :
# Version   : 3.5
#--------------------------------------------------------------------
sub show_default_help{
     my($i, $perl_dir, $arg_num_limit, $head ,$arg_num_limit, $key_press, $e, $help_item,
       @entries, @entries_I_want_write, $option_tb_found, $extension, $logname, $tmp );
     $logname=getlogin();
     my($pwd)=`pwd`;
     my($date)=`date`;
     chomp($date,$pwd);
     my($not_provided)="--- not provided ---\n";
     my($file_to_read) = $0;

     for($i=0; $i < @_; $i ++){
       if((ref($_[$i]) eq 'SCALAR')&&(${$_[$i]} =~ /^\d$/)){
               $arg_num_limit = ${$_[$i]};
       }elsif( (!(ref($_[$i]))) && ($_[$i] =~ /^\d$/)){
               $arg_num_limit = $_[$i];     }
     }
     my %entries = %{&read_head_box(\$file_to_read )};
     if($option_tb_found ==1){
        @option_tb=@{&read_option_table(\$file_to_read)};
     }

     @entries = keys %entries;
     for $help_item (@entries){
        if( $entries{$help_item}=~/^[\W]*$/  and  !defined($entries{$help_item}) ){
            $entries{$help_item}= $not_provided;
        }
     }

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
     #  Writing the format <<<<<<<<<<<
     #________________________________________
     $~ =HEADER_HELP;     write;   ## <<--  $~ is the selection operator
     $~ =DEFAULT_HELP_FORM;

     @entries_I_want_write=sort keys %entries;

     for ( @entries_I_want_write ){  write  }

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
		     \/_/\/_/    \/_____/   \/______/   \/_/ V 3.5`
_____________________________________________________________________
.

format DEFAULT_HELP_FORM =
 @<<<<<<<<<: @*
 $_,        $entries{$_}
.

}




#________________________________________________________________________
# Title     : open_fasta_files
# Usage     : %fasta_seq=%{&open_fasta_files($fasta_file, ['MJ0084'])};
#             if you put additional seq name as MJ0084 it will
#             fetch that sequence only in the database file.
#
#             %out=%{&open_fasta_files(@ARGV, \%index)};
#               while  %index has (seq indexpos seq2 indexpos2,,,)
#               In this case, the fasta file should have xxxx.fa format
#
# Function  : open fasta files and put sequences in a hash
#              If hash(es) is put which has sequence names and seek position
#              of the index file, it searches the input FASTA file to
#              fetch at that seek position. This is useful for Big fasta DBs
#             If the seq name has ranges like  XXXXXX_1-30, it will only
#              return 1-30 of XXXXXX sequence.
#
#             FASTA sequence file format is like this;
#
#             > 1st-seq
#             ABCDEFGHIJKLMOPABCDEFGHIJKLMOPABCDEFGHIJKLMOPABCDEFG
#             > 2nd.sequ
#             ABCDEFGHIJKLMOYYUIUUIUIYIKLMOPABCDEFGHIJKLMOPABCDEFG
#             >owl|P04439|1A03_HUMAN HLA CLASS I HISTOCOMPATIBILITY ANTIGEN, A-3 ALPHA CHAIN PRECURSOR....
#             MARGDQAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDT
#
#             This can also return the sizes of sequences rather than seqs.
#
#             This ignores any dup entrynames coming later.
#
# Example   : %out = %{&open_fasta_files(@ARGV)};
#             %out2=%{&open_fasta_files('seq.fa', \%index)};
#             %out3=%{&open_fasta_files('seq.fa', \%range)};
#             %seq=%{&open_fasta_files($PDB40_FASTA, \@seq_to_fetch)};
#
#             while @ARGV at prompt was: 'GMJ.pep MJ0084'
#
# Keywords  : open_fasta, open_fa_files, open_FASTA_files,
# Options   : Seq name to fetch the specified seq only.
#             as open_fasta_files.pl MY_SEQ_NAME Swissprot.fasta
#            -d  for giving back desc as well as the name. so it
#                gives  'HI0002 This is the description part'
#                as the key
#             If you put hash which is like('seq_name', ['20-30', '30-44',..])
#              it will produce hash which has got:
#              ( seq_name_20-30 'asdfasdfasdfasdfasd',
#                seq_name_30-44 'kljkljkjkjljkjljkll',
#                ....           .... )
#            -s for returning sequence size only
# Version   : 4.0
#--------------------------------------------------------------------
sub open_fasta_files{
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

   my (%sequence, %HASH, @Keys, $seq_found1, $S_start, $S_end, $seq_found,
	   $present_seq, @seq_Names, %Sizes, $bare_seq_name, $fasta_seq_idx_file,
	   %seq_fragments);

   if(@file<1){
	  print "\n \@file has less than 1 elem. There is no fileinput for open_fasta_files\n";
	  exit
   }

   print "\n# \@file are @file\n" if $verbose;

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  (1) When perl file INDEX pos info is given in hash, this speeds up
   #_________________________________________________________________________________
   for($d=0; $d < @hash; $d++){
	   my ($sequence, $NAME, $range_start, $range_leng);
	   %HASH=%{$hash[$d]};
	   my @Keys=keys %HASH;  ## <<< NOTE it is @Keys, not @keys
	   for($f=0; $f< @file; $f++){
		  #====== It must be xxxx.fa format =======
		  unless($file[$f]=~/\S\.fa[sta]?$/){
			  print "\n# open_fasta_files: \$file\[\$f\] does not have fasta extension, skipping"; next; }
		  open(FASTA, $file[$f]);
		  F0: for($e=0; $e< @Keys; $e++){
			 my $sequence;
			 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			 # When seq name has range attachment, it handles
			 #________________________________________________
			 if($Keys[$e]=~/^(\S+)_(\d+)\-(\d+)/){
				 $NAME=$1;
				 $range_start=$2-1;    ## to fit in substr function
				 $range_leng =$3-$2+1; ## to fit in substr
			 }else{
			     $NAME=$Keys[$e];
			 }
			 if($HASH{$Keys[$e]}=~/^(\d+)$/){
				 splice(@hash, $d, 1);
				 $d--;
				 splice(@file, $f, 1);
				 $f--;
				 seek(FASTA, $1-220, 0);  # -220 is necessary
				 while(<FASTA>){
					 if( /^\> *$NAME/  or
						 /^\> *owl\|\S+\|$NAME/){  # to handle ">owl|P04439|1A03_HUMAN HLA CLASS I HISTOCOMPATIBILITY
					        $seq_found1=1;
					 }elsif(/^(\w+)$/ and $seq_found1==1){	 $sequence .=$1;
					 }elsif(/^\> *\S+/ and $seq_found1==1){
						  #======= When range is defined, take only the ranged part==================
						  if($range_start =~/\d+/){
							  $sequence{$Keys[$e]}=substr($sequence, $range_start, $range_leng);
						  }else{	 $sequence{$Keys[$e]}=$sequence; }
						  $range_start='';
						  $sequence='';
						  $seq_found1=0; next F0;
					 }
				 }
			  }
		  }
	  }
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   # (2) opening FASTA files (NORmal, when no perl index pos number is given)
   #_______________________________________________________________________
   for($i=0; $i< @file; $i++){
	   unless(-s $file[$i]){ next; } ## this is essential as handle_arguments has a problem
	   my($entry_found, $name, $matched);
	   my($input_file) = ${$file[$i]} || $file[$i];

	   if($debug eq 1){ print "\n open_fasta_files: Inputfile is $input_file\n" };
	   unless (-e $input_file){
			print chr(7);
			print "\n\n\t This is sub open_fas_files in $0  \n\n";
			print "\t Fatal: The input file $input_file is not in the directory \n";
	   }
	   open(FILE_1,"$input_file");
	   if(@hash >=1){  ## if seq names are given in hash
		   for($h=0; $h< @hash; $h++){
			  @string=(@string, keys %{$hash[$h]});
		   }
	   }
	   @string=sort @string;
	   $num_of_seq_to_fetch=@string;
	   if(@string > 0){
		   print "\n# open_fasta_files(normal fasta fetch): \$num_of_seq_to_fetch is $num_of_seq_to_fetch\n";
	   }

	   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	   #  (2.1) when seq to fetch is given by \@sequences  format
	   #_______________________________________________________________________
	   if( @_ > 1  and  @string > 0 ){
		   print "\n#  open_fasta_files is fetching sequences from \$input_file= $input_file\n";
		   %sequence=%{&fetch_sequence_from_db($input_file, \@string)};
		   print "\n# $fasta_seq_idx_file file is made by open_fasta_files(fetch_sequence_from_db), you may remove it\n";
	   }
	   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	   #  (2.2) When seq names NOT given, fetches all (THE DEFAULT)
	   #____________________________________________________________
	   else{
		 while(<FILE_1>){                # file1 needs to be xxxx.fasta for the moment, automatic later
			if(/^> *gi\|\d+\|\S+\|(\S+)\|.*/){  ## for >gi|1669546|dbj|D84107|D84107 Human mRNA for Werner syndrome-1/type 1, complete cds
				 if($char_opt=~/[\-]?d/i){  # To add the description
					 $name=$_;  # entire line becomes the name of the seque.
				 }else{
					 if( $sequence{$1} ){
						 #------- To avoid identical entry reading repeatedly -----
						 print "\n# I am open_fasta_files: $1 seems to be the same as previous entry, ERROR??\n";
						 $entry_found=0;
					 }else{      $name=$1;   $entry_found=1;     }
				 }
            }elsif(/^> *owl\|\S*\|(\S+)/ or /^> *pir\|\S*\|(\S+)/){  ## for ">owl|P04439|1A03_HUMAN HLA CLASS I HISTOCOMPATIB
				 if($char_opt=~/[\-]?d/i){  # To add the description
					 $name=$_;  # entire line becomes the name of the seque.
				 }else{
					 if( $sequence{$1} ){
						 #------- To avoid identical entry reading repeatedly -----
						 print "\n# I am open_fasta_files: $1 seems to be the same as previous entry, ERROR??\n";
						 $entry_found=0;
					 }else{      $name=$1;   $entry_found=1;     }
				 }
			}elsif(/^> {0,5}([\w\-\.]+) *.*$/){
				 if($char_opt=~/[\-]?d/i){   $name=$_;  # To add the description
				 }else{
					 if( $sequence{$1} ){ # check if the entry already exists
						print "\n# $1 seems to be the same as previous entry, ERROR??\n";
						$entry_found=0;
					 }else{     $name=$1;   $entry_found=1;      }
				 }
			}elsif(/^([\w\.\- ]+)$/ and $entry_found == 1){
                 $matched=$1;    $matched=~s/ //g;
                 if($char_opt=~/s/){ # when SIZE(length of seq) return only option is set
                     $Sizes{$name}+=length($matched);
                     #print $Sizes{$name}, "\n";
                 }else{
                     $sequence{$name}.= $matched if defined($name);
                 }
			}elsif(/^$/){  next;
			}else{  $entry_found=0;  } ## this is when rubbish is matched
		 }# end of while
	   }
	   close FILE_1;
   }


   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~~~~~~~~~~~~~~~~~~~~~`
   # (3) When ranges information is given(via \@range), seq in those ranges are returned
   #______________________________________________________________________________________
   if(defined(@range)){
	   %seq_fragments=%{&get_seq_fragments(\%sequence, \@range)};
	   return(\%seq_fragments);
   }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   # (4) When only size is asked with -s option
   #_____________________________________________________________________________
   elsif($char_opt=~/s/){ # when SIZE(length of seq) return only option is set
	   #@seq_Names=keys %sequence;
	   #for($i=0; $i<@seq_Names; $i++){
	   #	  $Sizes{$seq_Names[$i]}=length($sequence{$seq_Names[$i]});
	   # }
	   return(\%Sizes);
   }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
   # (5) when hash which has range info is given(@range should not be defined)
   #_____________________________________________________________________________
   elsif(@hash >=1){
	   for($h=0; $h< @hash; $h++){
		   my %hash=%{$hash[$h]};
		   my @Keys=keys %hash;
		   for($k=0; $k<@Keys; $k++){
			   if(defined($hash{$Keys[$k]})){
				  ($S_start, $S_end)=$hash{$Keys[$k]}=~/(\d+)\-(\d+)/;
				  $sequence{$Keys[$k]}=substr($sequence{$Keys[$k]}, ($S_start-1), ($S_end-$S_start));
			   }
		   }
	   }
	   return(\%sequence);
   }else{
	   return(\%sequence);
   }
}
#________________________________________________________________________________
# Title     : tempname
# Usage     : $tmp=&tempname;
# Function  : Returns a unique temporary filename.
#             Reasonably robust but not completely immune to race conditions
#             with other processes simultaneously requesting a tempname.
# Example   :
# Keywords  :
# Options   :
# Version   : 1.0
#--------------------------------------------------------------------------------
sub tempname{
   foreach $suffix (0..99) {
 	if (! (-e "tmpxx$suffix")) {
		open(TMP,">tmpxx$suffix"); # Stamp it to reserve it.
		close(TMP);
		return "tmpxx$suffix";
	}
   }
}
#__________________________________________________________________________
# Title     : if_file_older_than_x_days
# Usage     : if( ${&if_file_older_than_x_days($ARGV[0], $days)} > 0){
# Function  : checks the date of last modi of file given and compares with
#             present time. Substracts diff and returns the actual diff days.
# Example   :
# Keywords  : how_old_file, how_old, is_file_older_than_x_days, file_age,
#             file_age_in_days, if_older_than_x_days,
# Options   :
# Returns   : the actual days older, so NON-ZERO, otherwise, 0
# Version   : 1.3
#----------------------------------------------------------------------------
sub if_file_older_than_x_days{
	if(@_ < 2){ print "\n# FATAL: if_file_older_than_x_days needs 2 args\n"; exit; }
	my $file=${$_[0]} || $_[0];
	my $days=${$_[1]} || $_[1];
	my ($new_idx_file, $how_old_days);
	unless(-s $file){
	    print "\n# FATAL, nearly!: if_file_older_than_x_days: $file does NOT exist !\n";
		$new_idx_file=${&make_seq_index_file($file)};
		print "        if_file_older_than_x_days called make_seq_index_file to make $new_idx_file\n";
        $how_old_days=(localtime(time- (stat($new_idx_file))[9]))[3];
	}else{
        $how_old_days=(localtime(time- (stat($file))[9]))[3];
    }
	if($how_old_days > $days){
		print "\n# if_file_older_than_x_days: $file is older than $days\n";
		return(\$days);
	}else{
		print "\n# if_file_older_than_x_days: $file is NOT older than $days\n";
		return(0);
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
#
# Version   : 2.8
#--------------------------------------------------------------------
sub read_file_names_only{
  my($in_dir, $i, $j, $x, $k, $dir, @final_files, @possible_dirs, $sort_opt, $ext, @extensions,
      $path_include, @in, $glob_given, @files_globed, @in_dir, $pwd, $extension_given,
      %target_file_names, @target_file_names, @read_files);
  $pwd=`pwd`; chomp($pwd);
  $in_dir=$pwd;
  @in=@_;

  print "\n# read_file_names_only: input are \"@in\"";

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  Directory entry and opts detection
  #_________________________________________
  for($k=0; $k < @in; $k++){
	 if   ( $in[$k] eq '.'){ push(@in_dir,$pwd); splice(@in, $k, 1);  $k--; next }
     if( !(ref($in[$k]))){
	    print "\n# read_file_names_only: $in[$k] is not a reference";
        if($in[$k]=~/D=(\S+)/){
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
            print "\n# (E) read_file_names_only: $in[$k] is not a file, nor a dir, STRANGE?";
            print "\n#     Was it a file extension?\n";
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
  sort @final_files if $sort_opt == 1;
  return(\@final_files);
}
#______________________________________________________________________
# Title     : sort_string_by_length
# Usage     : @output = @{&sort_string_by_length(@any_input_strings, [-r], @more)};
# Function  : sorts strings in array according to their sizes
#             bigger comes first.
# Example   :
# Warning   :
# Keywords  : sort_array_by_length, sort_str_by_length, sort_array_string_by
#             sort_string_by_leng, sort_by_length, sort_by_leng,
#             sort_array_by_string_length, sort_array_elements_by_string_length
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


#_____________________________________________________________________________
# Title     : fetch_sequence_from_db
# Usage     : %sequence=%{&fetch_sequence_from_db($input_file, \@string)};
# Function  : accept seq names (with or without ranges like _10-111 )
#              and produces hash ref.
#             As an option, you can write(xxxx.fa) the sequences in pwd
#              with the file names with sequence names.
#             The default database used is FASTA format OWL database.
#              You can change this by S (for Swissprot either fasta
#              or full format), P for PDB fasta format data.
#             If you give the path name of DB, it will look for the
#              DB given.
#
#             This automatically checks sequence family number as
#               in >d1bpi___7.6.1
#               and attaches the number in final %sequence output
#
# Example   : %seq=%{&fetch_sequence_from_db(\@input, seq.fa, seq.fa.idx)};
#              while @input=qw( 11S3_HELAN_11-31 A1AB_CANFA A1AT_PIG )
# Keywords  : fetch_seq_from_db, fetch_sequence_from_database
# Options   : _  or #  for debugging.
#     w       for write fasta file
#     s=      for putting source DB file name manually
#     d=p100  for PDB100 fasta database from ENV
#     d=p40   for PDB40  fasta database from ENV
#     d=p     for PDB database (usually p100) from ENV
#     d=s     for Swissprot database from ENV
#     d=o     for OWL database from ENV
#     i=      for index filename. If not specified, this looks for it in the same dir as fast     ˜
#     t=      for msp_threshold
#  msp_threshold=0.0005  # when MSP file is given as input for getting seq names
#
# Returns   : ref of hash
# Argument  : gets names of sequences
#             eg) \@array, \%hash, \$seq, while @array=(seq1, seq2), $seq='seq1 seq1'
#                                               %hash=(seq1, xxxx, seq2, yyyy);
#
# Version   : 3.5
#------------------------------------------------------------------------------
sub fetch_sequence_from_db{
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

        my(@DATABASE, @INDEX_FILE, %sequence, %seq_with_index, @input_seq_names,
	   %long_index, @Keys, $R_start, $NAME, $R_leng, $found_seq_count,
	   $eval_blastpgp_thresh, $seq_found1, $sequence, @keys, $index_file,
	   $source_DB_file, $matching_seq, $match_start, $match_stop);

        $eval_blastpgp_thresh=0.001; # default

        if($vars{'msp_threshold'}=~/(\S+)/ or $vars{'mt'}=~/(\S+)/ or $vars{'t'}=~/(\S+)/ ){
              $eval_blastpgp_thresh=$1;
              print "\n# (i) YOU have set the \$eval_blastpgp_thresh $eval_blastpgp_thresh\n\n";
        }
        if($vars{'s'}=~/(\S+)/ or $vars{'DB'}=~/(\S+)/ ){
              $source_DB_file=$1; push(@DATABASE, $source_DB_file);
        }

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# getting input seq names from all sources
	#________________________________________________________
	for(0..@hash){ # from the given hash (keys names)
	     push(@input_seq_names, keys %{$hash[$_]} );
	}
	for(0..@raw_string){  ## from given sequence names
	     push(@input_seq_names, split(/ +/, $raw_string[$_]) );
	}
        for($i=0; $i<@file; $i++){  ## From MSP file input (get only MATCHED sequences)
             if($file[$i]=~/\.msp/){
                  print "\n# (i) MSP file input is detected !\n";
                  my ($seq_with_range);
                  open(MSP, $file[$i]);
                  while(<MSP>){
                       if(/^ *\S+ +(\S+) +\S+ +\S+ +\S+ +\S+ +\S+ +\S+ +(nr_\S+_DROME_\S+) +/){
                            push(@input_seq_names, $2) if $1 < $eval_blastpgp_thresh;
                       }elsif(/^[\t ]*\S+[\t ]+(\S+)[\t ]+\S*[\t ]*(\d+)[\t ]+(\d+)[\t ]+\S+[\t ]+(\d+)[\t ]+(\d+)[\t ]+(\S+) */){
                            $matching_seq=$6; $match_start=$4; $match_stop=$5;
                            $evalue=$1;
                            if($matching_seq=~/^(\S+)_\d+\-\d+/){
                                 $seq_with_range=$matching_seq;
                            }else{
                                 $seq_with_range="$matching_seq\_$match_start\-$match_stop";
                            }
                            push(@input_seq_names, $seq_with_range) if $evalue < $eval_blastpgp_thresh;
                            print "\n# (i) pushing $seq_with_range" if $verbose;
                       }
                  }
                  close(MSP);
                  splice(@file, $i, 1);
                  $i--;
             }
        }

	print "\n# (1) fetch_sequence_from_db: \@raw_string has: ", scalar(@raw_string), " elements";
	print "\n# (2) fetch_sequence_from_db: No. of seq to fetch is:",scalar(@input_seq_names);
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# Choose the DBs and INDEX for fetching sequences. All files input must be DATABAE or INDEXfile
	#___________________________________________
	if(@file > 0){
             for($i=0; $i< @file; $i++){
                if(-T $file[$i] and $file[$i]=~/\.fa[sta]?$/){      push(@DATABASE, $file[$i]);   next}
                elsif((-T $file[$i]) and ($file[$i]=~/\.seq$/)){    push(@DATABASE, $file[$i]);   next}
                elsif((-T $file[$i]) and ($file[$i]=~/\.dat$/)){    push(@DATABASE, $file[$i]);   next}
                elsif(-T $file[$i] and $file[$i]=~/\.idx$/){        push(@INDEX_FILE, $file[$i]); next }
                if($file[$i] !~/\.idx/ and -T "$file[$i]\.idx"){    push(@INDEX_FILE, "$file[$i]\.idx"); }
                else{
                    print "\n#  WARN:  fetch_sequence_from_db:
                    You put a file-name-like which is not a fasta DB. Error. I am removing $file[$i]";
                    splice(@file, $i, 1);
                    $i--;
                }
             }
	}

	if($vars{'d'}=~/^p[100]*$/){
	   if( -T  $ENV{'PDB_FASTA'} ){             push(@DATABASE,   $ENV{'PDB_FASTA'} );     }
	   elsif(  -T $ENV{'PDB_SEQ_FASTA'} ){      push(@DATABASE,   $ENV{'PDB_SEQ_FASTA'}  ); }
	   elsif(  -T $ENV{'PDB100_FASTA'} ){       push(@DATABASE,   $ENV{'PDB100_FASTA'} ); }
	   if(  -T $ENV{'PDB_FASTA_INDEX'} ){       push(@INDEX_FILE, $ENV{'PDB_FASTA_INDEX'} ); }
	}elsif( $vars{'d'}=~/^p\d+d$/ ){
	   if(  -T $ENV{'PDB100D_FASTA'} ){         push(@DATABASE,   $ENV{'PDB100D_FASTA'});     }
	   elsif(  -T $ENV{'PDBD100_FASTA'}  ){     push(@DATABASE,   $ENV{'PDBD100_FASTA'}); }
	   elsif(  -T $ENV{'PDB100D_SEQ_FASTA'}  ){ push(@DATABASE,   $ENV{'PDB100D_SEQ_FASTA'}); }
	   elsif(  -T $ENV{'PDBD100_SEQ_FASTA'}  ){ push(@DATABASE,   $ENV{'PDBD100_SEQ_FASTA'}); }
	   if(  -T $ENV{'PDB100D_SEQ_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'PDB100D_SEQ_FASTA_INDEX'}); }
	   elsif(  -T $ENV{'PDBD100_SEQ_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'PDBD100_SEQ_FASTA_INDEX'}); }
	}elsif( $vars{'d'}=~/^p40/ ){
	   if(  -T $ENV{'PDB40_FASTA'} ){          push(@DATABASE,   $ENV{'PDB40_FASTA'});     }
	   elsif(  -T $ENV{'PDB40_SEQ_FASTA'}  ){  push(@DATABASE,   $ENV{'PDB40_SEQ_FASTA'}); }
	   if(  -T $ENV{'PDB40_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'PDB40_FASTA_INDEX'}); }
	}elsif( $vars{'d'}=~/^p90/ ){
	   if(  -T $ENV{'PDB90_FASTA'}  ){         push(@DATABASE,   $ENV{'PDB90_FASTA'}    ); }
	   elsif(  -T $ENV{'PDB90_SEQ_FASTA'} ){   push(@DATABASE,   $ENV{'PDB90_SEQ_FASTA'}); }
	   if(  -T $ENV{'PDB90_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'PDB90_FASTA_INDEX'}); }
	}
	if( $vars{'d'}=~/^s *$/){
	   if(  -T $ENV{'SWISS_FASTA'} ){          push(@DATABASE,   $ENV{'SWISS_FASTA'});     }
	   elsif(  -T $ENV{'SWISS_SEQ_FASTA'} ){   push(@DATABASE,   $ENV{'SWISS_SEQ_FASTA'}); }
	   elsif(  -T $ENV{"SWISS_DIR\/seq.fa"} ){ push(@DATABASE,   $ENV{"SWISS_DIR\/seq.fa"}); }
	   if(  -T $ENV{'SWISS_FASTA_INDEX'} ){    push(@INDEX_FILE, $ENV{'SWISS_FASTA_INDEX'}); }
	   elsif(  -T $ENV{'SWINDEX'} ){           push(@INDEX_FILE, $ENV{'SWINDEX'}); }
	}
	elsif( $vars{'d'}=~/^o *$/){
		if(  -T $ENV{'OWL_FASTA'} ){            push(@DATABASE,   $ENV{'OWL_FASTA'});     }
		elsif(  -T $ENV{'OWL_SEQ_FASTA'} ){     push(@DATABASE,   $ENV{'OWL_SEQ_FASTA'}); }
		elsif(  -T $ENV{"OWL_DIR\/owl.fa"} ){   push(@DATABASE,   $ENV{"OWL_DIR\/owl.fa"}); }
		if(  -T $ENV{'OWL_FASTA_INDEX'} ){      push(@INDEX_FILE, $ENV{'OWL_FASTA_INDEX'}); }
		print "\n# Fetching sequences from OWL\n";
	}
	elsif( $vars{'d'}=~/^n *$/){
	   if(  -T $ENV{'NRDB_FASTA'} ){            push(@DATABASE,   $ENV{'NRDB_FASTA'});     }
 	   elsif(  -T $ENV{'NRDB_SEQ_FASTA'} ){     push(@DATABASE,   $ENV{'NRDB_SEQ_FASTA'}); }
			 if(  -T $ENV{'NRDB_FASTA_INDEX'} ){      push(@INDEX_FILE, $ENV{'NRDB_FASTA_INDEX'}); }
			 elsif(  -T $ENV{'NRDB_FASTA_IDX'} ){     push(@INDEX_FILE, $ENV{'NRDB_FASTA_IDX'}); }
	   print "\n# Fetching sequences from OWL\n";
	}
	elsif( $vars{'d'}=~/^\S+\.\S+$/ and -T $vars{'d'} ){ push(@DATABASE, $vars{'d'} );     }
	if( $vars{'i'}=~/\S+\.\S+$/ and -T $vars{'i'} ){ push(@INDEX_FILE, $vars{'i'} );   }
	if(@INDEX_FILE > 0 and @DATABASE > 0){
             if( ${&if_file_older_than_x_days("$DATABASE[0]\.idx", 5)} > 0 ){
                  $index_file=${&make_seq_index_file(\@DATABASE)};
                  push(@INDEX_FILE, $index_file);
             }elsif((-s "$DATABASE[0]\.idx") > 50){
                     push(@INDEX_FILE, "$DATABASE[0]\.idx");
             }else{
                     print "\n# (E) fetch_sequence_from_db: Some weird error in pushing \$index_file to \@INDEX_FILE\n"; exit;
             }
	}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#  Final check for ALL the inputs
	#___________________________________________________
	if( @DATABASE  < 1){ print "\n# fetch_sequence_from_db: DATABASE file no found. Error\n"; exit     }
	if( @INDEX_FILE < 1){
             print "\n# fetch_sequence_from_db: \@INDEX_FILE has less than 1 elem. Error\n";
             push(@INDEX_FILE, ${&make_seq_index_file(@DATABASE)});
             print "     fetch_sequence_from_db called make_seq_index_file to make @INDEX_FILE\n";
	}
 	if($debug==1){
             print "\n# DATABASE used     : @DATABASE";
             print "\n# INDEX_FILE used   : @INDEX_FILE";
             print "\n# input_seq_names   : @input_seq_names";
	}


        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ##  Now I have @DATABASE, @INDEX_FILE, @input_seq_names
        ##_______________________________________________________________

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	#  Reading in index file to get 'seq' 'seek pos' to make %seq_with_index
	#__________________________________________________________________________
	print "\n#  fetch_sequence_from_db: \@INDEX_FILE @INDEX_FILE, \@DATABASE :@DATABASE\n";
	for($i=0; $i< @INDEX_FILE; $i++){
	   open(INDEX, "$INDEX_FILE[$i]");
	   while(<INDEX>){ if(/(\S+) +(\S+)/){  $long_index{$1}=$2;  }  }

	   for($j =0; $j < @input_seq_names; $j++){

                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~``
                 #  If DATABASE has sequence names with ranges already index the seq with ranges
                 #____________________________________________________________________________________
                 if($input_seq_names[$j]=~/(\S+\_\d+\-\d+)$/ and $long_index{$1}){
                     $seq_with_index{$1}=$long_index{$1};

                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~``
                 #  If DATABASE has sequence names without ranges index the seq without ranges
                 #____________________________________________________________________________________
                 }elsif($input_seq_names[$j]=~/(\S+)\_\d+\-\d+$/ and $long_index{$1}){

                         $seq_with_index{$input_seq_names[$j]}=$long_index{$1}; # !!!! <--- This line is critical

                 }elsif($input_seq_names[$j]=~/(\S+)\_\d+\-\d+/ and $long_index{"$1\_"}){ # to handle Tim's new pdb100.fa files

                     $seq_with_index{$input_seq_names[$j]}=$long_index{"$1\_"};
                     print "\n# Warning: $1 (from $input_seq_names[$j]) matched with $1\_ in $INDEX_FILE[$i],
                                   I hope this is correct!!\n";
                 }
                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~``
                 #  If input_seq_name has SCOP superfamily numbers
                 #____________________________________________________________________________________
                 elsif($input_seq_names[$j]=~/^(\S+)\_(\d+\.\d+\.\d+)[\.\d+\.\d+]*/ and $long_index{$1}){

                         $seq_with_index{"$1\_$2"}=$long_index{$1}; # !!!! <--- This line is critical

                 }elsif($input_seq_names[$j]=~/\S/ and $long_index{$input_seq_names[$j]}){
                         $seq_with_index{$input_seq_names[$j]}=$long_index{$input_seq_names[$j]}
                 }else{
                         print "\n#  $input_seq_names[$j](with, without range) have NO corresponding index in $INDEX_FILE[$i], ERR";
                 }
	   }
	   close INDEX;
	   if ( scalar(keys %seq_with_index) < 1){
                  print "\n# fetch_sequence_from_db: \%seq_with_index is too small, ERROR?\n";
	   }
	}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
	#  Fetching sequences from DATABASE
	#_______________________________________________________________
	print "\n# fetch_sequence_from_db: Fetching seqs from @DATABASE with  @INDEX_FILE ";
	@Keys= sort {$seq_with_index{$a} <=> $seq_with_index{$b} } keys %seq_with_index;        ## <<< NOTE it is @Keys, not @keys
	print "\n# (3) fetch_sequence_from_db: No. of seq indexed is:", scalar(@Keys);

	for($f=0; $f< @DATABASE; $f++){

           open(DB_FASTA, $DATABASE[$f]);

	   F0: for($e=0; $e< @Keys; $e++){
		  my ($seq_found1, $super_fam_class, $NAME, $R_leng, $R_start, $sequence);
		  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  # When seq name has range attachment, it handles
		  #________________________________________________
		  if($Keys[$e]=~/(\S+)_(\d+)\-(\d+)$/){
		      $NAME=$1;
		      $R_start=$2-1;      ## to fit in substr function
		      $R_leng =$3-$2+1; ## to fit in substr
		      print "\n# (4) fetch_sequence_from_db: Sequences have ranges ($R_start-$R_leng) only (not superfamily numb.) \n";
		  }
		  elsif($Keys[$e]=~/(\S+)_(\d+)\-(\d+)\_(\d+\.\d+\.\d+)[\.\d+\.\d+]*/){
		      $NAME=$1;
		      $R_start=$2-1;      ## to fit in substr function
		      $R_leng =$3-$2+1; ## to fit in substr
		      $super_fam_class=$4;
		      print "\n# (4) fetch_sequence_from_db: Sequences have ranges and superfamily numb.\n";
		  }
		  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		  # When superfamily (scop) number is attached
		  #___________________________________________________
		  elsif($Keys[$e]=~/(\S+)\_(\d+\.\d+\.\d+)[\.\d+\.\d+]*/){
                      $NAME=$1;
		      $super_fam_class=$2;
		      print "\n# (4) fetch_sequence_from_db: Sequences have SCOP superfamily numbers only \n";
		  }elsif($Keys[$e]=~/^ *(\S+)[\,]*$/){
		      print "\n# (4) fetch_sequence_from_db: Sequences DON't have ranges or SCOP superfam numb.\n";
							$NAME=$1;
		  }

					print "\n# BEFORE reading in DB file. \$NAME is  $NAME";

		  if($seq_with_index{$NAME}=~/(\d+)/        # It is importnt having $seq_with_index{$Keys[$e]}
			   or $seq_with_index{$Keys[$e]}=~/(\d+)/
			   or $seq_with_index{"$NAME\,"}=~/(\d+)/    # this is for overcoming '>xxxx,'  entry(the comma)
                           or $seq_with_index{"$NAME\_"}=~/(\d+)/    # to handle Tim's  >c1eru_ 3.30.1.1.4
                           or $seq_with_index{"pdb\_$NAME"}=~/(\d+)/ # to handle Sarah's pdb_xxxxx
                           or $seq_with_index{"nr\_$NAME"}=~/(\d+)/ # to handle Sarah's nr_xxxxx
			   ){
			   my $finding_position= $1-300;
			   if( $finding_position >= 0 ){   seek(DB_FASTA, $1-300, 0);  # -300 is necessary
			   }elsif($finding_position < 0){  seek(DB_FASTA, 0, 0); }      ## This is essential !!!

			   while(<DB_FASTA>){
				  if(!$seq_found1){
					if(/\> *$NAME[\,_]? *\d*/){
						  $seq_found1=1;
						  print "\n# $NAME  is found in DB, Good ";
					  }
				  }else{
					  if(/^ *(\w+) *$/ ){
					      $sequence .=$1;  ## you should use $1 to avoid including NEW line
						  unless(eof DB_FASTA){ next   ## This is critically important to prevent error.
						  }else{ goto PUT_SEQ }     ## If the last seq has only one single line seq string, it could be a problem
					  }elsif( (/^ *\> *\S+/)  or (eof DB_FASTA) ){
                                               #======= When range is defined ==================
                                               PUT_SEQ:
                                               if($R_start =~/\d+/){
                                                                                                                     $sequence{$Keys[$e]}=substr($sequence, $R_start, $R_leng); next; #
                                                                                                                     print "\n# $sequence{$Keys[$e]}  is put to \%sequence";
                                               }
                                               #======= To handle superfamily information given ==========
                                               if($super_fam_class){
                                                       $sequence{$Keys[$e]}=$sequence;
                                                       $acquired_seq_count++;
                                               }
                                               #======= When range is NOT defined ==================
                                               else{
                                                       $sequence{$Keys[$e]}=$sequence;
                                                                                                                     print "\n# $sequence  is put to \%sequence";
                                               }
                                               $R_start='';
                                               $sequence='';
                                               $seq_found1='';  ## reset $R_start, $seq_found1,,
                                               next F0;
					  }
				  }
			   }

		  }else{
			   print "\n# Error, the sequence pos for $NAME (from $Keys[$e]) in DB doesnt exist in xxxx.idx file?\n";
		  }
	   }
	   close DB_FASTA;
	}
	#print "\n# (6) fetch_sequence_from_db: counted fetched seqs: $found_seq_count, $acquired_seq_count";
	#print "\n# (7) fetch_sequence_from_db: Fetching seq has finished \n";

	return(\%sequence);
}

#_________________________________________________________________________________
# Title     : open_sso_files
# Usage     :  @sso=@{&open_sso_files(@file, $add_range, $add_range2, "u=$upper_expect_limit",
#			                            "l=$lower_expect_limit", "m=$margin", $new_format)};
# Function  : This reads the parseable( -m 10 option)
#              and non-parseable form of ssearch program output
#             If you give 5 files, it produces 5 hashes as a ref of array.
#             This understands xxxx.gz files.
#             This reads FASTA -m 10 output, too.
# Example   :
#  717    0         0.343  16    373    EC1260_16-373              74    434    YBL6_YEAST_74-434
#  348    9e-16     0.500  113   233    EC1260_113-233             27    146    YDBG_ECOLI_27-146
#  472    2.9e-08   0.271  13    407    EC1260_13-407              148   567    YHJ9_YEAST_148-567
#  459    1.9e-22   0.260  1     407    EC1260_1-407               65    477    YLQ6_CAEEL_65-477
#  452    4.5e-14   0.275  1     407    EC1260_1-407               103   537    YSCPUT2_103-537
#  1131   0         0.433  1     407    EC1260_1-407               112   519    ZMU43082_112-519
#
#  Input SSO file example)-> below
#
#   >>MG032 ATP-dependent nuclease (addA) {Bacillus subtilis  (666 aa)
#    Z-score: 88.3 expect()  1.9
#   Smith-Waterman score: 77;  27.143% identity in 70 aa overlap
#
#           30        40        50        60        70        80
#   MJ0497 RSAGSKGVDLIAGRKGEVLIFECKTSSKTKFYINKEDIEKLISFSEIFGGKPYLAIKFNG
#                                        : .. ...  . .:.:::. :: : ..:
#   MG032  HDKVRYAFEVKFNIALVLSINKSNVDFDFDFILKTDNFSDIENFNEIFNRKPALQFRFYT
#        200       210       220       230       240       250
#
#           90       100             110       120       130
#   MJ0497 EMLFINPFLLSTNGK------NYVIDERIKAIAIDFYEVIGRGKQLKIDDLI
#          .   ::   :: ::.      : ....... . ::. . :
#   MG032  K---INVHKLSFNGSDSTYIANILLQDQFNLLEIDLNKSIYALDLENAKERFDKEFVQPL
#        260          270       280       290       300       310
#
# Parseable form -m 10 option =========================================
#   >>>MJ0497.fa, 133 aa vs GMG.fa library
#   ; pg_name: Smith-Waterman (PGopt)
#   ; pg_ver: 3.0 June, 1996
#   ; pg_matrix: BL50
#   ; pg_gap-pen: -12 -2
#   >>MG032 ATP-dependent nuclease (addA) {Bacillus subtilis
#   ; sw_score:  77
#   ; sw_z-score: 88.3
#   ; sw_expect    1.9
#   ; sw_ident: 0.271
#   ; sw_overlap: 70
#   >MJ0497 ..
#   ; sq_len: 133
#   ; sq_type: p
#   ; al_start: 58
#   ; al_stop: 121
#   ; al_display_start: 28
#
# Warning   : By default, the SW score comes to the first
#             If expect value is not found, it becomes '0'
#             By default, the offset of seq match with a seq name like seq_30-40
#               will be 30 not 1.
#             It ignores special chars like , : .prot in the name (eg, AADF_FASDF: will be AADF_FASDF)
# Keywords  : open_ssearch_output_files, ssearch_output, ssearch, FASTA,
# Options   : _  for debugging.
#             #  for debugging.
#             u= for upper E value limit
#             l= for lower E value limit
#             r  for attaching ranges to out seq names (eg> HI0001_1-20 as a key)
#             U  for making the matched seqname to upppercase
#             L  for making the matched seqname to lowercase
#             R  for attaching ranges to out seq names for both TARGET and MATCH
#             n  for new format (msp2)
#             a  for getting alignments of the pair
#
# Category  :
# Version   : 4.5
#----------------------------------------------------------------------------
sub open_sso_files{
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
    my (@out_refs, @SSO, $create_sso, $parseable, @OUT, @temp_sso_lines,
            %match, $attach_range_in_names, $margin, $uppercase_seq_name,
            $lowercase_seq_name, $target_seq, $new_format, $get_alignment,
            $pvm_version_fasta_out, $original_target_seq, $big_msp_out_file);

    my ($upper_expect_limit, $lower_expect_limit)=(50,0);

    if($char_opt=~/R/){  $attach_range_in_names2=1; };
    if($char_opt=~/r2/){ $attach_range_in_names =1; $attach_range_in_names2=1 };
    if($char_opt=~/r/){  $attach_range_in_names =1; };
    if($char_opt=~/c/){  $create_sso   ='c' };
    if($char_opt=~/n/){  $new_format   ='n' };
    if($char_opt=~/a/){  $get_alignment='a' };
    if($char_opt=~/U[pperPPER]*/){ $uppercase_seq_name='U' };
    if($char_opt=~/L[owerOWER]*/){ $lowercase_seq_name='L' };
    if($vars{'u'}=~/(\.?\d+)/){ $upper_expect_limit = $vars{'u'} };
    if($vars{'l'}=~/(\.?\d+)/){ $lower_expect_limit = $vars{'l'} };
    if($vars{'m'}=~/\d+/){ $margin = $vars{'m'} };
    $attach_range_in_names2=$attach_range_in_names=1;

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # opening file input (can handle .gz  files)
    #_______________________________________________
    if(@file < 1 and @array > 0){
         for($i=0; $i< @array; $i++){
              @sso=@{$array[$i]};
         }
         print "\n# (INFO) \@sso has ", scalar(@sso), " lines. \n"  if $verbose;
         if(@sso > 3000){ # if @sso is very big, I remove the useless contents
             print "\n# (INFO) open_sso_files: size of \@sso for $file[$i] exceeds 3000 lines, ", scalar(@sso), " !!! \n";
         }
         push(@OUT, &read_sso_lines(\@sso, $create_sso,
                    "u=$upper_expect_limit",
                    "l=$lower_expect_limit",
                    $attach_range_in_names,
                    $attach_range_in_names2,
                    $new_format, $get_alignment) );
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Opening input FILE!
    #_______________________________________________
    else{
         print "\n# open_sso_files : processing @file \n\n";
         for($i=0; $i< @file; $i++){
              if($file[$i]=~/\S+\.msp *$/){ $big_msp_out_file=$file[$i]; splice (@file, $i, 1); $i--;
              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              # Opening zipped file
              #_______________________________________________
              }elsif($file[$i]=~/\S+\.\gz$/ or -B $file[$i]){  ## if file has xxxx.gz extension
                  my (@sso);
                  @sso=`gunzip -c $file[$i]`;
                  if(@sso < 30){  @sso=`zcat $file[$i]`; }      # if zcat fails to produce output use gunzip -c
                  if(@sso > 3000){ # if @sso is very big, I remove the useless contents
                      print "\n# open_sso_files: size of \@sso for $file[$i] exceeds 3000 lines, ", scalar(@sso), " !!! \n";
                  }
                  push(@OUT, &read_sso_lines(\@sso, $create_sso,
                      "u=$upper_expect_limit",
                      "l=$lower_expect_limit",
                      $attach_range_in_names,
                      $attach_range_in_names2,
                      $new_format, $get_alignment) );
              }
              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              # Opening plain file(not zipped)
              #_______________________________________________
              elsif($file[$i]=~/\S+\.[fsm]?sso/ or $file[$i]=~/\S+\.out/ or $file[$i]=~/\S+\.fso/){
                  print "\n# (INFO) openning text File format xxxx.[fms]sso $file[$i]";
                  open(SSO, "$file[$i]") or die "\n# (ERROR) open_sso_files: Failed to open $file[$i]\n";
                  my @sso=<SSO>;
                  if(@sso < 30){  @sso=`zcat $file[$i]`; }      # if zcat fails to produce output use gunzip -c
                  if(@sso > 3000){ # if @sso is very big, I remove the useless contents
                      print "\n# (INFO) open_sso_files: size of \@sso is for $file[$i] exceeds 3000 lines, ",
                             scalar(@sso), " !!! \n";
                  }
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  # Calling read_sso_lines sub
                  #_____________________________________
                  push(@OUT, &read_sso_lines([@sso], $create_sso,
                       "u=$upper_expect_limit",
                       "l=$lower_expect_limit",
                       $attach_range_in_names,
                       $attach_range_in_names2,
                       $new_format, $get_alignment) );
                  close SSO;
              }
         }
    }
    print "\n# \@OUT has ", scalar(@OUT), " elements \n" if $verbose;
    return(\@OUT); # @OUT has refs of hashes  (\%xxx, \%YYY, \%XXX,,,,)
}



#________________________________________________________________________
# Title     : write_fasta
# Usage     : many argments:  $seq_hash_reference  and $output_file_name
#             takes a hash which has got names keys and sequences values.
# Function  : writes multiple seqs. in fasta format (takes one or more seq.!!)
#             This needs hash which have 'name' 'actual sequence as value'
#
#             To print out each fasta seq into each single file, use write_fasta_seq_by_seq
#             This can rename seq names
#
# Example   : &write_fasta(\%in1, \$out_file_name, \%in2, \%in3,..., );
#             << The order of the hash and scalar ref. doesn't matter. >>
# Warning   : The default output file name is 'default_out.fa' if you do not
#             specify output file name.
#             OUTput file should have xxxxx.fa or xxxx.any_ext NOT just 'xxxxx'
# Keywords  : write_fasta_file, print_fasta_file, write fasta file, fasta_write
#             show_fasta, write_sequence_fasta, write_fasta_files,
# Options   : v for STD out.
#             r for rename the sequences so that Clustalw would not complain with 10 char limit
#               so result wuld be:  0 ->ASDFASDF, 1->ASDFASFASF, 2->ADSFASDFA
# Returns   :
# Argument  :
#   $sort_seq_names=s by s  ## in writing sorted sequences are written
#
# Version   : 2.6
#--------------------------------------------------------------------
sub write_fasta{
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
  #$| = 1;

  my($sort_seq_names, $string, $string_leng, $na,$out_file_name_provided);
  my($output_file) ='default_out.fa'; ### when no output file name is given, this is used
  if(@file>0){
	$output_file = $file[0];
	$out_file_name_provided=1;
  }else{ $output_file='default_out.fa'; }

  if($char_opt=~/s/){
      $sort_seq_names='s';
  }

  for ($n=0 ; $n < @hash; $n ++){
     my(@keys);
	 my %hash=%{$hash[$n]};
     if($sort_seq_names){
         @keys=sort keys %hash;
     }else{
         @keys= keys %hash;
     }
	 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 # When only one seq is given, use the seq name as output file
	 #________________________________________________________________
	 if(@hash==1 and @keys==1 and @file < 1){
	    $output_file="$keys[0]\.fa";
	 }elsif(@file < 1){
	    $output_file="default_fa_$n\.fa";
	 }

	 open (FASTAS_WRITE,">$output_file");      # $string is the seq string.

	 for ($i=0; $i < @keys; $i++){
		$na= $keys[$i];
		$string = "\U$hash{$na}";
		$string=~s/[\n \.-]//g;	    # replaces all non-chars to null. '_' is used for stop codon
		if($char_opt=~/r/){  # rename the seqeunces with '0, 1, 2, 3," etc for  clustalw
		   $na=$i;
		}

		if($debug == 1){
			print ">$na\n";
			print FASTAS_WRITE ">$na\n";
	    }elsif($char_opt=~/v/){
		    print ">$na\n";
		    print FASTAS_WRITE ">$na\n";
		}else{
		    print  FASTAS_WRITE ">$na\n";
		}

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#  Main algorithm of writing in 60 char leng line
		#_____________________________________________________
		$string_leng=length($string);
		for($j=0; $j< $string_leng; $j+=60){
			if($debug == 1){
				printf "%.60s\n", substr($string,$j,60);
				printf FASTAS_WRITE "%.60s\n", substr($string,$j,60);
			}elsif($char_opt=~/v/i){
				printf "%.60s\n", substr($string,$j,60);
				printf FASTAS_WRITE "%.60s\n", substr($string,$j,60);
			}else{
			   printf FASTAS_WRITE "%.60s\n", substr($string,$j,60);
			}
		}
	 }
	 close FASTAS_WRITE;
  }
  if( $out_file_name_provided != 1){
	  print "\n\n# You didnt give out file name, $output_file  used\n";
  }
  if( -s $output_file ){
	 if($verbose=~/\S/){ ## if v option is given, mesg is omitted to prevent comments to a redirected output
	    print "\n# Sequences were written in  $output_file ";
	 }
  }else{
	 print "\n# The size of written outfile \"$output_file\" is 0, error \n\n"
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
	if($vars{'f'}=~/\S+/){     $factor= $vars{'f'}; }
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
# Version   : 1.8
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
   my @seq_names_in_clu=@{$array[0]};
   $link_or_not=0;
   my $factor=7; # 7 means 70% now

   if($vars{'f'}=~/(\S+)$/){ $factor=$1 }
   if($char_opt=~/o/){ $optimize=1 }
   if($char_opt=~/S/){ $short_region='S'; }
   if($char_opt=~/L/){ $large_region='L';   }
   if($char_opt=~/A/){ $average_region='A'; }
   if($char_opt=~/v/){ $verbose=1 }

   if($verbose){ print "\n# (1) cluster_merged_seqlet_sets: Checking linkage and merging <<<<<>>>>>\n@seq_names_in_clu\n";   }

   F1: for($i=0; $i< @seq_names_in_clu; $i++){
      @splited1=split(/ +/, $seq_names_in_clu[$i]);
      for($j=0; $j< @seq_names_in_clu; $j++){
        if($seq_names_in_clu[$i] eq $seq_names_in_clu[$j]){ next  }
        @splited2=split(/ +/, $seq_names_in_clu[$j]);

        $link_or_not=${&check_linkage_of_2_similar_seqlet_sets(\@splited1, \@splited2, "f=$factor")};
        print "\n +++++ \$link_or_not is  $link_or_not +++" if $verbose;
        if($link_or_not==1){
             if($verbose){
                 print "\n# (2) cluster_merged_seqlet_sets: \n $seq_names_in_clu[$i]  \n and $seq_names_in_clu[$j] \n are linked \n";
             }

             if($optimize){ ##---- This will also remove similar seqlets, not only identical ones
                $seq_names_in_clu[$i]=join(' ', sort @{&remove_similar_seqlets( [@splited1, @splited2],
                                            $short_region, $large_region, $average_region)} );
             }else{
                $seq_names_in_clu[$i]=join(' ', sort @{&remove_dup_in_array( [@splited1, @splited2])} );
             }
             splice(@seq_names_in_clu, $j,1);
             $j--; $i--;
             next F1;
		 }
	  }
   }
   return(\@seq_names_in_clu);
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
# Version   : 2.4
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
   my $length_thresh =30; # sequence length threshold. overlap less than this will be ignored

   if($char_opt=~/v/){     $verbose = 'v'
   }if($char_opt=~/z/){    $optimize = 'z'
   }if($char_opt=~/S/){    $short_region='S';
   }if($char_opt=~/L/){	   $large_region='L';
   }if($char_opt=~/A/){	   $average_region='A'; }

   if($vars{'T'}=~/\d+/){
	  $length_thresh=$vars{'T'}; print "\n# merge_sequence_in_msp_chunk: Thresh is $length_thresh\n" if (defined $verbose);
   }if($vars{'f'}=~/\S+/){
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
	  if( ($R_leng1 < $length_thresh) || ($ssearch_score1 < $score) ){ splice(@seqlets, $i, 1); $i--; next; }
	  if( $evalue_found1 > $evalue){ splice(@seqlets, $i, 1); $i--; next; }

	  F2: for($j=0; $j < @seqlets; $j ++){
		 if($seqlets[$i] eq $seqlets[$j]){ next };
		 if($seqlets[$j]=~/^ *((\d+) +(\d+\.?[e\-\d]*) +(\d+) +(\d+) +(\S+) +(\d+) +(\d+)) +(\S+) *(.*)/){
			($long_match2, $enq_seq2, $mat_seq2, $R_start2, $R_end2)=($1, $6, $9, $4, $5);
			$R_leng2=$R_end2-$R_start2;			$ssearch_score2=$2;			$evalue_found2= $3;
	     }
		 if( ($R_leng2 < $length_thresh)||($ssearch_score2 < $score) ){ splice(@seqlets, $j, 1); $j--; next; }
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
#______________________________________________________________________________
# Title     : create_sorted_cluster
# Usage     : &create_sorted_cluster
# Function  : to make a "sorted_cluster_file" from the .clu files in a directory
# Example   :
# Keywords  : make_cluster_file, sort_clu_files
# Options   :
# Author    : Sarah A. Teichmann, modified by Jong
# Date      : 19th September 1997
# Version   : 1.7
#--------------------------------------------------------------------------------
sub create_sorted_cluster{
    my ($i, $q, $p, $n, $j, $clufile, @filecontent, $new_gene, $cluster_size,
        @clufiles, @cluster_sizes_new_unsorted, @cluster_sizes_new,
        $newclus_number, %hash);

    @clufiles=@{$_[0]};
    if(@clufiles < 1 ){
        @clufiles=@{&read_file_names_only('.','.clu')};
        print "\n# $0, create_sorted_cluster: \@_ is empty, reading PWD to get xxx.clu files\n";
        if(@clufiles < 1){
           print "\n# $0, create_sorted_cluster: I couldn\'t find any clu files, dying\n";
           exit;
        }
    }
    if($_[1]){  $good_cluster_file=${$_[1]}; }

    print "\n# (1) create_sorted_cluster: found \@clufiles are
                @clufiles\n";

    for ($i=0; $i < @clufiles; $i++) {
         open (CLU_FILE, "<$clufiles[$i]") or
            die "\n# $0: create_sorted_cluster: error opening $clufiles[$i]";

         my $cluster_size;
         while(<CLU_FILE>){
             if( /^ *Cluster +size +(\d+)/i){
                 $cluster_size=$1;
                 $hash{$cluster_size} .=$_;
             }elsif (/^ *Cluster +[number]* *\d+/i) {
                 $hash{$cluster_size} .=$_;
             }elsif (/^ *\d+ +\d+ +\S+/) {
                 $hash{$cluster_size} .=$_;
             }
        }
   }

   @sorted_by_size=sort { $a<=>$b } keys %hash;
   $good_cluster_file2="sorted_cluster_file\.gclu";

   open(GOODCLUS, ">$good_cluster_file") or die "\n# $0 create_sorted_cluster: I can not open $good_cluster_file\n";
   open(GOODCLUS2, ">$good_cluster_file2") or die "\n# $0 create_sorted_cluster: I can not open $good_cluster_file\n";
   for($i=0; $i< @sorted_by_size; $i++){
       print GOODCLUS $hash{$sorted_by_size[$i]};
       print GOODCLUS2 $hash{$sorted_by_size[$i]};
   }
   close(GOODCLUS);
   close(GOODCLUS2);
   return(\$good_cluster_file);
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

#_________________________________________________________________________________
# Title     : read_sso_lines
# Usage     : &read_sso_lines([@sso], $create_sso, $attach_range_in_names,
#                 $attach_range_in_names2, $new_format, $get_alignment) );
# Function  : Main subroutine for open_sso_files. This calls either machine
#              readable or unreadable form parsing subroutine
# Example   :
# Keywords  : read_sso_lines_in_array
# Options   : a c r r2 n
#             u= for upper E value limit
#             l= for lower E value limit
# Category  :
# Version   : 1.4
#----------------------------------------------------------------------------
sub read_sso_lines{
    my (@out_refs, $parseable, @SSO, $create_sso, $i, $j, $k, $attach_range_in_names);
    my ($upper_expect_limit, $lower_expect_limit)=(50,0); ##<<--- DEFAULT

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```
    # Processing input and options, don't worry!
    #__________________________________________
    for($i=0; $i< @_; $i++){
       if($_[$i]=~/u=(\S+)/){    $upper_expect_limit=$1 }
       elsif(ref($_[$i]) eq 'ARRAY'){ @SSO=@{$_[$i]};   }
       elsif($_[$i]=~/l=(\S+)/){ $lower_expect_limit=$1 }
       elsif($_[$i]=~/^c$/){     $create_sso = 'c' }
       elsif($_[$i]=~/^a$/){     $get_alignment='a'; }
       elsif($_[$i]=~/^r$/){   $attach_range_in_names='r' }
       elsif($_[$i]=~/^r2$/){   $attach_range_in_names2='r2' }
       elsif($_[$i]=~/^n$/){   $new_format='n' }
    }
    print "\n# (INFO) \$attach_range_in_names2 is $attach_range_in_names2\n" if $attach_range_in_names2;

    #~~~~~~ Checking if sso is a parseable form or not~~~~~~~~~~~~~
    TEMP:for($k=0; $k < @SSO; $k++){
          if($SSO[$k] =~ /\>\>\>/  or $SSO[$k] =~ /^ *\; \S+\:/ ){
                  $parseable++;  if($parseable >= 10){  last TEMP;     }
          }elsif($SSO[$k]=~/^  +\:+/){ $parseable--;
          }elsif($SSO[$k] =~ /^ +1\>\>\>(\S+)/){ $pvm_version_fasta_out=1; $parseable +=10; $original_target_seq=$1; last TEMP;
          }
    }
    if($parseable >= 10){
        @out_refs=@{&read_machine_readable_sso_lines(\@SSO,
                "u=$upper_expect_limit",
                "l=$lower_expect_limit",
                $get_alignment, $create_sso,
                $new_format, $attach_range_in_names,
                $attach_range_in_names2)};
    }else{
        @out_refs=@{&read_machine_unreadable_sso_lines(\@SSO,
                $get_alignment, $create_sso,
                "u=$upper_expect_limit",
                "l=$lower_expect_limit",
                $new_format, $attach_range_in_names,
                $attach_range_in_names2)};
    }
    return(@out_refs);
}



#______________________________________________________________________________
# Title     : sso_to_msp
# Usage     : &sso_to_msp(@ARGV, $single_out_opt);
# Function  : This takes sso file(s) and produces MSP file. It
#             concatenate sso file contents when more than one
#             sso file is given.
# Example   : &sso_to_msp(@ARGV, 'OUT.msp', $single_out_opt);
# Warning   : This capitalize all the input file names when
#              producing xxxxx.msp. xxxxx.sso -> XXXX.sso
# Keywords  : sso_file_to_msp_file, convert_sso_to_msp,
# Options   : _  for debugging.
#             #  for debugging.
#             v  for showing the MSP result to screen
#             s  for making single MSP file for each sso file
#                    as well as big MSP file which has all sso
#             u= for upper expectation value limit
#             l= for lower expect val limit
#             s= for single file name input eg. "s=xxxxx.msp"
#             n  for new format (msp2 format)
#             r  for adding range
#             r2 for adding ranges in all sequence names
#
# Returns   : the file names created (xxxx.msp, yyyy.msp,,,,)
# Argument  :
# Version   : 2.6
#-----------------------------------------------------------------------------
sub sso_to_msp{
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
   my ($upper_expect_limit, $lower_expect_limit)=(50, 0);
   my (%sso, @sso, @SSO, $big_out_msp1,  @final_out, $big_out_msp2,
	   $create_sso, $single_out_opt, $add_range, $add_range2, $big_out_msp,
	   $Evalue_thresh, $new_format, $Score_thresh, $margin, $single_file_name);
	if($vars{'u'}=~/([\.\d]+)/){ $upper_expect_limit = $vars{'u'} };
	if($vars{'l'}=~/([\.\d]+)/){ $lower_expect_limit = $vars{'l'} };
	if($vars{'t'}=~/(\d+)/){ $Score_thresh  = $vars{'t'} };
	if($vars{'m'}=~/(\d+)/){ $margin  = $vars{'m'} };
	if($vars{'s'}=~/\S/){ $single_file_name  = $vars{'s'} };
	if($char_opt=~/r2/){  $add_range='r'; $add_range2='r2' }
	if($char_opt=~/r/){   $add_range = 'r' }
	if($char_opt=~/c/){   $create_sso = 'c' }
	if($char_opt=~/s/){   $single_out_opt='s' }
	if($char_opt=~/n/){   $new_format='n' }
   print "\n# File given to sso_to_msp is \"@file\", Normally xxx.sso file names\n";

   if($single_file_name=~/\S/){
	   $big_out_msp=$single_file_name;
   }else{
	   for($i=0; $i < @file; $i++){
		   if($file[$i]=~/\.msp$/){ ## when output file name is given
			   $big_out_msp=$file[$i];
			   splice(@file, $i, 1);
			   $i--;
		   }elsif($file[$i]=~/^(\d+\-\d+)([_\d]*)\.[mfs]?sso/){  ## creates xxxx.msp file name from xxxx.sso
			   $big_out_msp1="\U$1"."$2"."\.msp";
			   $big_out_msp2="\U$1".".msp";
		   }elsif($file[$i]=~/^(\S+)\.[mfs]?sso$/){
			   $big_out_msp1="\U$1"."\.msp";
			   $big_out_msp2="\U$1"."_all".".msp";
			   print "\n# sso_to_msp: File matched  xxxx.sso  format \n";
		   }elsif($file[$i]=~/^(\S+)\.out$/){
			   $big_out_msp1="\U$1"."\.msp";
			   $big_out_msp2="\U$1"."_all".".msp";
			   print "\n# sso_to_msp: File matched  xxxx.out  format \n";
		   }elsif($file[$i]=~/^(\S+)\.p[rot\,]*\.ts\.gz/){
			   $big_out_msp1="\U$1".".msp";
			   $big_out_msp2="\U$1"."_all".".msp";
		   }elsif($file[$i]=~/^(\S+)\.ts\.gz/){
			   $big_out_msp1="\U$1".".msp";
			   $big_out_msp2="\U$1"."_all".".msp";
		   }elsif($file[$i]=~/^(\S+)\.out\.gz/ or $file[$i]=~/^(\S+)\.[mfs]?sso\.gz/){
			   $big_out_msp1="\U$1".".msp";
			   $big_out_msp2="\U$1"."_all".".msp";
		   }
	   }
   }
   if(defined($big_out_msp)){
	   $big_out_msp1=$big_out_msp2=$big_out_msp;
	   print "\n# \$big_out_msp is defined as \'$big_out_msp\'\n";
   }else{
	   print "\n# sso_to_msp: You did not define the big MSP file out format, so $big_out_msp1 \n";
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  (1) When File was given to this sub routine
   #__________________________________________
   if(@file == 1){   ## ONE single file input??
	  print "# one file @file is given, OUT will be: $big_out_msp1 \n";
	  @sso=@{&open_sso_files(@file, $add_range, $add_range2,
	          "u=$upper_expect_limit",
			  "l=$lower_expect_limit",
			  "m=$margin",
			  $new_format,
			  "s=$big_out_msp")};
	  push(@final_out, &write_msp_files(@sso, $big_out_msp1,
	        $single_out_opt, $add_range) );

   }elsif(@file > 1){ ## MOre than 1 file input??
	  @sso=@{&open_sso_files(@file, $add_range, $add_range2,
	        "l=$lower_expect_limit",
	        "u=$upper_expect_limit",
	        "m=$margin",
	        $new_format)};
	  push(@final_out, @{&write_msp_files(@sso, $big_out_msp2,
			$single_out_opt, $add_range)} ); ## concatenates all the hash ref to one
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  (2) When NO File but ARRAY is given
   #      Here, you can have SSO files created
   #__________________________________________
   elsif(@array >=1){
	  print "\n# In sso_to_msp, \@array is given rather than \@file";
	  @sso=@{&open_sso_files(@array, "u=$upper_expect_limit", $add_range2,
			  "l=$lower_expect_limit", $add_range, $create_sso,
			  "m=$margin", $new_format)};
	  push(@final_out, @{&write_msp_files(@sso, $big_out_msp,
						  $single_out_opt, $add_range)} );
   }
   return(\@final_out);
}
#________________________________________________________________________________
# Title     : make_clustering_summary
# Usage     : &make_summ($sorted_cluster_file)
# Function  : to make a summary file of a sorted cluster file
# Example   :
# Keywords  : summary, make_cluster_summary, subclustering summary
# Options   :
# Author    : Sarah A. Teichmann
# Date      : 19th September 1997
# Version   : 1.5
#--------------------------------------------------------------------------------
sub make_clustering_summary{
    my ($good_cluster_file, $summary_file, @filecontent, $i, $filecontent,
        $cluster_size, @cluster_sizes, $cluster_number, $number_of_clusters,
         $summary_file, @filecontent, %hash, @keys, @temp_clu, $make_separate_summary);
    $good_cluster_file=${$_[0]} || $_[0];

    $make_separate_summary=0;
    for($i=0; $i< @_; $i++){
        if($_[$i]=~/u *$/){  $make_separate_summary='u';
                             splice(@_, $i, 1)
        }
    }

    $summary_file="$good_cluster_file".".summary";

    open(CLU_FILE, "$good_cluster_file");
    while(<CLU_FILE>){
          push(@temp_clu, $_);  ## copying the content to ;
          if( /^ *Cluster +size +(\d+)/i){
              $cluster_size=$1;
          }elsif (/^ *Cluster +[number]* *(\d+)/) {
              $hash{$cluster_size} ++;
          }
    }
    close(CLU_FILE);

    open(CLU_FILE, ">$good_cluster_file"); # now overwrting it.

    if($make_separate_summary){
        open (SUMM, ">$summary_file");
        print SUMM "Cluster size    No. of clusters\n";
    }

    print CLU_FILE "Cluster size    No. of clusters\n";
    @keys=sort {$a<=>$b} keys %hash;
    for ($i=0; $i<@keys; $i++){
        print SUMM "     $keys[$i]               $hash{$keys[$i]}\n" if $make_separate_summary;
        print CLU_FILE "     $keys[$i]               $hash{$keys[$i]}\n";
    }
    close (SUMM) if $make_separate_summary;
    print CLU_FILE "\n# This file is created by $0 with make_clustering_summary sub, Details below\n\n";
    for(@temp_clu){  print CLU_FILE $_ }
    close (CLU_FILE);
    return(\$summary_file);
}


#________________________________________________________________________________
# Title     : read_machine_unreadable_sso_lines
# Usage     : @out_refs=@{&read_machine_unreadable_sso_lines(\@SSO, $get_alignment,
#                           $create_sso, $upper_expect_limit,$new_format, $lower_expect_limit,
#                           $attach_range_in_names, $attach_range_in_names2)};
# Function  :
# Example   : output will look=>
#  ZFH1_DROME 60 d1ad3a_ 446 d1ad3a_ 0.9 71 34 3-37 253-287
#  ZFH1_DROME 60 d1ahdp_ 68 d1ahdp_ 0.00018 100 56 2-58 3-59
#  ZFH1_DROME 60 d1crka2 282 d1crka2 8.4 58 50 5-55 73-123
#  ZFH1_DROME 60 d1dkza_ 215 d1dkza_ 4.9 59 40 1-41 112-152
#  ZFH1_DROME 60 d1ecra_ 305 d1ecra_ 3.2 63 47 9-56 201-248
#
# Keywords  : read_normal_sso_lines
# Options   : a c r r2 n
#             u= for upper E value limit
#             l= for lower E value limit
# Version   : 1.3
#--------------------------------------------------------------------------------
sub read_machine_unreadable_sso_lines{
   my ($upper_expect_limit, $lower_expect_limit)=(50,0);
   my (@SSO, @out_refs, $match_seq, $match_evalue, $alignment_found,
       $match_zscore, $query_seq_len, $space, %match, $entry_found,
       $match_gap_len, $query_seq_name, $i, $j, $sw_score, $match_evalue,
       $sw_ident, $query_seg_start, $query_seg_end, $match_seg_start,
       $match_seg_end, $query_gap_len, $match_gap_len, $overlap,
       $gap_start, $start, $match_range);
   for($i=0; $i< @_; $i++){
       if($_[$i]=~/u=(\S+)/){    $upper_expect_limit=$1 }
       elsif(ref($_[$i]) eq 'ARRAY'){ @SSO=@{$_[$i]};   }
       elsif($_[$i]=~/l=(\S+)/){ $lower_expect_limit=$1 }
       elsif($_[$i]=~/^c$/){     $create_sso = 'c' }
       elsif($_[$i]=~/^a$/){     $get_alignment='a'; }
       elsif($_[$i]=~/^r$/){   $attach_range_in_names='r' }
       elsif($_[$i]=~/^r2$/){   $attach_range_in_names2='r2' }
       elsif($_[$i]=~/^n$/){   $new_format='n' }
   }

   print "\n# (i)read_machine_unreadable_sso_lines : You have put non-parseable format of xxxx.sso\n";
   print "#      Did you set \'M\' option in do_sequence_search? \n";

   for($j=4; $j< @SSO; $j++){
       if($SSO[$j]=~/^ *\S+\: +(\d+) +\w+ *$/){ ## matching 'hbox.fa: 60 aa'
            $query_seq_len=$1;
            print "\n (INFO) target seq len is  $query_seq_len \n";

       }elsif($query_seq_len and $SSO[$j]=~/^ \>(\w[\w\-\.\/\\]+)/){ ##  matching  >MJ0497
            $query_seq_name=$1;
            $j+=4; ## jumping to skip the stat bars
            print "\n# (i)read_machine_unreadable_sso_lines : Found Query seq=> $query_seq_name ";
            next;
       }

       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # Matching the Match sequences
       #_________________________________________________
       if($query_seq_name and $SSO[$j]=~/^ {0,4}\>\> *(\S+) +.+\((\d+) aa\) *$/){ # # matching >>MG032 ATP-d (addA) Bacillus subtilis  (666 aa)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
            # This is the assignment to OUT hash
            #________________________________________________
            if($alignment_found > 0){ ## If I have already read >>xxxx
                $match{$match_seq}=sprintf("%-5s %-8s %-6s %-4s %-5s %-30s %-4s %-5s %s\n",
                   $sw_score, $match_evalue, $sw_ident, $query_seg_start, $query_seg_end, $query_seq_name,
                   $match_seg_start, $match_seg_end, $match_seq);
            }
            $match_seq=$1;
            $entry_found=1;    $alignment_found=0;
            $query_gap_len=$match_gap_len=$match_seg_start=$match_seg_end=0;
            $query_seg_start=$query_seg_end=$overlap=0;
            undef($space);
            print "\n# (INFO) read_machine_unreadable_sso_lines : Found MATCHed seq $match_seq\n" if $verbose;
       }elsif($entry_found and $SSO[$j]=~/ *Z\-score: *(\S+) *expect\( *\) +(\S+)/){ ## getting Evalue
            $match_zscore=$1;   # just in case needed
            $match_evalue=$2;
            if($upper_expect_limit < $match_evalue or $lower_expect_limit > $match_evalue){
                $match_seq=$alignment_found=$entry_found=$query_gap_len='';
            }
       }elsif($entry_found and $SSO[$j]=~/Smith\-Waterman +score\: +(\d+)\; *(\S+) *\% *identity *in *(\d+) *aa overlap/i){
            $sw_score=$1;
            $sw_ident=$2;
            $overlap=$3;
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # Measuring the gap before position info "^         50        60$"
       #_______________________________________________________________________
       }elsif( $entry_found and !$alignment_found  and $SSO[$j]=~/^( +)(\d+) +\d+/  ){
            $start=$2;                                ## $start is the digit
            $gap_start=length($1)+length($start)-1;
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # $alignment_found is set HERE!
            #_____________________________________
            $alignment_found=1;
            print "\n# $query_seq_name $match_seq \$gap_start is  $gap_start, \$start: $start" if $verbose;
       }elsif($alignment_found==1 and $SSO[$j]=~/^( +)[\.\:]/ ){ ### matching    .: .: : ::     :.:..: :.. .. ..
            $space=length($1);
            $query_seg_start=$space-$gap_start+$start;
            $query_seg_end=$query_seg_start+$overlap;
            $query_range="$query_seg_start-$query_seg_end";
       }elsif($space and $alignment_found==1 and  $SSO[$j]=~/^( +)(\d+)/ ){ ##  !! $alignment_found == 1 !!
            $alignment_found++; ## <<--- $alignment_found++;  should increment!
            $match_gap_start=length($1)+length($2)-1;
            $match_start=$2;
            $match_seg_start=$space-$match_gap_start+$match_start;
            $match_seg_end=$match_seg_start+$overlap;
            $match_range ="$match_seg_start-$match_seg_end";
            #$match{$match_seq}.="$query_range $match_range ";
            #print "\n $query_seq_name $match_seq $match_evalue $overlap $query_range $match_range";
       }
    }# end of for $j
    if( ($create_sso=~/c/) && (@file < 1) ){
       open (SSO2, ">$query_seq_name\.sso");
       print SSO2 @SSO, "\n";
       print "\n# (INFO) $query_seq_name\.sso  is created";
       close SSO2;
    }
    push(@out_refs, \%match);
    return(\@out_refs);
}# end of for $i



#________________________________________________________________________________
# Title     : convert_clu_to_sso_to_msp  (do not use this)
# Usage     : &clu_to_sso_to_msp(\$clu);
# Function  : reads in a big single linkage cluster file(or normal cluster file)
#              and creates a big msp file which contains all the entries in the
#              cluster file (usually with the extension of sclu or clu)
#             This normally reads in xxxx.mso, xxxx.sso like files, but if the
#              corresponding  xxx.msp file already exists, it concatenates them to
#              make a bigger one.
# Example   :
# Keywords  : clu_2_sso_2_msp, cluster_to_msp, cluster_to_sso_to_msp
#              clu_to_sso_to_msp
# Options   :
# Category  :
# Version   : 1.8
#--------------------------------------------------------------------------------
sub convert_clu_to_sso_to_msp{
     my($i, $j, $k, $s, $u, $p, $m, $n, $y, @possible_extensions, @list,
          @final_files, @U_L_case, $file, @file, @written_msp_files);

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # Opening cluster file (xx.clu)
     # %clus looks like this:  2-507     YGR041W YLR353W
     #                         3-308     YDR222W YDR346C YLR225C
     #                         2-184     YCL066W YCR040W
     #______________________________________________________________
     my $clu=${$_[0]} || $_[0];
     if($verbose){
           print "\n# clu_to_sso_to_msp : \"$clu\" is given
                and I am processing it with clu_to_sso_to_msp\n" if defined $clu;
     }
     my %clus=%{&open_clu_files(\$clu)};
     my @keys= keys %clus;
     my $num_of_cluster=@keys=@{&sort_by_cluster_size(\@keys)};

     print "# $0: clu_to_sso_to_msp: No. of cluster=$num_of_cluster after open_clu_files \n" if $verbose;

     &show_array(\@keys) if $verbose;
     &show_hash(\%clus) if $verbose;
     @possible_extensions=('msp', 'sso', 'msso', 'msso.gz','fsso', 'ssso', 'fso', 'out', 'prot.sso', 'prot.ts');
     @U_L_case=('\U', '\L');

     for($i=0; $i< @keys; $i++){
         my (@list, @final_files, $clus_name, $big_out_msp, @msp_hashes);
         $clus_name=$keys[$i];
         unless($single_file_name=~/\S/){
              $big_out_msp="$clus_name\_cluster\.msp"; #<<<----- final output name
         }else{
              $big_out_msp=$single_file_name;
         }
         push(@written_msp_files, $big_out_msp); ## This is the output of this sub

         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #  If $clus_name.msp is already there, skip
         #_____________________________________________
         if( (-s $big_out_msp) > 100  and !$over_write ){
             print "\n# clu_to_sso_to_msp : $big_out_msp MSP file already exists, skipping\n";
             print "#    Use  \$over_write option \'o\' to start all over again or \n";
             print "#    delete clustering files like XX-XX_cluster.clu to go on\n";
             next ;
         }
         $num_of_seq_member=@list=split(/ +/, $clus{$keys[$i]}); # @list has (HIU001, HI002, HI333, MJ111, etc)
         print "# $0: clu_to_sso_to_msp: No. of seq member=$num_of_seq_member after split \n" if $verbose;

         FOR0: for($j=0; $j < @list; $j++){
                  my($sub_dir_head, $file_name_low, $file_name_up, $file_name_prot_low, @sub_dir_heads,
                     $file_name_prot_up, $file_name_low_gz, $file_name_up_gz,
                     $file_name_prot_low_gz, $file_name_prot_up_gz);

                  $each_seq_name=$list[$j];
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  #  Here I take chars from the sequ names, as dirs have fragments of chars
                  #_______________________________________________________________________________
                  for($s=1; $s <=2 ; $s++){  ## here, number 2 indicates, I check single or 2 char sub dir names
                       $sub_dir_head= substr($list[$j], 0, $s);
                       push(@sub_dir_heads, "\L$sub_dir_head") if (-d "\L$sub_dir_head" );
                       push(@sub_dir_heads, "\U$sub_dir_head") if (-d "\U$sub_dir_head" );
                  }
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  #  Checking all the possible subdirectories to crop all the sso files
                  #_______________________________________________________________________________
                  FOR1: for($p=0; $p < @sub_dir_heads; $p++){
                       $subd=$sub_dir_heads[$p];
                       FOR2 : for($e=0; $e < @possible_extensions; $e++){
                            $ext=$possible_extensions[$e];
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            #  This makes all the possible lower upper case names
                            #______________________________________________________
                            for( $u=0; $u < @U_L_case; $u++){
                                  if($U_L_case[$u]=~/U/){  $each_seq_name="\U$each_seq_name";
                                  }else{                   $each_seq_name="\L$each_seq_name"; }

                                  if(-s "$each_seq_name\.$ext"){   push(@final_files, "$each_seq_name\.$ext" ) ; next FOR0 }
                                  elsif(-s "$each_seq_name\.$ext\.gz"){ push(@final_files, "$each_seq_name\.$ext\.gz" ) ; next FOR0 }
                                  else{
                                       $file_wanted="\.\/$subd\/$each_seq_name\.$ext";
                                       if(-s $file_wanted){
                                                       push( @final_files, $file_wanted); next FOR0 }
                                       elsif(-s "$file_wanted\.gz"){
                                                       push( @final_files, "$file_wanted\.gz");
                                                       next FOR0
                                       }
                                  }
                            }
                       } # FOR2
                  } # FOR1

         } # FOR0

         print "\n# @final_files \n=============> $big_out_msp  \n\n" if $verbose;

         if(@final_files < 1){
              print "\n# clu_to_sso_to_msp :LINE no.: ", __LINE__, " ERROR: \@final_files is empty. Serious error\n";
              print "\n If you have sub dir which have more than 2 chars as names, you may increase the default 2 to 3 in the above\n";
              next;
         }
         # $write_each_msp_to_disk='w';

         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #  Check if small msp files have already made in previous steps
         #________________________________________________________________
         if($final_files[0]=~/\.msp *$/){ ##  concatenate msp into big_msp
             open(BIG_MSP_FILE, ">$big_out_msp");
             for($y=0; $y< @final_files; $y++){
                 open(SINGLE_MSP, "$final_files[$y]");
                 while(<SINGLE_MSP>){
                                 print BIG_MSP_FILE $_;
                 }
             }
             close(BIG_MSP_FILE);
             close(SINGLE_MSP);
             push(@written_msp_files, $big_out_msp);

         }else{
             if($write_each_msp_to_disk){
                  print "\# $0 : going to run open_sso_files with $write_each_msp_to_disk opt\n";
                  $big_out_msp=${&open_sso_files(\@final_files, $uppercase_seq_name, $write_each_msp_to_disk,
                                                                          "u=$upper_expect_limit", $new_format, $add_range, $add_range2, $big_out_msp, $over_write)};
                  if(-s $big_out_msp > 200){  print "\n# $0: SUCCESS to create $big_out_msp :) :) :-) :-) ?\n"; }
             }else{
                  print "\n# clu_to_sso_to_msp: I am running open_sso_files. \n";
                  @msp_hashes=@{&open_sso_files(\@final_files, $uppercase_seq_name, $write_each_msp_to_disk,
                                                                          "u=$upper_expect_limit", $new_format, $add_range, $add_range2, $big_out_msp, $over_write)};

                  &write_msp_files(@msp_hashes, $big_out_msp); ## concatenates all the hash ref to one
             }
         }
     }
     return(\@written_msp_files);
}# end of


#________________________________________________________________________________
# Title     : convert_clu_to_msp
# Usage     : @written_msp_files=@{&convert_clu_to_msp(\$single_linkage_file)};
# Function  : reads in a big single linkage cluster file(or normal cluster file)
#              and creates a big msp file which contains all the entries in the
#              cluster file (usually with the extension of sclu or clu)
#             This normally reads in xxxx.mso, xxxx.sso like files, but if the
#              corresponding  xxx.msp file already exists, it concatenates them to
#              make a bigger one.
# Example   :
# Keywords  : clu_2_sso_2_msp, cluster_to_msp, cluster_to_sso_to_msp
#              clu_to_sso_to_msp
# Options   :
# Category  :
# Version   : 2.3
#--------------------------------------------------------------------------------
sub convert_clu_to_msp{
     my($i, $j, $k, $s, $u, $v, $p, $m, $n, $y, @possible_extensions, $single_file_name,
        @seq_names, @final_files, @U_L_case, $file, @file, @name_types,
        @poss_sub_dir_heads, @written_msp_files, $Lean_output, $subdir_char_size,
        $search_file_base, $found_real_subdir_name, $found_search_prog_exention_used);

     $subdir_char_size=2; # default
     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # Opening cluster file (xx.clu)
     # %clus looks like this:  2-507     YGR041W YLR353W
     #                         3-308     YDR222W YDR346C YLR225C
     #                         2-184     YCL066W YCR040W
     #______________________________________________________________
     my $clu=${$_[0]} || $_[0];
     $Lean_output=${$_[1]} || $_[1];

     if($verbose){
           print "\n# convert_clu_to_msp : \"$clu\" is given
                and I am processing it with clu_to_sso_to_msp\n" if defined $clu;
     }
     my %clus=%{&open_clu_files(\$clu)};
     my @clusters= keys %clus;
     my $num_of_cluster=@clusters=@{&sort_by_cluster_size(\@clusters)};

     print "# (i) $0: convert_clu_to_msp: No. of cluster=$num_of_cluster after open_clu_files \n" if $verbose;

     &show_array(\@clusters) if $verbose;
     &show_hash(\%clus) if $verbose;
     @possible_extensions=('msp', 'msp.gz', 'msso', 'msso.gz','fsso', 'pbla', 'pbla.gz',
                                  'ssso', 'fso', 'out', 'prot.sso', 'prot.ts');
     @U_L_case=('\U', '\L');

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
     # Making each SINGLE linkage clu to MSP file format to be ready for divclus
     #______________________________________________________________________________
     for($i=0; $i< @clusters; $i++){
         my (@seq_names, @final_files, $clus_name, $big_out_msp, @msp_hashes);
         $clus_name=$clusters[$i];
         unless($single_file_name=~/\S/){
             $big_out_msp="$clus_name\_cluster\.msp"; #<<<----- final output name
         }else{
             $big_out_msp=$single_file_name;
         }
         push(@written_msp_files, $big_out_msp); ## This is the output of this sub

         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #  If $clus_name.msp is already there, skip
         #_____________________________________________
         if( (-s $big_out_msp) > 100  and !$over_write ){
             print "\n# (i) convert_clu_to_msp : $big_out_msp MSP file already exists, skipping\n";
             print "#    Use  \$over_write option \'o\' to start all over again or \n";
             print "#    delete clustering files like XX-XX_cluster.clu to go on\n";
             next ;
         }
         $num_of_seq_member=@seq_names=split(/ +/, $clus{$clusters[$i]}); # @seq_names has (HIU001, HI002, HI333, MJ111, etc)
         print "# $0: convert_clu_to_msp: No. of seq member=$num_of_seq_member after split \n" if $verbose;

         FOR0: for($j=0; $j < @seq_names; $j++){
               my($sub_dir_head, $file_name_low, $file_name_up, $file_name_prot_low,
                  $file_name_prot_up, $file_name_low_gz, $file_name_up_gz,
                  $file_name_prot_low_gz, $file_name_prot_up_gz);
                  $each_seq_name=$seq_names[$j];
               my @poss_sub_dir_heads=('.'); ## <<<<------- This is critically important, when 'D' opt is not used!

               if($each_seq_name=~/(\S+)_\d+\-\d+$/){
                   $each_seq_name_range=$each_seq_name;
                   $each_seq_name=$1;
                   @name_types=($each_seq_name, $each_seq_name_range);
               }else{
                   @name_types=($each_seq_name);
               }
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               #  Here I take chars from the sequ names, as dirs have fragments of chars
               #_______________________________________________________________________________
               for($s=1; $s <= $subdir_char_size ; $s++){  ## here, number 2 indicates, I check single or 2 char sub dir names
                   $sub_dir_head= substr($seq_names[$j], 0, $s);
                   push(@poss_sub_dir_heads, "\L$sub_dir_head") if (-d "\L$sub_dir_head" );
                   push(@poss_sub_dir_heads, "\U$sub_dir_head") if (-d "\U$sub_dir_head" );
               }
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               #  Checking all the possible subdirectories to crop all the sso files
               #_______________________________________________________________________________
               FOR1: for($p=0; $p <= @poss_sub_dir_heads; $p++){ ## Default has '.' will make things like '././fam_8_8.pbla.gz'
                    $subd=$poss_sub_dir_heads[$p];               ## Also, the '<=' not '<' cures the same problem.
                    FOR2 : for($e=0; $e <  @possible_extensions; $e++){
                         $ext=$possible_extensions[$e];
                         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                         #  This makes all the possible lower upper case names
                         #______________________________________________________
                         for( $u=0; $u < @U_L_case; $u++){
                            for($v=0; $v <@name_types; $v++){
                               $each_seq_name=$name_types[$v];
                               if($U_L_case[$u]=~/U/){  $each_seq_name="\U$each_seq_name";
                               }else{                   $each_seq_name="\L$each_seq_name"; }

                               if(-s "$each_seq_name\.$ext"){
                                    push(@final_files, "$each_seq_name\.$ext" ) ;
                                    $found_search_prog_exention_used=$ext;
                                    $found_real_subdir_name=$subd; ## This is to report the name of the actual subd found
                                    $found_search_prog_exention_used=$ext;
                                    next FOR0
                               }elsif(-s "$each_seq_name\.$ext\.gz"){
                                    push(@final_files, "$each_seq_name\.$ext\.gz" ) ;
                                    $found_search_prog_exention_used=$ext;
                                    $found_real_subdir_name=$subd; ## This is to report the name of the actual subd found
                                    $found_search_prog_exention_used=$ext;
                                    next FOR0
                               }else{
                                    $file_wanted="\.\/$subd\/$each_seq_name\.$ext";
                                    if(-s $file_wanted){
                                        push( @final_files, $file_wanted);
                                        $found_real_subdir_name=$subd; ## This is to report the name of the actual subd found
                                        $found_search_prog_exention_used=$ext;
                                        next FOR0
                                    }elsif(-s "$file_wanted\.gz"){
                                        push( @final_files, "$file_wanted\.gz");
                                        $found_search_prog_exention_used=$ext;
                                        $found_real_subdir_name=$subd; ## This is to report the name of the actual subd found
                                        next FOR0;
                                    }
                               }
                            }
                         }
                    } # FOR2
               } # FOR1
               print @final_files, "\n";
         } # FOR0

         #print "\n# @final_files \n=============> $big_out_msp  \n\n";

         if(@final_files < 1){
              print "\n# convert_clu_to_msp :LINE no.: ", __LINE__, " ERROR: \@final_files is empty. Serious error\n";
              print "\n If you have sub dir which have more than 2 chars as names, you may increase the default 2 to 3 in the above\n";
              next;
         }
         $write_each_msp_to_disk='';

         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #  Check if small msp files have already made in previous steps
         #________________________________________________________________
         if($final_files[0]=~/(\S+)\.msp/){ ##  concatenate msp into big_msp
             $search_file_base=$1;
             print "\n# $search_file_base $found_real_subdir_name $found_search_prog_exention_used\n" if $verbose;

             if($final_files[0]=~/\S\.gz$/){
                 print "\n# $final_files[0] is gzipped \n";
                 system("gzip -d  $final_files[0]");
                 $final_files[0]=~s/\.gz//;
             }
             open(BIG_MSP_FILE, ">$big_out_msp");
             for($y=0; $y< @final_files; $y++){
                 open(SINGLE_MSP, "$final_files[$y]");
                 while(<SINGLE_MSP>){
                     print BIG_MSP_FILE $_;
                 }
             }
             close(BIG_MSP_FILE);
             close(SINGLE_MSP);
             push(@written_msp_files, $big_out_msp);

             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
             # Doing something for L option ($Lean_output)
             #___________________________________________________
             if($Lean_output and -d $found_real_subdir_name and $found_search_prog_exention_used){
                 for($y=0; $y< @seq_names; $y++){
                     unlink("$found_real_subdir_name\/$seq_names[$y]\.$found_search_prog_exention_used");
                     unlink("$found_real_subdir_name\/$seq_names[$y]\.$found_search_prog_exention_used\.gz");
                 }
             }elsif($Lean_output){
                 for($y=0; $y< @seq_names; $y++){
                     unlink("$seq_names[$y]\.$found_search_prog_exention_used");
                     unlink("$seq_names[$y]\.$found_search_prog_exention_used\.gz");
                 }
             }
         }else{
             if($write_each_msp_to_disk){
                  print "\# $0 : going to run open_sso_files with $write_each_msp_to_disk opt\n";
                  $big_out_msp=${&open_sso_files(\@final_files, $uppercase_seq_name, $write_each_msp_to_disk,
                                                                          "u=$upper_expect_limit", $new_format, $add_range, $add_range2, $big_out_msp, $over_write)};
                  if(-s $big_out_msp > 200){  print "\n# $0: SUCCESS to create $big_out_msp :) :) :-) :-) ?\n"; }
             }else{
                  print "\n# convert_clu_to_msp: I am running open_sso_files. \n";
                  @msp_hashes=@{&open_sso_files(\@final_files, $uppercase_seq_name, $write_each_msp_to_disk,
                                                 "u=$upper_expect_limit", $new_format, $add_range,
                                                 $add_range2, $big_out_msp, $over_write)};
                  &write_msp_files(@msp_hashes, $big_out_msp); ## concatenates all the hash ref to one
             }
         }
     }## end of  for($i=0; $i< @clusters; $i++){
     return(\@written_msp_files);
}# end of


#________________________________________________________________________________
# Title     : clu_to_sso_to_msp
# Usage     : &clu_to_sso_to_msp(\$clu);
# Function  : reads in a big single linkage cluster file(or normal cluster file)
#              and creates a big msp file which contains all the entries in the
#              cluster file (usually with the extension of sclu or clu)
#             This normally reads in xxxx.mso, xxxx.sso like files, but if the
#              corresponding  xxx.msp file already exists, it concatenates them to
#              make a bigger one.
# Example   :
# Keywords  : clu_2_sso_2_msp, cluster_to_msp, cluster_to_sso_to_msp
# Options   :
# Version   : 1.7
#--------------------------------------------------------------------------------
sub clu_to_sso_to_msp{
     my($i, $j, $k, $s, $u, $p, $m, $n, $y, @possible_extensions, @seq_names,
        @final_files, @U_L_case, $file, @file, @written_msp_files);

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # Opening cluster file (xx.clu)
     # %clus looks like this:  2-507     YGR041W YLR353W
     #                         3-308     YDR222W YDR346C YLR225C
     #                         2-184     YCL066W YCR040W
     #______________________________________________________________
     my $clu=${$_[0]} || $_[0];
     if($verbose){
         print "\n# clu_to_sso_to_msp : \"$clu\" is given
             and I am processing it with clu_to_sso_to_msp\n" if defined $clu;
     }
     my %clus=%{&open_clu_files(\$clu)};
     my @keys= keys %clus;
     my $num_of_cluster=@keys=@{&sort_by_cluster_size(\@keys)};

     print "# $0: clu_to_sso_to_msp: No. of cluster=$num_of_cluster after open_clu_files \n";

     &show_array(\@keys) if $verbose;
     &show_hash(\%clus) if $verbose;
     @possible_extensions=('msp', 'sso', 'msso', 'msso.gz',
              'pbla.gz', 'pbla', 'fsso', 'ssso', 'fso', 'out', 'prot.sso', 'prot.ts');
     @U_L_case=('\U', '\L');

     for($i=0; $i< @keys; $i++){
        my (@seq_names, @final_files, $clus_name, $big_out_msp, @msp_hashes);
        $clus_name=$keys[$i];
        unless($single_file_name=~/\S/){
            $big_out_msp="$clus_name\_cluster\.msp"; #<<<----- final output name
        }else{
            $big_out_msp=$single_file_name;
        }
        push(@written_msp_files, $big_out_msp); ## This is the output of this sub

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  If $clus_name.msp is already there, skip
        #_____________________________________________
        if( (-s $big_out_msp) > 100  and !$over_write ){
            print "\n# clu_to_sso_to_msp : $big_out_msp MSP file already exists, skipping\n";
            print "#    Use  \$over_write option \'o\' to start all over again or \n";
            print "#    delete clustering files like XX-XX_cluster.clu to go on\n";
            next ;
        }
        $num_of_seq_member=@seq_names=split(/ +/, $clus{$keys[$i]}); # @seq_names has (HIU001, HI002, HI333, MJ111, etc)
        print "# $0: clu_to_sso_to_msp: No. of seq member=$num_of_seq_member after split \n" if $verbose;

        FOR0: for($j=0; $j < @seq_names; $j++){
           my($sub_dir_head, $file_name_low, $file_name_up, $file_name_prot_low, @sub_dir_heads,
               $file_name_prot_up, $file_name_low_gz, $file_name_up_gz,
               $file_name_prot_low_gz, $file_name_prot_up_gz);

           $each_seq_name=$seq_names[$j];
           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           #  Here I take chars from the sequ names, as dirs have fragments of chars
           #_______________________________________________________________________________
           for($s=1; $s <=2 ; $s++){  ## here, number 2 indicates, I check single or 2 char sub dir names
               $sub_dir_head= substr($seq_names[$j], 0, $s);
               push(@sub_dir_heads, "\L$sub_dir_head") if (-d "\L$sub_dir_head" );
               push(@sub_dir_heads, "\U$sub_dir_head") if (-d "\U$sub_dir_head" );
           }
           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           #  Checking all the possible subdirectories to crop all the sso files
           #_______________________________________________________________________________
           FOR1: for($p=0; $p < @sub_dir_heads; $p++){
               $subd=$sub_dir_heads[$p];
               FOR2 : for($e=0; $e < @possible_extensions; $e++){
                    $ext=$possible_extensions[$e];
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    #  This makes all the possible lower upper case names
                    #______________________________________________________
                    for( $u=0; $u < @U_L_case; $u++){
                        if($U_L_case[$u]=~/U/){  $each_seq_name="\U$each_seq_name";
                        }else{                   $each_seq_name="\L$each_seq_name"; }

                        if(-s "$each_seq_name\.$ext"){   push(@final_files, "$each_seq_name\.$ext" ) ; next FOR0 }
                        elsif(-s "$each_seq_name\.$ext\.gz"){ push(@final_files, "$each_seq_name\.$ext\.gz" ) ; next FOR0 }
                        else{
                            $file_wanted="\.\/$subd\/$each_seq_name\.$ext";
                            if(-s $file_wanted){
                                push( @final_files, $file_wanted); next FOR0 }
                            elsif(-s "$file_wanted\.gz"){
                                push( @final_files, "$file_wanted\.gz");
                                next FOR0
                            }
                        }
                    }
               } # FOR2
           } # FOR1

        } # FOR0

        print "\n# @final_files \n=============> $big_out_msp  \n\n" if $verbose;

        if(@final_files < 1){
           print "\n# clu_to_sso_to_msp :LINE no.: ", __LINE__, " ERROR: \@final_files is empty. Serious error\n";
           print "\n If you have sub dir which have more than 2 chars as names, you may increase the default 2 to 3 in the above\n";
           next;
        }
        $write_each_msp_to_disk='w';

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #  Check if small msp files have already made in previous steps
        #________________________________________________________________
        if($final_files[0]=~/\.msp *$/){ ##  concatenate msp into big_msp
             open(BIG_MSP_FILE, ">$big_out_msp");
             for($y=0; $y< @final_files; $y++){
                 open(SINGLE_MSP, "$final_files[$y]");
                 while(<SINGLE_MSP>){
                     print BIG_MSP_FILE $_;
                 }
             }
             close(BIG_MSP_FILE);
             close(SINGLE_MSP);
             push(@written_msp_files, $big_out_msp);

        }else{
            if($write_each_msp_to_disk){
                 print "\# $0 : going to run open_sso_files with $write_each_msp_to_disk opt\n";
                 $big_out_msp=${&open_sso_files(\@final_files, $uppercase_seq_name, $write_each_msp_to_disk,
                               "u=$upper_expect_limit", $new_format, $add_range, $add_range2, $big_out_msp, $over_write)};
                 if(-s $big_out_msp > 200){  print "\n# $0: SUCCESS to create $big_out_msp :) :) :-) :-) ?\n"; }
            }else{
                 print "\n# clu_to_sso_to_msp: I am running open_sso_files. \n";
                 @msp_hashes=@{&open_sso_files(\@final_files, $uppercase_seq_name, $write_each_msp_to_disk,
                               "u=$upper_expect_limit", $new_format, $add_range, $add_range2, $big_out_msp, $over_write)};

                 &write_msp_files(@msp_hashes, $big_out_msp); ## concatenates all the hash ref to one
            }
        }
     }
     return(\@written_msp_files);
}# end of




#________________________________________________________________________
# Title     : read_dir_names_only
# Usage     : @all_dirs_list = @{&read_dir_names_only(\$absolute_path_dir_name, ....)};
# Function  : read any dir names and and then put in array. If no argument
#             for the target directory, it opens PWD automatically
#             You can specify the length of dir names to choose.
# Example   : @files=@{&read_dir_names_only('n', "s=1", '.')};
# Warning   : This does not report '.', '..'
#             Only file names are reported. Compare with &read_any_dir
# Keywords  : read_dir_only, get_dir_names, get_dir_names_only, get_subdir_names,
# Options   : n   for names only reading(not the full path) , default is full path
#             s=  for the size of dirs name. If you want all the dir names
#                   with a size of 1 char, s=1
# Returns   : one ref. of array.
# Argument  : takes one or more scaler references. ('.', \$path, $path, ... )
# Version   : 3.4
#--------------------------------------------------------------------
sub read_dir_names_only{
  my($in_dir, $i,$k, @possible_dirs, @chopped_pwd_path, @args,
	  @final_files, $full_dir, $pwd, $path,@read_files,
	  $size_of_dir_name);
  $pwd=`pwd`;
  chomp($pwd);
  $full_dir=1;
  @args=@_;

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Checking option
  #__________________________________________________
  for($k=0; $k < @args; $k++){
	 if(    $args[$k] eq 'n' or ${$args[$k]} eq 'n'){
	     $full_dir=0;
	     print "\n# read_dir_names_only: You put \'n\' option \n";
		 splice(@args, $k, 1); $k--;
	 }elsif( $args[$k] =~/s=(\d+)/ or ${$args[$k]} =~/s=(\d+)/){
		 $size_of_dir_name=$1;
		 print "\n# read_dir_names_only : You have put the size of dir names : $size_of_dir_name\n";
		 splice(@args, $k, 1); $k--;
	 }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # When no arg, this opens PWD automatically
  #_________________________________________________
  if(@args == 0){
	 if($full_dir==1){
		 $in_dir=$pwd;
	 }else{
		 $in_dir='.';
	 }
	 @final_files=@{&open_and_read_dir_names_only(\$in_dir)};
  }elsif(@args > 0){
	 for($k=0; $k < @args; $k++){
		if(!(ref($args[$k]))){    $in_dir=$args[$k];
		}elsif(ref($args[$k])){   $in_dir=${$args[$k]};    }

		if($in_dir ne '..' and $in_dir !~ /\// ){
			push(@final_files, @{&open_and_read_dir_names_only(\$in_dir)} );
		}elsif($in_dir eq '..' and $full_dir==1){
			print "\n# read_dir_names_only: \"..\" is given to open\n";
			@chopped_pwd_path=split(/\//, $pwd);
			pop(@chopped_pwd_path);
			$in_dir=join('/', @chopped_pwd_path);
			push(@final_files, @{&open_and_read_dir_names_only(\$in_dir)} );
		}elsif($in_dir eq '..'){
		    $in_dir = '..';
			push(@final_files, @{&open_and_read_dir_names_only(\$in_dir)} );
			for(@final_files){ $_=~s/\.//; }
		}
		##########  Main READING PART ##########

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Embedded subroutine
		#_________________________________________________
		sub open_and_read_dir_names_only{
			my ($i, @final_files);
			my $in_dir=${$_[0]};
			opendir(DIR1,"$in_dir");
			my @read_files = readdir(DIR1);
			if($size_of_dir_name){
				for($i=0; $i < @read_files; $i ++){
					unless(length($read_files[$i]) == $size_of_dir_name){
						next;
					}
					if($full_dir ==1){
						$read_files[$i]="$in_dir\/$read_files[$i]";
					}
					if( ($read_files[$i] !~ /\/\.\.?$/) && ( -d $read_files[$i]) ){
						 $read_files[$i]=~s/[\.\/]*//; ## removing ./ in front of dirs (in bash)
						 push(@final_files, "$read_files[$i]");
					}
				}
				return([@final_files]);
			}else{
				for($i=0; $i < @read_files; $i ++){
					if($full_dir ==1){
						$read_files[$i]="$in_dir\/$read_files[$i]";
					}
					if( ($read_files[$i] !~ /\/\.\.?$/) && ( -d $read_files[$i]) ){
						 $read_files[$i]=~s/[\.\/]*//; ## removing ./ in front of dirs (in bash)
						 push(@final_files, "$read_files[$i]");
					}
				}
				return([@final_files]);
			}
  		}
	 }
  }
  return([sort @final_files]);
}
#________________________________________________________________________________
# Title     : get_overlapping_seq_match_size
# Usage     : $ovlapsize=${&get_overlapping_seq_match_size($st1, $en1, $st2, $en2)
# Function  :
# Example   :
# Keywords  : CF: get_overlapping_range, get_overlapping_seq_match
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
#___________________________________________________________
# Title     : get_seq_fragments
# Usage     : @seq_frag=&get_seq_fragments(\%msf, @RANGE);
# Function  : gets sequence(string) segments with defined
#             ranges.
# Example   :
#  %test=('seq1', '1234AAAAAAAAAAAaaaaa', 'seq2', '1234BBBBBBB');
#  @range = ('1-4', '5-8');
#
#  %out = %{&get_seq_fragments(\%test, \@range)};
#  %out => (seq1_5-8   AAAAA
#           seq2_5-8   BBBBB
#           seq1_1-4    1234
#           seq2_1-4    1234 )
#
# Warning   :
# Keywords  : get_sequence_fragments,
# Options   : _  for debugging.
#             #  for debugging.
#             l=  for min seqlet length
#             r  for adding ranges in the seq names
#
# Returns   :
# Argument  :
# Version   : 1.8
#-------------------------------------------------------
sub get_seq_fragments{
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
   my $min_seqlet_size=10;
   my @vars=keys %vars;
   my $no_range_in_name=1;
   for($i=0; $i< @vars; $i++){
	   if($vars[$i] eq 'l'){
		  $min_seqlet_size=$vars{$vars[$i]};
	   }
   }
   if($char_opt=~/v/){ print "\n \$char_opt is $char_opt  @char_opt\n"; }
   if($char_opt=~/n/){ $no_range_in_name = 1 }
   if($char_opt=~/r/){ $no_range_in_name = 0 }

   print "\nget_seq_fragments \$no_range_in_name is $no_range_in_name \n";
   for($i=0; $i< @hash; $i++){
	 my (%out_frag, $frag_name, $range_start, $range_end, @out_hash);
	 my %seqs = %{$hash[$i]};
	 my @names = keys %seqs;
	 if(@names==1){
	    for($j=0; $j < @names; $j++){
		   my $seq_name = $names[$j];
		   my $seq = $seqs{$seq_name};
		   for($k=0; $k< @range; $k++){
			  my $range = $range[$k];
			  if($no_range_in_name==1){
				 $frag_name = "$seq_name";
			  }else{
			     $frag_name = "$seq_name\_$range";
			  }
			  #if(length($frag_name)>14 ){
			  #	 $frag_name ='x'."${j}_${range}";
		      #}
			  ($range_start, $range_end)=$range=~/(\d+\.?\d*)\-(\d+\.?\d*)/;
			  my $frag_len = $range_end-$range_start+1;
			  if($frag_len < $min_seqlet_size){
			     next;
			  }
			  my $fragment = substr($seq, $range_start-1, $frag_len);
			  $out_frag{$frag_name}=$fragment;
		   }
		}
		push(@out_hash,  \%out_frag);
	 }elsif(@names > 1){
	    for($k=0; $k< @range; $k++){
		  my %out_frag=();
	      my $range=$range[$k];
		  ($range_start, $range_end)=$range=~/(\d+\.?\d*)\-(\d+\.?\d*)/;
	      my $frag_len = $range_end-$range_start+1;
		  if($frag_len < $min_seqlet_size){
		     next;
		  }
	      for($j=0; $j < @names; $j++){
	         my $seq_name=$names[$j];
			 my $seq = $seqs{$seq_name};
		     if($no_range_in_name==1){
				 $frag_name = "$seq_name";
			 }else{
			     $frag_name = "$seq_name\_$range";
			 }
			 #if(length($frag_name)>15 ){
			 #	$frag_name ='x'."${j}_${range}";
		     #}
			 if($range_start==0){ $range_start++; } ## This is a bugfix
			 my $fragment = substr($seq, $range_start-1, $frag_len);
			 $out_frag{$frag_name}=$fragment;
		  }
		  push(@out_hash, \%out_frag);
		}
	 }
   }
   if(@out_hash > 1){ return(@out_hash)
   }elsif(@out_hash==1){ return($out_hash[0]) }
}
#________________________________________________________________________________
# Title     : print_clusfile_from_hash
# Usage     : &print_clusfile_from_hash(\%hash)
# Function  : To print out a file in cluster file format from an input hash containing the genes as keys and the cluster number as values.
# Example   :
# Keywords  : print_single_linkage_cluster, print_cluster_file
# Options   :
# Author    : Sarah A. Teichmann
# Version   : 1.2
#--------------------------------------------------------------------------------
sub print_clusfile_from_hash {
    my ($i, $j, $k, $gene, @subclus, %hash, $single_linkage_cluster);
    $single_linkage_cluster="single_linkage.sclu";
    open(SING, ">$single_linkage_cluster") or die "\n# $0: print_clusfile_from_hash: failed to open $single_linkage_cluster\n";
    if( @_==1 and ref($_[0]) eq 'HASH'){%hash=%{$_[0]};}

    my @clusters=values(%hash);

    @clusters=@{&remove_dup_in_array(\@clusters)};
    my @genes=keys(%hash);
    for ($i=0; $i<@genes; $i++){
        my $clus=$hash{"$genes[$i]"};
        push(@{"$clus"},$genes[$i]);
        next;
    }

    my (%sizes);
    for ($i=0; $i<@clusters; $i++){
        my $cluster=$clusters[$i];
        @{"$cluster"}=@{&remove_dup_in_array(\@{"$cluster"})};
        my $size=@{"$cluster"};
        $sizes{"$size"}.="$cluster\n";
        next;
    }


    my @clus_sizes=keys(%sizes);
    @clus_sizes=sort {$a<=>$b} (@clus_sizes);

    for ($i=0; $i<@clus_sizes; $i++){
        my $clus_size=$clus_sizes[$i];
        unless ($clus_size>1){next;}
        print SING "Cluster size $clus_size\n";
        my @subclus=split(/\n/,$sizes{"$clus_size"});
        for ($j=0; $j<@subclus; $j++){
            my $clus=$subclus[$j];
            print SING "Cluster $clus\n";
            for ($k=0; $k<@{"$clus"}; $k++){
                $gene=${"$clus"}[$k];
                print SING " 1 1 $gene\n";
                next
                }
            next;
       }
       next;
    }
    close SING;
    return(\$single_linkage_cluster);
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
# Returns   :
# Argument  :
# Version   : 1.1
#--------------------------------------------------------------
sub get_internal_dup_in_a_cluster{
	$cluster_line=$_[0] || ${$_[0]};
	my($i, $j, @seq, %out, $seq_name);
	my $overlap_factor=40;
	my $min_inside_dom_size=40;
	@seq=split(/ +/, $cluster_line);  ## These sequence are single seq with different regions
	@seq= map{$_->[0]} sort{$a->[1] cmp $b->[1] or $a->[2] <=> $b->[2] }
			             map {/^((\S+)_(\d+)\-(\d+) *.*)$/ && [$1, $2, $3, $4]} @seq;

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
	#@out=@{&remove_similar_seqlets(\@temp, "f=2")};
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
    print "\n# import_ENV_vars: ALL the ENV settings are imported to $0 program\n" if $verbose;
}

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
# Title     : sort_by_cluster_size
# Usage     : @out=@{&sort_by_cluster_size(\@input_line_array)};
# Function  : it sorts by the 1st digit before '-'  as in 2-183_cluster, 2-140_cluster,
#               etc.
# Example   :
# Keywords  : sort_by_columns, sort_by_text_columns, sort_by_column_numerically
#             sort_by_pattern
# Options   :
# Version   : 1.2
#----------------------------------------------------------------------------
sub sort_by_cluster_size{
   my (@in, @M, $col);
   if(@_ < 1  ){ print "\n# FATAL: sort_by_cluster_size needs 1 argument\n"; exit }
   if(ref $_[0] eq 'ARRAY'){ 	  @in = @{$_[0]};      }else{ 	  @in = @_;    }
   $col=0;
   @in= map {$_->[0]} sort { $a->[1] <=> $b->[1] } map { [$_, ($_=~/^(\S+)\-/)[$col] ] } @in;
   return(\@in);
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
               print "#  They are: @seqlets\n";
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
#______________________________________________________________________________
# Title     : make_seq_index_file
# Usage     : @idx_files_made=@{&make_seq_index_file(\@file)};
# Function  : creates xxxx.fa.idx file and makes a link to pwd. If @file contains
#              names with .idx extension already, it will not put another idx
#              index to it.
# Example   :
# Keywords  : make_fasta_seq_index_file, create_seq_index_file, make_idx_file,
#             create_idx_file, create_seq_idx_file, make_index_file, create_index_file
#             make_sequence_index_file, create_sequene_index_file
# Options   :
# Version   : 1.4
#----------------------------------------------------------------------------
sub make_seq_index_file{
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
	my(@index_files_made, $fasta_db_input, $fasta_db_idx, %index);
	print "\n# make_seq_index_file : input \@file was @file\n";

	for($i=0; $i< @file; $i++){
		$fasta_db_input=$file[$i];
		if($fasta_db_input !~/\S+\.idx$/){
			$fasta_db_idx="$fasta_db_input.idx";
		}

		open(FASTA_DB, "$fasta_db_input");
		open(FASTAIDX, ">$fasta_db_idx");

 	    print FASTAIDX "# fasta_index for $fasta_db_input\n";
		while(<FASTA_DB>){
			if(/^\> {0,4}(\S+) */){
				$index{$1}=tell(FASTA_DB);
				print FASTAIDX "\n$1 $index{$1}";
			}
		}
		close(FASTA_DB);
		close(FASTAIDX);
		if(-s $fasta_db_idx){
			print "\n# The size of $fasta_db_idx is more than 0, looks O.K. \n";
			push(@index_files_made, $fasta_db_idx);
			system("ln -s $fasta_db_idx .");
		}else{
		    print "\n# The size of $fasta_db_idx is less than 0, ERROR??\n";
		}
	}
	if(@file < 2){
	   return( \$fasta_db_idx );
	}else{
	   return(\@index_files_made);
	}
}
#______________________________________________________________
# Title     : sort_words_in_string
# Usage     :
# Function  : sort words in strings sperated by ' ' or "\n"
# Example   :
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
# Keywords  : sort_words_in_sequences, sort_sequences_in_string,
#             sort_strings_in_string, sort_string_by_words, sort_elements_in_string
# Options   : _  for debugging.
#             #  for debugging.
# Returns   :
# Argument  :
# Version   : 1.1
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
#________________________________________________________________________________
# Title     : msp_single_link_hash
# Usage     : %hash=%{&msp_single_link_hash(\@msp_files, E-value);
# Function  : To make a hash with all the genes in the msp files as the keys,
#             which are linked at or below the E-value threshhold,
#             with the values denoting the cluster number
# Example   :
# Keywords  : single_linkage, msp_single_linkage, msp_single_linkage_hash
#              make_msp_single_link_hash
# Options   :
# Author    : Sarah A. Teichmann with thanks to Alex Bateman,
#               sat@mrc-lmb.cam.ac.uk, jong@salt2.med.harvard.edu
# Version   : 1.4
#--------------------------------------------------------------------------------
sub msp_single_link_hash{
    my (@msp_files, $i, $j, $k, $e_val, $gene_1, $gene_2,
        @mspcont, $gene_1, $gene_2, $E_cut, %hash, $array);

    $E_cut=0.001; ## Default Eval cut

    if( @_ =>2 and ref($_[0]) eq 'ARRAY'){
        @msp_files=@{$_[0]};
        $E_cut=${$_[1]} || $_[1];
    }else{
        print "Subroutine msp_single_link_hash takes one input array and the E-value as its arguments!!" && die;
    }
    print "\n# (i) msp_single_link_hash : \$E_cut is $E_cut \nwith @msp_files\n" if $verbose;

    for ($i=0; $i<@msp_files; $i++){
        if ($msp_files[$i]=~/\S+\.msp$/){
            open(MSP, "$msp_files[$i]");
            @mspcont=<MSP>;
            close(MSP);
        }elsif ($msp_files[$i]=~/\S+\.msp\.gz$/){
            @mspcont=`gunzip -c $msp_files[$i]`;
        }

        $array++;
        for ($j=0; $j<@mspcont; $j++){
            if ($mspcont[$j]=~/^\d+ +(\S+) +\S+ +\d+ +\d+ +(\S+) +\d+ +\d+ +(\S+)/){
                $e_val=$1;
                unless($e_val<=$E_cut){next;}
                $gene_1=$2;
                $gene_2=$3;
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
                # Removing REGION information as Sarah's can not handle regions
                #_______________________________________________________________
                if($gene_1=~/(\S+)_\d+\-\d+/){ $gene_1=$1 }
                if($gene_2=~/(\S+)_\d+\-\d+/){ $gene_2=$1 }
                if ($gene_1 eq $gene_2){next;}
                if ( ! $hash{"$gene_1"} and $gene_1){
                    $hash{"$gene_1"}="$array";
                    push (@{"$array"}, $gene_1);
                }
                if (! $hash{"$gene_2"} and $gene_2){
                    $hash{"$gene_2"}="$array";
                    push (@{"$array"}, $gene_2);
                }
                if ($hash{"$gene_1"}==$hash{"$gene_2"}){next;}
                if ( $hash{"$gene_1"} gt  $hash{"$gene_2"}){
                    push (@{"$hash{$gene_2}"}, @{"$hash{$gene_1}"});
                    for ($k=0; $k<@{"$hash{$gene_2}"}; $k++){
                         $hash{${"$hash{$gene_2}"}[$k]}=$hash{"$gene_2"};
                         next;
                    }
                    next;
                }
                if ( $hash{"$gene_2"} gt  $hash{"$gene_1"}){
                    push (@{"$hash{$gene_1}"}, @{"$hash{$gene_2}"});
                    for ($k=0; $k<@{"$hash{$gene_1}"}; $k++){
                         $hash{${"$hash{$gene_1}"}[$k]}=$hash{"$gene_1"};
                         next;
                    }
                    next;
                }
                next;
            }else{ next; }
        }
    } ## enf of for loop
    return (\%hash);
} ## msp_single_link_hash


#________________________________________________________________________
# Title     : write_fasta_seq_by_seq
# Usage     : &write_fasta_seq_by_seq(\%hash, [$extension], [\$output_filename]);
# Function  : accepts one hash of multiple sequences and writes many files
#             of single sequences by using the names as file names.
#             If $extension is provided, it writes an output as in
#             the below example (seq1_sc.fasta). If not, it just attach
#             'fa' to files.
#             This needs, hash of 'name', 'actual sequence as value'
# Example   : with >xxxx
#                  ASDFASDFASDFASDFASDFASDFASDF
#                  >yyyy
#                  ASDFASDFASDFASDFASDFASDFSDAFSD
#
#             You will get two files (xxxx.fa, yyyy.fa)
# Keywords  : write_each_fasta, write_single_fasta, write_fasta_single
#             single_fasta_write, write_fasta_files_seq_by_seq,
#             write_single_fasta_files,
# Options   : can specify extension name.
#             e  for checking fasta file exists or not and skipps if so
#             r for rename the sequences so that Clustalw would not complain with 10 char limit
#               so result wuld be:  0 ->ASDFASDF, 1->ASDFASFASF, 2->ADSFASDFA
# Returns   : nothing. default OUTPUT file name is '$key.fa' !!
# Version   : 2.0
#--------------------------------------------------------------------
sub write_fasta_seq_by_seq{
	 my ($i, $exists_opt, $rename_seq_opt, $out_file_name_given);
	 for($i=0; $i< @_; $i++){
		if($_[$i]=~/e$/){
		   $exists_opt=1;
		   splice(@_, $i, 1);
		   $i--;
		}elsif($_[$i]=~/r$/){
		   $rename_seq_opt='r';
		   splice(@_, $i, 1);
		   $i--;
		}elsif( $_[$i] =~/\.fa/ or -e $_[$i] ){
		   $out_file_name_given=1;
		   $out_file_name = $_[$i];
		   splice(@_, $i, 1);
		   $i--;
		}elsif( ref ($_[$i]) eq 'SCALAR'){
		   if( ${$_[$i]} =~/\.fa/ or -e ${$_[$i]} ){
		      $out_file_name_given=1;
		      $output_file=${$_[$i]};
		      splice(@_, $i, 1);
		      $i--;
		   }
		}
	 }
	 my(%temp_hash, $key, $output_file);
	 my(%input)     =%{$_[0]};
	 my($extension) =${$_[1]} || $_[1];
	 for $key (keys %input){
		my %temp_hash=();
		$temp_hash{$key}=$input{$key};
		if (($key=~ /\_$extension$/)||($#_ == 0)){
			$output_file = "$key\.fa";
		}else{
			$output_file = "$key\_$extension\.fa";
		}
		if( ($exists_opt==1)&&(-e $output_file)){
		   print "# write_fasta_seq_by_seq: File $output_file exists, NO write\n";
		}elsif( $out_file_name_given == 1){
		   &write_fasta(\%temp_hash, \$output_file, $rename_seq_opt);
		}else{
		   &write_fasta(\%temp_hash, \$output_file, $rename_seq_opt);
		}
	 }
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
# Title     : reverse_sequences
# Usage     : %out = %{&rev_sequence_one_hash(\%input_seq_hash, \%hash2,...)};
# Function  : gets ref. of strings, reverses the elems.
# Example   :
# Warning   :
# Keywords  : reverse_sequence, reverse_sequence_hash, rev_sequence_hash
# Options   :
# Returns   : one or more hash references.
# Argument  : hash, eg(1, 'skdfj', 2, 'kdfjkdj', 3, 'kdfjk');
#             Input example:
#             ..
#             >HI0256
#             FLSANVLPIAPIINGGRTAVDNITQSVSDKPFVKDIGTKIKEAIALSKYSTQPQYISTTN
#             >HI0094
#             DILRTFVKMETGLKFPKKFKLKANLALFMNRRNKRPDTIMTAVADAGQKISEAKLNTTAK
#             ..
#
#             Output example: (Reversed :-)
#             ..
#             >HI0256_rv   <<-- note the added extension
#             ALDJFLKAJFJALSDJFLAJSLFJAKLSDFJLASJDFLAJSLDFJASJDFLJSDFJSDLJ
#             >HI0094_rv
#             LASJDFLKAJFJALSDJFLKSDJLFAJLKDJFLASJDFLKDFJKDJFKDJFKDJFKJDLJ
#             ..
#
# Version   : 1.2
#--------------------------------------------------------------------
sub reverse_sequences{
	my(%rev_hash, @rev_hash_refs, $name, $name_with_ext, $i);
	for($i=0; $i < @_; $i++){
	    my %in_hash = %{$_[$i]};
		my @keys    = keys %in_hash;
		for $name (@keys ){
		    $name_with_ext = "$name\_rv";
			$rev_hash{$name_with_ext} = reverse($in_hash{$name});
		}
		push(@rev_hash_refs, \%rev_hash);
	}
	if(@rev_hash_refs ==1){ return($rev_hash_refs[0]);}
	else{ return(@rev_hash_refs);}
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
# Title     : open_clu_files
# Usage     : %clus=%{&open_clu_files(\$input)};
# Function  :
# Example   : Clu file eg)
#
#  Cluster 7360103
#    1  1 SLL1058         7-255       2   Origin: 3   736   Sub:3
#    1  1 MJ0422          17-283      2   Origin: 3   736   Sub:3
#    1  1 HI1308          3-245       2   Origin: 3   736   Sub:3
#
# Warning   : You MUST NOT delete '# options : ..' entry
#              as it is read  by various subroutines.
#              This automatically converts lower to upper letters
# Keywords  : open_cluster_files,
# Options   : _  for debugging.
#             #  for debugging.
#             b  for to get just names ($simple_clu_reading)
#             r  for adding ranges in the names
#             U  for makeing sequence names upppercase
#
# Returns   : a ref of hash of $clus{"$clus_size\-$id"}.=$m."\n";
#             Actual content:
#             3-133 => 'HI00111 HI00222 MG1233 '
# Argument  :
# Version   : 1.9
#--------------------------------------------------------------
sub open_clu_files{
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
   my($simple_clu_reading, $possible_range, $add_ranges,
	  $id, $name_range, %clus, $found, $upper_case_seq_name);
   my $file=$file[0];
   if($char_opt=~/b/){ $simple_clu_reading= 'b' };
   if($char_opt=~/U/){ $upper_case_seq_name='U'; };

   my $clus_size=1;
   open(CLU, "$file");
   while(<CLU>){
	  if($simple_clu_reading=~/\S/){ ## to get just names
		  if(/^ *\d+ +\d+ +\d+ +\d+ +\d+/){  ## To skip the very first summary columns
			 next;
		  }elsif(/^ *#/ ){ next;
		  }elsif(/^ *\d+ +\d+ +(\S+) +(\S+)/){
			 $seq_name=$1;
			 $possible_range=$2;
			 if($2=~/\d+\-\d+/ and $char_opt=~/r/){
				$name_range="$seq_name\_$possible_range";
				$clus{$name_range} = $name_range;
			 }else{
			    $clus{$seq_name}=$seq_name;
			 }
		  }
	  }else{
		  if(/^ *\d+ +\d+ +\d+ +\d+ +\d+/){  ## To skip the very first summary columns
			 next;
		  }elsif(/^ *#/ ){
			 next;
		  }elsif(/^ *Cluster +size +(\d+)/i ){
			 $clus_size=$1;
			 $found=1;
		  }elsif(/^ *Cluster +([_\d]+) *size:? *(\d+)/i){  # to match 'Cluster 14313'  or  'Cluster 234_1234_1'
			 $id  =$1;
			 $found=1;
			 $clus_size=$2; # if defined($2);
		  }elsif(/^ *Cluster +[number]* *([\w]+)/i){  # to match 'Cluster 14313'  or  'Cluster 234_1234_1'
			 $id  = $1;
			 $found=1;
		  }elsif(($found==1)&&(/^ *\S* *\S* *(\S+)\.prot\,? *.*/)){ ## this is to correct MP genome names
			 $m=$1;
             if($upper_case_seq_name){
                 $clus{"$clus_size\-$id"}.="\U$m ";
             }else{
                 $clus{"$clus_size\-$id"}.="\U$m ";
             }
		  }elsif(($found==1)&&(/^ *(\d+) *\d* *(\S{2,32}) *(\S*)/)){          # general clu match
			 $clus_size=$1 unless ($clus_size);
			 $m=$2;
			 $possible_range=$3;
			 if($2=~/\d+\-\d+/ and $char_opt=~/r/){
				$name_range="$m\_$possible_range";
                if($upper_case_seq_name){
                    $clus{"$clus_size\-$id"}.="\U$name_range ";
                }else{  $clus{"$clus_size\-$id"}.="$name_range "; }
			 }else{
                if($upper_case_seq_name){
                    $clus{"$clus_size\-$id"}.="\U$m ";
                }else{  $clus{"$clus_size\-$id"}.="$m ";  }
			 }
		  }
	  }
   }
   return(\%clus);
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
#             remove_same_array_elements, remove_redundancy, remove_redundant_elements
#             remove_duplication_in_array
# Options   :
#   s  for sorting the array output
# Returns   : one or more references.
# Argument  : one or more refs for arrays or one array.
# Version   : 1.6
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


#________________________________________________________________________________
# Title     : read_machine_readable_sso_lines
# Usage     : @out_refs=@{&read_machine_readable_sso_lines(\@SSO, $get_alignment,
#                           $create_sso, $upper_expect_limit,$new_format, $lower_expect_limit,
#                           $attach_range_in_names, $attach_range_in_names2)};
# Function  :
# Example   :
# Keywords  : read_m10_sso_lines read_msso_lines
# Options   : a c r r2 n
#             u= for upper E value limit
#             l= for lower E value limit
# Category  :
# Version   : 1.5
#--------------------------------------------------------------------------------
sub read_machine_readable_sso_lines{
   my ($upper_expect_limit, $lower_expect_limit)=(50,0);
   my (%match, @out_refs, $query_found, $query_sq_stop, $query_sq_statrt, $match_found,
      $match_seq, $match_found2, $i, $j,$match_found3, $overlap, $sw_score,
      $match_sq_stop, $match_seq2, $sw_ident, $name_range, $query_seq,
      $al_display_start, $match_seq_count);
   for($i=0; $i< @_; $i++){
       if($_[$i]=~/u=(\S+)/){    $upper_expect_limit=$1 }
       elsif(ref($_[$i]) eq 'ARRAY'){ @SSO=@{$_[$i]};   }
       elsif($_[$i]=~/l=(\S+)/){ $lower_expect_limit=$1 }
       elsif($_[$i]=~/^c *$/){     $create_sso = 'c'; print "\n# read_machine_readable_sso_lines: \$create_sso is set"; }
       elsif($_[$i]=~/^a *$/){     $get_alignment='a'; }
       elsif($_[$i]=~/^ *r *$/){     $attach_range_in_names='r' }
       elsif($_[$i]=~/^r2 *$/){    $attach_range_in_names2='r2' }
       elsif($_[$i]=~/^n *$/){     $new_format='n' }
   }
   $attach_range_in_names='r';

   print "\n# read_machine_readable_sso_lines : You put PARSEABLE form of sso file";
   for($j=0; $j< @SSO; $j++){
	  if($SSO[$j]=~/\>\>\> *(\S+)\,? +(\d+) +/){  ## >>>  line
               $query_found=1;  $query_seq_leng=$2;  ## Ignoring the $1, as file name can be different from rea seq names
                   $j+=8;
	  }elsif( $query_found==1 and $SSO[$j]=~/\>\>(\w[\w\-\.]+)([\.prot\,\:]*) */ ){ ##
               $match_found=1;
               $match_seq_count++;
               $al_display_start=0;
               if(length($2)>0){  print "\n# read_machine_readable_sso_lines: Seq name has this special char \"$2\". I ignore it"; }
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               #  Changing the CASE according to the option
               #_____________________________________________
               if($uppercase_seq_name eq 'U'){
                       $match_seq="$1"; $match_seq="\U$match_seq";  ## make it uppercase
               }elsif($lowercase_seq_name eq 'L'){
                       $match_seq="$1"; $match_seq="\L$match_seq"; ## make it lowercase
               }else{ $match_seq="$1"; } ## make it uppercase
               next;
	  }elsif($match_found and $SSO[$j]=~/^\; +\w+_expect\:? +(\S+)/){
               #~~~~~~~~~~~~~~~~~~~~~~~
               # Filtering by E val
               #_______________________
               if( $1 > $upper_expect_limit or $1 < $lower_expect_limit ){
                       $match_found=0; next;
               }else{ $expect =$1; }
	  }elsif($match_found and $SSO[$j]=~/^ *\; +sw_score *\: +(\S+)/i){  $sw_score =$1;
	  }elsif($match_found and $SSO[$j]=~/^\; +sw_ident\: +(\S+)/){  $sw_ident =$1;
	  }elsif($match_found and $SSO[$j]=~/^ *\; +sw_overlap\: +(\S+)/){  $overlap=$1;
	  }elsif($match_found and $SSO[$j]=~/^ *\>(\w[\w\-\.]+)([\.prot\,\:]*) *[\d+]*/){
               $match_found2=1;        $match_found=0;
               if( length($2)>0 ){  print "\n# read_machine_readable_sso_lines: Seq name has this special char \"$2\". I ignore it"; }
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               #  Changing the CASE according to the option
               #_____________________________________________
               if($uppercase_seq_name eq 'U'){
                       $match_seq2="$1"; $match_seq2="\U$match_seq2"; ## make it uppercase
               }elsif($lowercase_seq_name eq 'L'){
                       $match_seq2="$1";  $match_seq2="\L$match_seq2"; ## make it lowercase
               }else{ $match_seq2="$1";  }
               $query_seq=$match_seq2;
	  }elsif($match_found2==1 and $SSO[$j]=~/\; +sq_len\: +(\S+)/){
		     $query_sq_len=$1;
	  }elsif($match_found2==1 and $SSO[$j]=~/\; +al_start\: +(\S+)/){
		     $query_sq_start=$1;
	  }elsif($match_found2==1 and $SSO[$j]=~/\; +al_stop\: +(\S+)/){
		     $query_sq_stop=$1;
	  }elsif($SSO[$j]=~/\; +al_display_start/ and $al_display_start < 1){
             $al_display_start ++;
	  #------------------------------------------------------------
	  }elsif($match_found2 and $SSO[$j]=~/\>(\w[\w\-\.]+)([\.prot\,\:]*) *[\d+]*/){
             $match_found3=1; $match_found2=0;
             if(length($2)>0){  print "\n# open_sso_files: Seq name has this special char \"$2\". I ignore it"; }
	  }elsif($match_found3 and $SSO[$j]=~/\; +sq_len\: +(\S+)/){
		     $match_sq_len=$1;
	  }elsif($match_found3 and $SSO[$j]=~/\; +al_start\: +(\d+)/){
		     $match_sq_start=$1;
	  }elsif($match_found3 and $SSO[$j]=~/\; +al_stop\: +(\d+)/){
             $match_sq_stop=$1;
	  }elsif($match_found3 and $SSO[$j]=~/\; +al_display_start/){
             $match_found3=0;          $al_display_start++;
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # To prevent Null SWscore and Eval
             #___________________________________________________
             if($expect=~/^$/){
                 $expect='0.0';
                 if((!$sw_score and $expect==0) or $sw_ident > 0.8){
                     $sw_score=1000;
                 }
             }
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # adding the offset for names with ranges
             #__________________________________________________
             if($query_seq=~/^\S+_(\d+)\-(\d+)/){ $query_sq_start +=$1-1; $query_sq_stop +=$1-1;  }

             #~~~~~~~~~~~~~~~~~~~~~~~~~
             # Attaching the ranges  (under NO e option)
             #_________________________

             if($attach_range_in_names){
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  # Checks margin opt and adds it
                  #__________________________________
                  if($margin=~/\d+/){
                      if($match_sq_start < $margin){  $match_sq_start=1;
                      }else{          $match_sq_start-=$margin;   }
                      $match_sq_stop += $margin;
                  }
                  $name_range="$match_seq\_$match_sq_start\-$match_sq_stop";

                  #~~~~~~~~ If 'rr' opt is set, put ranges for both target and match seqs ~~~~~~~
                  if($attach_range_in_names2 and $query_seq !~/^\S+_(\d+)\-(\d+)/){
                      $query_seq="$query_seq\_$query_sq_start\-$query_sq_stop";
                  }
                  if($original_target_seq=~/\S+/){ $query_seq=$original_target_seq } # for PVM version out
                  if($new_format=~/n/){  # under NO e option
                      $match{$name_range}=
                         sprintf("%s %s %s %s %s %s %s %s %s\n",
                         $query_seq, $query_sq_start, $query_sq_stop, $sw_score, $expect, $sw_ident,
                         $match_sq_start, $match_sq_stop, $name_range);
                  }else{
                      $match{$name_range}=
                         sprintf("%-5s %-9s %-6s %-4s %-5s %-30s %-4s %-5s %s\n",
                         $sw_score, $expect, $sw_ident, $query_sq_start, $query_sq_stop, $query_seq,
                         $match_sq_start, $match_sq_stop, $name_range);
                  }
             }else{
                 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 # Checks margin opt and adds it
                 #__________________________________
                 if($margin=~/\d+/){
                     if($match_sq_start < $margin){  $match_sq_start=1;
                     }else{                          $match_sq_start-=$margin; }
                     $match_sq_stop += $margin;
                 }
                 if($original_target_seq=~/\S+/){ $query_seq=$original_target_seq } # for PVM version out
                 if($new_format=~/n/){
                     $match{$match_seq}=
                        sprintf("%s %s %s %s %s %s %s %s %s\n",
                        $query_seq, $query_sq_start, $query_sq_stop, $sw_score, $expect, $sw_ident,
                        $match_sq_start, $match_sq_stop, $match_seq);
                 }else{
                    $match{$match_seq}=sprintf("%-5s %-9s %-6s %-4s %-5s %-30s %-4s %-5s %s\n",
                       $sw_score, $expect, $sw_ident, $query_sq_start, $query_sq_stop, $query_seq,
                       $match_sq_start, $match_sq_stop, $match_seq);
                 }
             }
	  }elsif($get_alignment and $al_display_start==1 and $SSO[$j]=~/^([\w\-]+) *$/){
		  ${"match_alignment\_$match_seq_count"}{$match_seq2} .= $1;
	  }elsif($get_alignment and $al_display_start==2 and $SSO[$j]=~/^([\w\-]+) *$/){
		  ${"match_alignment\_$match_seq_count"}{"$match_seq"} .= $1;
	  }elsif($get_alignment and $SSO[$j]=~/^ *\;al_cons\:/){
		  $al_display_start=0;
		  my %temp=%{"match_alignment\_$match_seq_count"};
		  push(@out_refs, \%temp );
		  %{"match_alignment\_$match_seq_count"}=();
	  }
   } ## <-- for($j=0; $j< @SSO; $j++)

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # If create sso option is set, it creates SSO files(array input case)
   #________________________________________________________________________
   if( $create_sso and !$get_alignment){
	   open (SSO2, ">$query_seq\.msso");
	   print SSO2 @SSO, "\n";
	   print "\n# read_machine_readable_sso_lines : $query_seq\.msso file  is created by \"c\" opt ";
	   close SSO2
   }
   unless($get_alignment){
	   push(@out_refs, \%match);
   }
   return(\@out_refs);
}



#________________________________________________________________________________________
# Title     : merge_sequence_in_msp_file
# Usage     :
# Function  :
# Example   : INPUT: (MSP file) ===>
#  59     2.6        47    64     d2pia_3        10    30     d1erd___10-30
#  161    1.1e-07    24    91     d2pia_3        16    85     d1frd___16-85
#
#  722    0          1     106    d1put__        1     106    d1put___1-106
#  66     4.9        2     68     d1put__        43    106    d2lbp___43-106
#  69     1.3        12    49     d1put__        81    120    d1cgo___81-120
#
#  60     3.3        13    38     d1frd__        32    57     d1orda1_32-57
#  65     1.7        21    58     d1frd__        40    69     d2mtac__40-69
#
#   ==== OUTPUT ===>
#    d1frd___1-98 d1frd___1-98_1-98 d1frd___16-85 d2pia_3_24-91_24-91
#    d1frd___16-85_16-85 d2pia_3_24-91
#    d1put___1-106 d1put___1-106_1-106
#    d2pia_3_1-98 d2pia_3_1-98_1-98
#
# Keywords  : mergr_seq_in_msp_file, merge_sequence_in_msp, merge_sequences_in_msp_file
# Options   :
#  $dynamic_factor =  y by y -y   # adjusting factor value dynamically(more seq higher factor)
#  $short_region   =  S by S -S  # taking shorter region overlapped in removing similar regions
#  $large_region   =  L by L -L  # taking larger  region overlapped in removing similar regions
#  $average_region =  A by A -A # taking average region overlapped in removing similar regions
#
# Version   : 2.9
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
        $length_thresh, $merge, $factor, $evalue, $score, $dynamic_factor, $score_match,
        $eval_match, $query_seq, $query_start, $query_stop, $match_seq, $match_start,
        $short_region, $large_region, $average_region, $original_clu_size, $match_stop);
    $factor=$default_factor=7; #~~~~ default connection factor U, 7 means 70% now!
    $length_thresh=30;
    $evalue=10;
    $score =75;
    if(@file < 1){ print "\n# (E) merge_sequence_in_msp_file needs at least 1 MSP file\n"; exit }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Following changes the defaults with given parameters
    #_____________________________________________________________
    if($char_opt=~/z/i){       $optimize='z';    ## This will cause using remove_similar_seqlets than remove_dup_in_array !
    }if($char_opt=~/m/){       $merge='m';
    }if($char_opt=~/y/){       $dynamic_factor='y';
    }if($char_opt=~/v/){       $verbose='v';
    }if($char_opt=~/S/){       $short_region='S';
    }if($char_opt=~/L/){       $large_region='L';
    }if($char_opt=~/A/){       $average_region='A';
    }if($vars{'T'}=~/\d+/){    $length_thresh=$vars{'T'};
    }if($vars{'f'}=~/\S+/){    $factor=$vars{'f'};  ## Here I give a generous $factor !
    }if($vars{'s'}=~/\d+/){    $score = $vars{'s'};
    }if($vars{'e'}=~/\d+/){    $evalue= $vars{'e'}; }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #  Just to inform what parameters have been chosen
    #_____________________________________________________________
    if($verbose){
        print "\n# (1) merge_sequence_in_msp_file : default \$score      : $score";
        print "\n#                                : default \$evalue     : $evalue";
        print "\n#                                : used \$length_thresh : $length_thresh";
        print "\n#                                : default \$factor     : $default_factor";
        print "\n#                                : used    \$factor     : $factor";
        print "\n#                                : \$dynamic_factor     : $dynamic_factor\n";
    }

    for($c=0; $c< @file; $c++){
       open(MSP, "$file[$c]");
       $base=${&get_base_names($file[$c])};
       $clu_out="$base\_F${factor}.clu"; # <-- This is the most important output. Sarah's program will process this
       $sat_out="$base\_F${factor}.sat";
       print "# (2) merge_sequence_in_msp_file : processing $file[$c] for $clu_out\n" if $verbose;
       my @msp1=<MSP>;

       for($i=0; $i< @msp1; $i++){
          #~~~~~~~~~~ Include range or NOT in the seq name ~~~~~~~~~~~~~~~~~~~~~~~~~~`
          # %temp_hash is just to get the chunk of MSP block. As msp file uses empty line as a delimiter
          #____________________________________________________________________________
          if($char_opt=~/r/){
             if($msp1[$i]=~/^ *(\d+) +(\S+) *\S* +(\d+) +(\d+) +(\S+)[ \t]+(\d+)[ \t]+(\d+)[ \t]+(\S+)/){
                  $score_match=$1;
                  $eval_match=$2;
                  $query_seq=$5;
                  $query_start=$3;
                  $query_stop=$4;
                  $match_seq=$8;
                  $match_start=$6;
                  $match_stop=$7;
                  if($query_seq=~/\S+_\d+\-\d+$/){ $new_seq1=$query_seq }else{ $new_seq1="$query_seq\_$query_start\-$query_stop"; }
                  if($match_seq=~/\S+_\d+\-\d+$/){ $new_seq2=$match_seq }else{ $new_seq2="$match_seq\_$match_start\-$match_stop"; }

                  if($new_seq1 eq $new_seq2){
                     print "# (3.0) merge_sequence_in_msp_file: Skipped = $msp1[$i]"  if $verbose;
                     next
                  };

                  if($score_match < $score or $eval_match > $evalue){
                      print "# (3.1) merge_sequence_in_msp_file: Skipped = $msp1[$i]"  if $verbose;
                      next
                  };
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  # Modifying $msp1[$i] line !!!
                  #______________________________
                  $msp1[$i]=sprintf("%s %-3s %s %s %s %s %s %s",
                                $score_match, $eval_match, $query_start, $query_stop, $new_seq1, $match_start,
                                $match_stop, $new_seq2);
                  $temp_hash{$query_seq}.="$msp1[$i]\n";
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
        print "\n# ### \$factor: $factor\n";
        $factor += (log($original_clu_size)*5)/10 - 1; ## This is a simplistic.
        if($factor > 9.5){ $factor=9.5 } # this is the very upper limit for any factor.
        print "\n# ### \$factor: $factor\n";
    }

    if($verbose){
        print "# (4) merge_sequence_in_msp_file: The msp chunks used are:\n";
        for($i=0; $i< @msp_chunks; $i++){
            print "$msp_chunks[$i]\n";
        }
    }

    for($i=0; $i< @msp_chunks; $i++){
        print "\n# (5) merge_sequence_in_msp_file: Processing eash chunk with merge_sequence_in_msp_chunk\n" if $verbose;
        @arr=@{&merge_sequence_in_msp_chunk($msp_chunks[$i], $verbose, $optimize,
                  "$merge", "e=$evalue", "s=$score", "f=$factor", "T=$length_thresh",
                  $short_region, $large_region, $average_region)};
        push(@all_seqlets,  @arr);
    }

    #~~~~~~~~~ sorting inner sequences in strings ~~~~~~~~~
    #______________________________________________________
    @all_seqlets=@{&sort_words_in_string(@all_seqlets)}; ## This speeds up about 2 times !!!

    #~~~~~~~ Sort by the _digit-  in seqlet names ~~~~~~~~~
    @all_seqlets= map{$_->[0]} sort{$a->[1] cmp $b->[1] or $a->[2] <=> $b->[2]  }
                  map {/^ *((\S+)_(\d+)\-(\d+).*)/ && [$1, $2, $3, $4]} @all_seqlets;
    if(  $verbose){
        print "\n# (6) merge_sequence_in_msp_chunk: Showing the very final result before merging\n";
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
        @all_seqlets=@{&remove_dup_in_array(\@all_seqlets)};

    }else{
        @all_seqlets=@{&remove_dup_in_array(\@all_seqlets)};
    }
    return(\@all_seqlets);
}

#______________________________________________________________________________
# Title     : run_fasta_sequence_search
# Usage     : $gzipped_msp_file=${&run_fasta_sequence_search("a=$algorithm",
#                        "O=$out_file_msp_name", "File=$temp_file_name", "e=$E_val",
#                        "DB=$sequence_DB", "k=$k_tuple", "$machine_readable")};
#
# Function  :
# Example   :
# Keywords  : run_ssearch_sequence_search, do_fasta_sequence_search
# Options   :
#             Query_seqs=  for enquiry sequences eg)  "Query_seqs=$ref_of_hash"
#             DB=   for target DB  "DB=$DB_used"
#             File= to get file base(root) name.  "File=$file[0]"
#             i= to get file base(root) name. same as File=
#             m  for MSP format directly from FASTA or Ssearch result than through sso_to_msp to save mem
#             s  for the big single output (msp file output I mean)
#             s= for the single big msp file name
#             O= for Out file name, same as s=
#             o  for overwrite existing xxxx.fa files for search
#             c  for create SSO file (sequence search out file)
#             d  for very simple run and saving the result in xxxx.gz format in sub dir starting with one char
#             r  for reverse the query sequence
#             R  for attaching ranges of sequences
#             k= for k-tuple value. default is 1 (ori. FASTA prog. default is 2)
#             u= for $upper_expect_limit
#             l= for $lower_expect_limit
#             a= for choosing either fasta or ssearch algorithm
#             d= for defining the size of subdir made. 2 means it creates
#                    eg, DE while 1 makes D
#             d  for $make_gz_in_sub_dir_opt, putting resultant sso files in gz format and in single char subdir
#             D  for $make_msp_in_sub_dir_opt, convert sso to msp and put in sub dir like /D/, /S/
#             n  for new format to create new msp file format with sso_to_msp routine
#          PVM=  for PVM run of FASTA (FASTA only)
#             M  for machine readable format -m 10 option
#             M= for machine readable format -m 10 option
#             N  for 'NO' do not do any processing but, do the searches only.
#       FILE_AGE for defining the age of file in days to be overwritten.
# Author    : Jong Park, jong@salt2.med.harvard.edu, for commercial use, ask me.
# Category  :
# Version   : 1.1
#------------------------------------------------------------------------------
sub run_fasta_sequence_search{
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
    my($algorithm, $E_val, $sequence_DB, $k_tuple, $machine_readable,$create_sso,
       @msp_from_temp, $gzipped_msp_file, @temp, $temp_file_name, $add_range,
       $out_file_msp_name, $input_file_name, $out_file_sso_name,
       $gzipped_sso_file, $machine_readable);
    $algorithm='fasta';    $out_file_msp_name='default_out_file.msp';
    $E_val=1;              $sequence_DB=$ENV{'NRDB_FASTA'};
    $k_tuple=1;            $machine_readable='M';
    $upper_expect_limit=1; $lower_expect_limit=0;
    $Score_thresh=100;

    if($vars{'a'}=~/\S+/){ $algorithm          = $vars{'a'}            };
    if($vars{'O'}=~/\S+/){ $out_file_msp_name  = $vars{'O'}            };
    if($vars{'e'}=~/\S+/){ $E_val              = $vars{'e'}            };
    if($vars{'u'}=~/\d+/){ $upper_expect_limit = $vars{'u'}            };
    if($vars{'l'}=~/\d+/){ $lower_expect_limit = $vars{'l'}            };
    if($vars{'k'}=~/\d+/){ $k_tuple            = $vars{'k'}            };
    if($vars{'t'}=~/\d+/){ $Score_thresh       = $vars{'t'}            };
    if($vars{'r'}=~/\S+/){ $add_range          = 'r'                   };
    if($vars{'s'}=~/\S+/){ $single_big_msp     = 's'                   };
    if($vars{'DB'}=~/\S+/){            $sequence_DB=$vars{'DB'} ;
        if(-s $sequence_DB){
        }elsif(-s "../$sequence_DB"){  $sequence_DB= "../$sequence_DB"
        }elsif(-s "../../$sequence_DB"){  $sequence_DB= "../../$sequence_DB";
        }else{
            print "\n# (E) do_sequence_search: You set DB param, but I can\'t find $sequence_DB\n";
            exit;
        }
    }else{  print "\n# (E) do_sequence_search: I need DB param defined, sorry, aborting\n"; }

    if($vars{'FILE'}=~/\S+/){ $input_file_name = $vars{'FILE'}; };
    if($vars{'File'}=~/\S+/){ $input_file_name = $vars{'File'}; };
    if($vars{'FILE_AGE'}=~/\S+/){ $age_in_days_of_out_file= $vars{'FILE_AGE'};  };
    if($vars{'Query_seqs'}=~/\S+/){ %seq_input = %{$vars{'Query_seqs'}}};
    if($vars{'Query'}=~/\S+/){      %seq_input = %{$vars{'Query'}}};
    if($vars{'u'}    =~/\S+/){ $E_val          = $vars{'u'}            };
    if($vars{'PVM'}  =~/\S+/){ $PVM_FASTA_run  = $vars{'PVM'}; print "\n# PVM opt is set\n";     };
    if($vars{'M'}  =~/\S+/){ $machine_readable = $vars{'M'};           };

    if($char_opt=~/r/){    $add_range          = 'r' }
    if($char_opt=~/o/){    $over_write         = 'o' }
    if($char_opt=~/c/){    $create_sso         = 'c' }
    if($char_opt=~/s/){    $single_big_msp     = 's'; print "\n# Single file opt is set\n"; }
    if($char_opt=~/M/){    $machine_readable   = 'M' }
    if($char_opt=~/N/){    $No_processing      = 'N'; $create_sso='c'; }
    print "\n# (i) run_fasta_sequence_search: $algorithm,$out_file_msp_name,$input_file_name,$E_val,DB=$sequence_DB,$k_tuple,$machine_readable\n";

    if(!$algorithm or !$input_file_name or !$sequence_DB){
        print "\n# (E) One of \$sequence_DB, \$input_file_name, \$algorithm is missing\n";
        exit;
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # (1) Machine readable opt for fasta and ssearch set ??
    #____________________________________________________
    if($machine_readable=~/M/){
       print "\n# (i) Running $algorithm -m 10 -H  -E $E_val $input_file_name $sequence_DB $k_tuple\n";
       @temp=`$algorithm -m 10 -H  -E $E_val $input_file_name $sequence_DB $k_tuple`;
    }else{
       @temp=`$algorithm -H -E $E_val $input_file_name $sequence_DB $k_tuple`;  # -H is for NO histogram
    }
    print "\n# (i) \@temp has ",scalar(@temp), " lines @temp!\n" if $verbose;

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # (2) @temp is parsed directly by open_sso_files to make MSP hashes
    #_______________________________________________________________________
    @msp_hashes_from_temp = @{&open_sso_files(\@temp, $add_range, "u=$upper_expect_limit", "l=$lower_expect_limit")};
    if(@msp_hashes_from_temp < 1){
        print "\n# (E) do_sequence_search : Error, something is wrong with open_sso_files, LINE=", __LINE__, "\n";
        exit;
    }else{   print "\n# (i) Good, \@msp_from_temp has ",scalar(@msp_hashes_from_temp), " lines !\n";   }
    @msp_from_temp= values %{$msp_hashes_from_temp[0]};
    print "\n# (i) @msp_from_temp\n" if $verbose;

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # (3) Writing MSP file, and gzipping if possible
    #________________________________________________
    open(MSP, ">$out_file_msp_name") or die "\n# (E) $out_file_msp_name could not be open\n";
    for(@msp_from_temp){    print MSP $_;  }
    close MSP;
    $gzipped_msp_file=${&compress_files_by_gzip($out_file_msp_name)};

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # (4) making SSO file if asked
    #__________________________________________________
    if($create_sso){
       $seq_name=${&get_base_names($input_file_name)};
       if($algorithm=~/ssearch/){
          $out_file_sso_name="$seq_name\.ssso";
          $out_file_sso_name_gz="$seq_name\.ssso\.gz";
       }elsif($algorithm=~/fasta/){
          $out_file_sso_name="$seq_name\.fsso";
          $out_file_sso_name_gz="$seq_name\.fsso\.gz";
       }
       open(SSO, ">$out_file_sso_name");
       for(@temp){  print SSO $_;  }; close (SSO);
       if($machine_readable){
           $out_file_msso_name="$seq_name\.msso";
           &cp( $out_file_sso_name, $out_file_msso_name);
       }
       @gzipped_sso_files=@{&compress_files_by_gzip($out_file_sso_name, $out_file_msso_name)};

       if(@gzipped_sso_files > 0){  print "\n# (i) @gzipped_sso_files are created"; }
    }
    if(-s $gzipped_msp_file and !(-s $out_file_msp_name)){
       print "\n# (i) $gzipped_msp_file is found and will be returned";
       return(\$gzipped_msp_file);
    }elsif(-s $out_file_msp_name){
       print "\n# (i) $out_file_msp_name is found and will be returned";
       return(\$out_file_msp_name);
    }
}## end of run_fasta_sequence_search


#______________________________________________________________________________
# Title     : check_if_defined
# Usage     : $defined=&check_if_defined($var, $file);
# Function  : checks if all the args are defined
# Example   :
# Keywords  :
# Options   :
# Author    : jong@salt2.med.harvard.edu
# Category  :
# Version   : 1.0
#------------------------------------------------------------------------------
sub check_if_defined{
    my($i, $input_arg_count);
    for($i=0; $i< @_; $i++){
       $input_arg_count++;
       if(-s $_[$i]){
          print "\n# (INFO) check_if_defined: $input_arg_count does not exist as a file in PWD\n";
       }elsif(! $_[$i]){
          print "\n# (INFO) check_if_defined: $input_arg_count is not defined\n";
          return(0);
       }
    }
    return(1);
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
#   $take_only_the_last_iteration=l by l
#   $accumulative_hits_eval_thresh= by e=
#   $genome_seq_only=g      by g
#   $nrdb_seq_only=n        by n
#   $evalue_thresh=         by E=
#   $Accumulate_matches=A   by A -A
#   $Lean_output=L          by L -L  # to remove search output to unclutter
#
# Author    : Sarah Teichmann and Jong Park, jong@salt2.med.harvard.edu
# Version   : 3.9
#--------------------------------------------------------------------------------
sub convert_bla_to_msp{
   my($i, $j, $k, @lines, $match_string_count,  $line_count, $query_string_count,
      $match_length, $Lean_output,
      $lines, $duplicated_match_count, $new_sorted_name, $sorted_name, $verbose,
      $pdbd_seq_only, $entry_found, $which_iteration, $report_only_the_best,
      $genome_seq_only, $all_seq, $header_found, $accumulative_hits_eval_thresh,
      $take_only_the_last_iteration, $original_query, $nrdb_seq_only,
      $get_the_final_iteration, $read_entry_lines, $verbose, $Accumulate_matches);
   my $match_leng_thresh=10;
   ### This localization is critial NOT my, as I use a sub which relies on this
   local(%hash_out, %accumulative_hits, $file, $score, $score_ori, $evalue,
         $evalue_ori, $seq_id, $query_range_start, $query_range_stop,
         $query, $match_string_start, $match_string_stop, $matched,
         $read_point_found);
   $duplicated_match_count=0;

   my $evalue_thresh=$accumulative_hits_eval_thresh=1; ## default eval threshes
   $query='query_seq'; ## default query seq name, to avoid blank name

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Processing the input arguments to get file and options etc
   #_____________________________________________________________
   for (@_){
       if(ref $_ eq 'ARRAY'){ @lines =@{$_};
       }elsif( ref $_ eq 'SCALAR' and -s ${$_} ){ $file=${$_};
       }elsif( -s $_ ){            $file=$_;
       }elsif(/^ *d *$/){          $pdbd_seq_only='d'; $all_seq=''; $genome_seq_only='';
           print "\n $0: convert_bla_to_msp,  You set \$pdbd_seq_only option, I will skip others.\n";
       }elsif(/^ *[i=]*(\d+) *$/){ $which_iteration=$1;
       }elsif(/^ *b *$/){          $report_only_the_best='b';
       }elsif(/^ *a *$/){          $all_seq='a'; $genome_seq_only=''; $pdbd_seq_only=''; $nrdb_seq_only='';
       }elsif(/^ *g *$/){          $genome_seq_only='g'; $all_seq=''; $pdbd_seq_only='';$nrdb_seq_only='';
       }elsif(/^ *n *$/){          $nrdb_seq_only='n'; $all_seq=''; $pdbd_seq_only=''; $genome_seq_only='';
       }elsif(/^ *l *$/){          $take_only_the_last_iteration='l';
       }elsif(/^ *v *$/){          $verbose='v';
       }elsif(/^ *L *$/){          $Lean_output='L';
       }elsif(/e=(\S+)/){          $accumulative_hits_eval_thresh=$1;
       }elsif(/E=(\S+)/){          $evalue_thresh=$1;
       }elsif(/A$/){               $Accumulate_matches='A';  }
   }

   unless($which_iteration){  $get_the_final_iteration=1 }
   print "\n\n# (W) convert_bla_to_msp: NO \$Accumulate_matches opt \'A\' set. Sure??\n\n", chr(7) if !$Accumulate_matches;

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
   # If the input file is gzipped, uncompress it to text file and then open
   #__________________________________________________________________
   if($file=~/\.gz *$/){
       print "\n# (INFO) >>>>> Running gunzip to open $file";
       open(BLA_FILE, "gunzip -c $file|") || die "\n# $0: Failed to open $file\n";
       print "\n# (INFO) $0: bla_to_msp, input file was $file\n";
       if($file=~/^([de]*\d\d*\w\w\w\w\w)\./){         $query=$1;
       }
   }else{
       open(BLA_FILE, "$file") || die "\n# $0: Failed to open $file\n";
       print "\n# (i) $0: convert_bla_to_msp, input file was $file\n";
       if($file=~/^([de]*\d\d*\w\w\w\w\w)\./){         $query=$1;
       }
   }
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # UP to NOW is frivalous option handling stuff
   #_______________________________________________________


   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # (1) Main reading in .pbla file (or any extension)
   #____________________________________________________________________________
   while(<BLA_FILE>){
      $line_count++;      $lines=$_;  ## putting $_ to $lines var
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # (1) Matching "BLASTP 2.0.4 [Feb-24-1998]"
      #____________________________________________________________
      if(!$header_found and $lines=~/^ *BLASTP +\d/){  # blastp header as   BLASTP 2.0.2 [Sep-3-1997] or BLASTP 2.0.4 [Feb-24-1998]
           $header_found=1;
           print "\n# (1) \$header_found becomes $header_found\n";
      }elsif($header_found and $lines=~/^ *BLASTP +\d/){  print "\n# (E) \$header_found is already SET and I found \"BLASTP\", ERROR!!!\n"; exit }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # (1.1) If I reach the end of the opened file, I use &put_msp_lines_to_hash_from_bla sub to write the final msp line and finish
      #________________________________________________________________________________________________________________
      if( eof ){
           if( $read_point_found <= $which_iteration  or  $get_the_final_iteration){
               print "     # (i) <<<< The end of file reached, writing $sorted_name \n\n";
               @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out, \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                        $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                        $match_string_stop, $read_point_found, $accumulative_hits_eval_thresh,
                                        $take_only_the_last_iteration, $accumulative_hits_eval_thresh, $evalue_thresh)};
               %hash_out=         %{$out_from_put_msp_lines[0]};
               %accumulative_hits=%{$out_from_put_msp_lines[1]};
               $read_point_found= $out_from_put_msp_lines[2];
               last;
           }
      }

      if($lines=~/^ *$/ or $lines=~/^  +Length +\= +\d+ *$/){      next     }  ## skipping some junk lines

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # (2) Extracting query seq name(this is the only place to get it)
      #____________________________________________________________
      if($lines=~/^ *Query= +(\S+)/){ $query=$original_query=$1;    next    }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~
      # (3) 'Searching......done'  line indicates new search step(iteration)
      #_________________________________________________________________________
      if( $lines=~/^ *Searching\.\.+[done]/i ){
          $read_point_found++;
          if($verbose){
             print "\n==========================================================";
             print "\n# (3.1) Searching........ line found, so \$read_point_found is $read_point_found\n";
          }
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # (3.2) This is just error checking
          #__________________________________________________________
          if(!$header_found){ ##  header is  'BLASTP 2.0.4 [Feb-24-1998]',  the very first line in bla output
              print chr(7); print chr(7); print chr(7);print chr(7);print chr(7);
              print "\n# (ERROR)\'Sequences producing\' line is matched while BLASTP header is not found\n";
              print "\n# (ERROR) convert_bla_to_msp:  It means the file is truncated !!!!\n"; exit;
          }

          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
          #  (3.3) Following is the KEY part for controlling iteration
          #__________________________________________________________
          if( $which_iteration and $read_point_found < $which_iteration){
               print "\n# (INFO) skipped, \$which_iteration: ($which_iteration), \$read_point_found: ($read_point_found)" if $verbose;
               $match_string_count=$query_string_count=$score=$evalue=$seq_id=$score_ori=$evalue_ori='';
               $query_range_stop=$query_range_start=$match_string_stop=$msp_line=$new_sorted_name='';
               $entry_found=$duplicated_match_count=0;

               if( !$Accumulate_matches){
                   %hash_out=(); ## this is to remove any discarded pairs in the iteration
               }
               print "\n# (INFO) ===> New iteration ====\$read_point_found: $read_point_found, \$which_iteration:$which_iteration\n";
               next;
          }elsif( $which_iteration and $read_point_found == $which_iteration){
               $read_entry_lines=1;  next;
          }elsif( $which_iteration and $read_point_found >  $which_iteration){
               @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out, \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                        $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                        $match_string_stop, $read_point_found, $accumulative_hits_eval_thresh,
                                        $take_only_the_last_iteration, $accumulative_hits_eval_thresh, $evalue_thresh)};
               %hash_out=         %{$out_from_put_msp_lines[0]};
               %accumulative_hits=%{$out_from_put_msp_lines[1]};
               $read_point_found= $out_from_put_msp_lines[2];
               last;
          }elsif(!$which_iteration){
               #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
               # (3.4) Default situation
               #____________________________________________________________
               print "\n# (WARN) You did not set \$which_iteration option \n\n" if $verbose;
               if($read_point_found > 1){
                       print "\n (3.3) Writing the last entry $sorted_name BEFORE next Searching........ line\n";
                       @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out, \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                                $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                                $match_string_stop, $read_point_found, $accumulative_hits_eval_thresh,
                                                $take_only_the_last_iteration, $accumulative_hits_eval_thresh, $evalue_thresh)};
                       %hash_out=         %{$out_from_put_msp_lines[0]};
                       %accumulative_hits=%{$out_from_put_msp_lines[1]};
                       $read_point_found= $out_from_put_msp_lines[2];
               }
               $match_string_count=$query_string_count=$score=$evalue=$seq_id=$score_ori=$evalue_ori='';
               $query_range_stop=$query_range_start=$match_string_stop=$msp_line=$new_sorted_name='';
               $entry_found=$duplicated_match_count=0;
               if( !$Accumulate_matches){  %hash_out=(); $entry_found=0; $duplicated_match_count=0;     }
               $read_entry_lines=1; ## this is set by 'Searching......' line
               next;
          }
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # (4) '> xxxx '  New sequence entry, '>' starts
      #__________________________________________________________
      elsif($read_entry_lines and $lines=~/^\> *(\S+)/){
          $temp_match=$1;

          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # (4.0) To get the pdbd seq names only 'dxxx__' sort of thing
          #________________________________________________________
          if($pdbd_seq_only and ($temp_match !~/^pdb_\S+/ and $temp_match !~/^[cde]\d\S+/)  ){
              $entry_found=0; print "\n# NOT pdb seq\n";    next;
          }elsif($genome_seq_only and $temp_match !~/^gn_\S+/){
              $entry_found=0; print "\n# NOT genome seq\n"; next
          }elsif($nrdb_seq_only and $temp_match !~/^nr_\S+/){
              $entry_found=0; print "\n# NOT nrdb\n";      next
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
          # (4.1) This is the DEFAULT
          #_____________________________________________________________
          }else{ ## This is default and equivalent to have the $all_seq option on.
              $entry_found=1;
              if($match_string_count){ ## $match_string_count is incremented only by 'Sbjct' line
                    print "              (4.1) Writing $sorted_name ent to \%hash_out\n" if $verbose;
                    @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out, \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                             $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                             $match_string_stop, $read_point_found, $accumulative_hits_eval_thresh,
                                             $take_only_the_last_iteration, $accumulative_hits_eval_thresh, $evalue_thresh)};
                    %hash_out=         %{$out_from_put_msp_lines[0]};
                    %accumulative_hits=%{$out_from_put_msp_lines[1]};
                    $read_point_found= $out_from_put_msp_lines[2];
                    $match_string_count=0;
                    $duplicated_match_count=0;
                    print "              (4.1) reset \$match_string_count\n\n" if $verbose;
              }else{ print "              (4.1) \$match_string_count or \$duplicated_match_count is not right\n" if $verbose; }

              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              #  Only with new seq entry, I count the pair occurrances
              #__________________________________________________________________
              $query=$original_query; $query_string_count='';
              $matched=$temp_match; ## this should be here, after if
              $sorted_name=join(' ', sort($query, $matched) );
              print "\n# (4) >>> NEW \$sorted_name is  $sorted_name <<------ \$entry_found ($entry_found)\n" if $verbose;
          }
          $match_string_count=0;
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # (5) Matching  Score =  325 bits (824), Expect = 6e-89           << 2 >>
      #_________________________________________________________________
      elsif( ($entry_found and $lines=~/^[\t ]*Score[\t ]*\=[\t ]*(\S+)[\t ]*bits +\(\S+\)\,[\t ]*Expect[\t ]*=[\t ]*(\S+)/i)
          or ($entry_found and $lines=~/^[\t ]*Score[\t ]*\=[\t ]*(\S+)[\t ]*bits.+\,[\t ]*Expect[\t ]*=[\t ]*(\S+)/i)){

          $score_ori=$1;
          $evalue_ori=$2;
          if($evalue_ori=~/^e\-\d\d\d/){ $evalue_ori="1".$evalue_ori; } ## bug fix for short eval in blast distribution

          if($match_string_count){ # $match_string_count is increased when Sbjct word is found
              if($evalue > $evalue_thresh){ $evalue=$evalue_ori; $score=$score_ori; next }
              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
              # When Only the first match(best evalue) is required, write msp line and reset $entry_found var
              #_________________________________________________________________________________________________
              if($report_only_the_best){
                  print "      (5)  \$report_only_the_best is set\n" if $verbose;
                  @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out, \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                           $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                           $match_string_stop, $read_point_found, $accumulative_hits_eval_thresh,
                                           $take_only_the_last_iteration, $accumulative_hits_eval_thresh, $evalue_thresh)};
                  %hash_out=         %{$out_from_put_msp_lines[0]};
                  %accumulative_hits=%{$out_from_put_msp_lines[1]};
                  $read_point_found= $out_from_put_msp_lines[2];
                  $entry_found=0; next;
              }else{
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`~`
                  # duplicated match count means, query matched more than one region of a match seq
                  #__________________________________________________________________________________
                  $duplicated_match_count++;
                  $sorted_name="$sorted_name $duplicated_match_count";
                  print " ====(5) Multiple region for \"$new_sorted_name\" is found =========== e= $evalue\n" if $verbose;
                  @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out, \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                           $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                           $match_string_stop, $read_point_found, $accumulative_hits_eval_thresh,
                                           $take_only_the_last_iteration, $accumulative_hits_eval_thresh, $evalue_thresh)};
                  %hash_out=         %{$out_from_put_msp_lines[0]};
                  %accumulative_hits=%{$out_from_put_msp_lines[1]};
                  $read_point_found= $out_from_put_msp_lines[2];
              }
              $score=$score_ori; $evalue=$evalue_ori;
          }else{
              print "     (5) \$match_string_count is not set NO write \$evalue_ori $evalue_ori\n" if $verbose;
              $evalue=$evalue_ori; $score=$score_ori;
          } ## to next line

          sub reset_all_the_vars{
             print "            !!!!  Reseting all the vars !!!!\n" if $verbose;
             $query_string_count=$score=$evalue=$seq_id=$query_range_stop=$query_range_start='';
             $match_string_stop=$msp_line=$new_sorted_name=$match_string_count='';
          }
      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
      # (6) Matching   Identities = 158/158 (100%), Positives = 158/158 (100%)    ,
      #____________________________________________________________________________________
      elsif( $entry_found and $lines=~/^ *Identities += +\S+\/(\S+) +\( *(\S+) *\%\)/i){
          $query_string_count=$match_string_count=0;
          $seq_id=$2/100;
          $match_length=$1;
          if($match_length < $match_leng_thresh){
              print "     (6) \$match_leng_thresh $match_leng_thresh > \$match_length $match_length" if $verbose;
              $entry_found=0;
              $match_string_count=1;
              next;
          }else{
              print "     (6) $sorted_name : ABOVE leng thresh. \$seq_id= $seq_id, \$match_length= $match_length\n" if $verbose;
          }
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # (7) Matching 'Query: 2 GIRAATSQEINELT..' line    ,
      #_________________________________________________________________
      elsif($entry_found and $lines=~/^ *Query\:?[\t ]+(\d+) +\D+ +(\d+)/){
          $query_string_count++;
          $query_line_found=1;
          if($query_string_count==1){      $query_range_start=$1;   $query_range_stop =$2;
          }elsif($query_string_count > 1){ $query_range_stop=$2;     }
          print "        (7) Query: line found: $query\_$query_range_start\-$query_range_stop\n" if $verbose;
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # (8) Matching 'Sbjct: 2 GIRAATSQEINELT..' line
      #_________________________________________________________________
      elsif($entry_found and $query_line_found and $lines=~/^ *Sbjct\:? +(\d+) +[\w\-]+ +(\d+)/i){
          $match_string_count++;
          $subject_line_found=1;
          if($match_string_count==1){      $match_string_start=$1;
                                           $match_string_stop =$2;
          }elsif($match_string_count > 1){ $match_string_stop=$2;      }
          print "        (8) Sbjct: line found: $temp_match\_$match_string_start\-$match_string_stop\n" if $verbose;
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # (9) Matching '   Database: ' line    ,                << END >>
      #_________________________________________________________________
      elsif( ($entry_found and  $lines=~/^ +Database: +\S+/) or eof){ # the very last write
          if($evalue > $evalue_thresh){ last
          }else{
              print "        <<<< The end of file reached, writing $sorted_name\n" if $verbose;
              @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out, \%accumulative_hits, $query,$matched,$evalue, $score, $seq_id,
                                       $sorted_name, $query_range_start, $query_range_stop,$match_string_start,
                                       $match_string_stop, $read_point_found, $accumulative_hits_eval_thresh,
                                       $take_only_the_last_iteration, $accumulative_hits_eval_thresh, $evalue_thresh)};
              %hash_out=         %{$out_from_put_msp_lines[0]};
              %accumulative_hits=%{$out_from_put_msp_lines[1]};
              $read_point_found= $out_from_put_msp_lines[2];
              last;
          }
      }
   }
   close(BLA_FILE);
   unless( $take_only_the_last_iteration){
       print "\n# >> ACCUMULATIVE HITS are reported as you did not set \$take_only_the_last_iteration opt!!\n";
       %hash_out=(%hash_out, %accumulative_hits);
   }
   print "\n Congratulations!\n";
   print ">>>>>>>---- Finished reading in xxxx.pbla file Now sorting \%hash_out ----<<<<<<<<\n\n";
   &show_hash(\%hash_out) if $verbose;
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
}


#______________________________________________________________________________
# Title     : put_msp_lines_to_hash_from_bla
# Usage     : @out_from_put_msp_lines=@{&put_msp_lines_to_hash_from_bla(\%hash_out,
#                                        $query,$matched,$evalue, $score, $seq_id,
#                                        $sorted_name, $query_range_start,
#                                        $query_range_stop,$match_string_start,
#                                        $match_string_stop, $read_point_found,
#                                        $accumulative_hits_eval_thresh,
#                                        $take_only_the_last_iteration)};
# Function  :
# Example   :
# Keywords  :
# Options   :
# Author    : jong@salt2.med.harvard.edu,
# Category  :
# Version   : 1.2
#------------------------------------------------------------------------------
sub put_msp_lines_to_hash_from_bla{
     my (@finale_out, $sorted_name, $msp_line, $evalue, $score, $matched, $seq_id, $query_range_start,$accumulative_hits_eval_thresh,
         $query_range_stop, $query, $match_string_start, $match_string_stop, $read_point_found, %hash_out, %accumulative_hits, $evalue_thresh);
     %hash_out=%{$_[0]};         %accumulative_hits=%{$_[1]};
     $query=$_[2];               $matched=$_[3];
     $evalue=$_[4];              $score=$_[5];
     $seq_id=$_[6];              $sorted_name=$_[7];
     $query_range_start=$_[8];   $query_range_stop =$_[9];
     $match_string_start=$_[10]; $match_string_stop=$_[11];
     $read_point_found=$_[12];   $accumulative_hits_eval_thresh=$_[13];
     $take_only_the_last_iteration=$_[14];
     $accumulative_hits_eval_thresh=$_[15];
     $evalue_thresh=$_[16];

     $query  ="$query\_$query_range_start\-$query_range_stop";

     if($matched !~/^\S+\_\d+\-\d+ *$/){         $matched="$matched\_$match_string_start\-$match_string_stop";
     }elsif($matched =~/^(\S+)\_\d+\-\d+ *$/){   $matched="$1\_$match_string_start\-$match_string_stop";     }

     if($score=~/\S/ and $evalue=~/\S/ and $match_string_start=~/\S/ and $evalue_thresh > $evalue){
         $msp_line=sprintf("%-6s %-9s %-5s %-5s %-5s %-32s %-5s %-5s %-38s %-3s\n",
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
         if(!$take_only_the_last_iteration and $evalue <= $accumulative_hits_eval_thresh ){
            if($accumulative_hits{$sorted_name}){
               if($accumulative_hits{$sorted_name}=~/^[\t ]*\S+[\t ]+(\S+)[\t ]/){
                   if($evalue < $1){
                       $accumulative_hits{$sorted_name}=$msp_line;   }   }
            }else{ $accumulative_hits{$sorted_name}=$msp_line;     }
         }
     }else{     }
     @finale_out=(\%hash_out, \%accumulative_hits, $read_point_found, $query, $matched, $evalue, $score, $seq_id, $sorted_name,
                  $query_range_start, $query_range_stop, $match_string_start, $match_string_stop  );
     return(\@finale_out);
}



#______________________________________________________________________________
# Title     : do_psi_blast_search
# Usage     : &do_psi_blast_search(\@files, "d=$source_DB_file",
#                     "i=$input_seq_file",
#                     $over_write,
#                     $make_msp_in_sub_dir_opt);
# Function  :
# Example   :
# Keywords  :
# Options   :
#    $source_DB_file=           by d= s=
#    $input_seq_file=           by i=
#    $Eval_limit=               by E=
#    $iteration_limit=          by j=
#    $step_evalue=              by h= e=
#    $over_write=o              by o
#    $make_msp_in_sub_dir_opt=D by D
#    $Lean_output=L             by L -L
# Author    :
# Version   : 1.7
#------------------------------------------------------------------------------
sub do_psi_blast_search{
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
    $|=1;
    my(@final_msp_out, $iteration_limit, $source_DB_file, $input_seq_file,$over_write,
       $step_evalue, $Eval_limit, %seq_input, $base, @seq_names, $final_out_pbla_file,
       $make_msp_in_sub_dir_opt, $sub_dir_size, $answer_for_self_self,
       $defined_all_ok, $ENV_BLASTDB_var, $Lean_output);

    my $blast_env_var='BLASTDB';


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Default parameters for blastpgp
    #_____________________________________________
    $Eval_limit=100;
    $iteration_limit=20;
    $step_evalue=0.0005;
    $sub_dir_size      =2;  # the default char number taken from seq name to make sub dirs


     if($ENV{'BLASTDB'}=~/^ *$/){
         print "\n# (i) Your ENV variable setting  \"BLASTDB\" is not set\n";
         print "\n# (i) BLASTDB is for the directory for various blast DBs\n";
         print "\n# (i) You can set it by \"export BLASTDB=/your/dir/for/db\" in bash shell\n";
         print "\n# (i) or \"set BLASTDB /your/dir/for/db\" in cshell or in .bashrc or .cshrc file\n";
         print "\n# (i) Do you want to set the directory path right NOW?(y/n)\n";
         $answer_for_blastenv=getc;
         if($answer_for_blastenv=~/^y/i){
             &ask_for_ENV_vars($blast_env_var);
             print "\n# (Warn) You need to source your .bashrc or .cshrc file to activate the ENV var\n";
         }elsif($answer_for_blastenv=~/^n/i){
             print "\n# (i) I hope it is O.K., to run Blast without BLASTDB env var\n";
         }else{
             print "\n# (ERROR) You pressed none of y or n, I am exiting \n";
             exit;
         }
     }elsif($ENV{'BLASTDB'}=~/\S/){
         print "\n# (i) $blast_env_var is set to ", $ENV{'BLASTDB'}, "\n";
         $ENV_BLASTDB_var=$ENV{'BLASTDB'};
     }

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # getting options
     #________________________________________________________
     if($vars{'d'}){    $source_DB_file =$vars{'d'} };
     if(-s $vars{'i'} ){ $input_seq_file =$vars{'i'} };
     if($vars{'E'}){    $Eval_limit     =$vars{'E'} };
     if($vars{'j'}){    $iteration_limit=$vars{'j'} };
     if($vars{'e'}=~/(\S+)/ or $vars{'h'}=~/(\S+)/){ $step_evalue=$1 };
     if($char_opt=~/o/){ $over_write='o' }
     if($char_opt=~/D/){$make_msp_in_sub_dir_opt= 'D' } ## You'd better set this as default
     if($char_opt=~/L/){    $Lean_output        = 'L' }
     #print "\n# \$make_msp_in_sub_dir_opt $make_msp_in_sub_dir_opt\n"; exit;

     if(! $source_DB_file and ! $input_seq_file){
         print "\n# (Warn) You did not set d= and i= options for DB and Inquiry, I set automatically\n";
         $source_DB_file=${&get_largest_file(@file)};
         $input_seq_file=${&get_smallest_file(@file)};

         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~```
         # If one of the files is missing, ask if it is for self self search and
         #  make query=target, else die
         #_______________________________________________________________________
         $defined_all_ok=&check_if_defined($source_DB_file, $input_seq_file);
         unless($defined_all_ok){
             print "\n# (Q) $source_DB_file $input_seq_file Did you want to do self self search? ->(y/n) ";
             $answer_for_self_self=<STDIN>;
             if($answer_for_self_self =~/y/i){
                 if($input_file_name){  $sequence_DB=$input_file_name }
                 else{ $input_file_name=$sequence_DB };
                 close(STDIN);
             }else{
                 print "\n# (E) You seemed made a mistake, O.K., I will kill myself!\n\n";
                 print chr(7);  exit;
             }
         }
         print "\n# (i) \$source_DB_file is $source_DB_file, \$input_seq_file is $input_seq_file\n";
     }else{
         print "\n# (i) I have \$input_seq_file : $input_seq_file, \$source_DB_file : $source_DB_file\n";
         print "\n# (i) $0 uses \"gzip\" program, just warning :-) \n";
     }

     #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     # Check if formatdb has to be been run
     #_______________________________________________________________
     $blast_db_phr_file="$source_DB_file\.phr";
     $blast_db_pin_file="$source_DB_file\.pin";
     $blast_db_psq_file="$source_DB_file\.psq";
     unless(-s $blast_db_phr_file and -s $blast_db_pin_file and -s $blast_db_psq_file ){
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
         # Checking in the ENV BLASTDB for xxxx.phr, xxx.pin, xxxx.psq
         #_____________________________________________________
         unless( -s "$ENV_BLASTDB_var\/$blast_db_phr_file" and
                 -s "$ENV_BLASTDB_var\/$blast_db_pin_file" and
                 -s "$ENV_BLASTDB_var\/$blast_db_psq_file"){
              if(&find_program_in_path("formatdb")){
                  system("formatdb -t $source_DB_file -i $source_DB_file");
              }else{
                  print "\n# (E) I am trying to run formatdb to make blastdb, but can not find it\n";
                  print "\n# (E) Without running formatdb over seq DB, you can not run psi-blast\n";
                  exit;
              }
         }else{
              print "\n# (i) do_psi_blast_search: I found ENV BLASTDB ($ENV_BLASTDB_var) set";
         }
     }

     $base=${&get_base_names(\$source_DB_file)};

     %seq_input=%{&open_fasta_files(\$input_seq_file)};

     @seq_names= keys %seq_input;
     @seq_names= @{&scramble_array(\@seq_names)};

     print "\n# (i) Number of sequences to search is: ", scalar(@seq_names), " \n";

     for($i=0; $i< @seq_names; $i++){
         my(%msp_format_from_bla, $first_char);
         my $seq_name=$seq_names[$i];
         print "\n# (i) $seq_name is being processed, please wait! \n";
         my $temp_file_name="$seq_name\.fa";
         my $out_pbla_file="$seq_name\.pbla";
         my $out_msp_file="$seq_name\.msp";
         my $gzipped_pbla="$out_pbla_file\.gz";
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         # If D option is set
         #____________________________________________
         if($make_msp_in_sub_dir_opt){
             $first_char= substr("\U$seq_name", 0, $sub_dir_size);
             unless(-d $first_char){
                 mkdir ("$first_char", 0777);
                 print "\n# (i) do_psi_blast_search: You set \'d\' or \'D\' opt\n";
                 print "# (i) making subDIRs ($first_char) with $seq_name $sequence_DB to store MSP files\n";
             }
         }
         #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         # Let's make each fasta file for each seq to be used in searching
         #_____________________________________________________________________
         %single_seq=($seq_name, $seq_input{$seq_name});
         &write_fasta(\%single_seq, $temp_file_name ); ## e for writing each file

         if( !(-s $out_pbla_file) or $over_write){
             if($make_msp_in_sub_dir_opt){
                 $final_out_pbla_file="$first_char/$out_pbla_file";
             }else{
                 $final_out_pbla_file="$out_pbla_file";
             }
             print "\n# (i) Running Blast \n";
             if(&find_program_in_path("blastpgp")){
                 print "\n# (i) I am just running blastpgp which is in the PATH\n";
                 system("blastpgp -d $source_DB_file -i $temp_file_name -j $iteration_limit -e $Eval_limit -h $step_evalue -b 2500 -v 1000 -o $final_out_pbla_file ");
             }elsif(-s "/usr/local/bin/blastpgp"){
                 system("/usr/local/bin/blastpgp -d $source_DB_file -i $temp_file_name -j $iteration_limit -e $Eval_limit -h $step_evalue -b 2500 -v 1000 -o $final_out_pbla_file ");
             }elsif(-s "/local/bin/blastpgp"){
                 system("/local/bin/blastpgp -d $source_DB_file -i $temp_file_name -j $iteration_limit -e $Eval_limit -h $step_evalue -b 2500 -v 1000 -o $final_out_pbla_file");
             }

             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
             # Going into SUB directory if 'D' is set and Writing MSP
             #________________________________________________________
             if($make_msp_in_sub_dir_opt){
                 chdir("$first_char");
             }
             #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             # Converting PBLA to MSP file format, for each run (for geanfammer each conversion is O.K.
             #___________________________________________________________________________________________________
             %msp_format_from_bla=%{&convert_bla_to_msp($out_pbla_file, $Lean_output, 'A')};
             &compress_files_by_gzip($out_pbla_file) unless $Lean_output;
             print "\n# (i) Writing MSP file\n";
             &write_msp_files(\%msp_format_from_bla, $out_msp_file);
             push(@final_msp_out, $out_msp_file);
             if($make_msp_in_sub_dir_opt){  chdir('..');  }
         }else{
             print "\n# (Warn) It seems $out_pbla_file or $gzipped_pbla is already present, skipping psi run\n";
         }
         unlink($temp_file_name);
    }
    return(\@final_msp_out);
} # do_psi_blast_search


#______________________________________________________________________________
# Title     : die_if_file_not_present
# Usage     : &die_if_file_not_present($var, $file);
# Function  : checks if all the args are present
# Example   :
# Keywords  : die_unless_present, die_unless_file_present
# Options   :
# Author    : jong@salt2.med.harvard.edu
# Category  :
# Version   : 1.0
#------------------------------------------------------------------------------
sub die_if_file_not_present{
    my($i, $input_arg_count);
    for($i=0; $i< @_; $i++){
       $input_arg_count++;
       unless(-s $_[$i]){
           print "\n# (i)die_if_file_not_present: ${input_arg_count}th exist as a file in PWD\n";
       }
    }
}





#________________________________________________________________________
# Title     : cp
# Usage     : &cp( \$srcFile, \$dstFile); or  &cp( $srcFile, $dstFile);
#             or &cp(FILEHANDLE1, FILEHANDLE2),  or  &cp(FILEHANDLE1, $output)
# Function  : copies files fast, replacement of 'system("cp xxx xxxx"); '
# Example   : cp("cp.pl", *STDOUT);  # This will print cp.pl contents to your screen.
# Keywords  : copy files fast. cp_file, cp_files, copy_files, copy_file
# Options   :
# Argument  : 2 references of file name or  2 file names.
# Warning   : 27 times slower than 'cp' at prompt.  using system is 32 times slower
# Category  :
# Author    : Larry Wall, Jong
# Version   : 1.4
#--------------------------------------------------------------------
sub cp{
    croak(" Usage: cp ($file1, $file2)" ) unless(@_ >= 2);

    my($in)=$_[0];
    my($out)=$_[1];

    if( (ref($in) eq 'GLOB') || (ref($in) eq 'FileHandle')){
                    *IN_CP = *$in; }  ## for  \*STDOUT like input
    elsif( ref(\$in) eq 'GLOB'){
                    *IN_CP = $in;  }
    else{ open(IN_CP, "<$in") or die "Can't open output $in: $!\n";  }

    if( (ref($out) eq 'GLOB') || (ref($_[1]) eq 'FileHandle')){
                    *OUT_CP = *$out;   }
    elsif( ref(\$out) eq 'GLOB'){  print "\n2\n";
                    *OUT_CP = $out;  ## for  *STDOUT like input
    }else{  open(OUT_CP, ">$out") or die "Can't open output $out: $!\n";  }

    my ($access,$mod) = (stat IN_CP)[8,9];
    syswrite(OUT_CP, $buf, $len) while $len = sysread(IN_CP, $buf, 8192);
    close IN_CP;
    close OUT_CP;
    utime $access, $mod, $dstFile;
}


