#!/usr/bin/perl
# Last Update by /gn0/jong/Perl/update_subroutines.pl: Thu Jul 10 14:33:16 BST 1997
#_________________________________________________________________________
# Title     : convert_sso_to_msp.pl
# Usage     : convert_sso_to_msp.pl <1bia_1-63.sso>  or  <xxxxxx.clu>
# Function  : converts FASTA or SSEARCH -m 10  file output format to
#             matching sequence pair (MSP format)
#             WARN: This will concatenate all the output into ONE msp file
#               unless you set 'e' option at the prompt for each
#               xxx.msso to xxx.msp
#
# Example   : convert_sso_to_msp.pl steve_clufile.clu
# Keywords  : convert_sso_to_msp, convert_ssso_to_msp, convert_fsso_to_msp
# Options   :
#   e  for each msp file is made for each input xxxxxx.sso file
#   s  for single file out format (one big concatenated file will be made)
#   u= for upper E value limit (above this are discarded)
#   l= for lower E value limit (below this are discarded)
#   r  for adding ranges
#   $single_out_opt         =  s by s
#   $verbose                =  v by v -v     # to show some info
#   $uppercase_seq_name     =  U by U -U
#   $upper_expect_limit     =    by u=
#   $lower_expect_limit     =    by l=
#   $single_file_name       =    by s=
#   $make_each_msp          =  e by e -e
#   $over_write             =  o by o
#   $add_range              =  r by r -r
#   $add_range2             = r2 by r2 R
#   $new_format             =  n by n -n
#   $write_each_msp_to_disk =  w by w -w  # for writing each msp files to save memory
#   $eval_threshold=  by t=
#
# Returns   : Writes a file automatically
# Author    : Sarah Teichmann, Jong Park, sat@mrc-lmb.cam.ac.uk,jong@salt2.med.harvard.edu
# Version   : 2.0
#--------------------------------------------------------------
$|=1;
my (@sub_dir_heads, $s, $i, $j, $k, $clu);
local $eval_threshold=0.001;
my @sso_files=@{&parse_arguments(1)}; ## SSO files are (xx.ssso, xx.fsso, xx.sso)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Default parameter setting
#__________________________________
unless($upper_expect_limit=~/\S/){   $upper_expect_limit=0.081; }
unless($lower_expect_limit=~/\S/){   $lower_expect_limit=0; }
if($eval_threshold){ $upper_expect_limit=$eval_threshold; }

print "\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ";
print "\n# Input files are                 opt: @sso_files ";
print "\n# \$add_range                      opt: $add_range";
print "\n# \$upper_expect_limit             opt: $upper_expect_limit";
print "\n# \$uppercase_seq_name             opt: $uppercase_seq_name";
print "\n# \$single_out_opt                 opt: $single_out_opt";
print "\n# \$verbose                        opt: $verbose";
print "\n# \$over_write                     opt: $over_write";
print "\n# \$write_each_msp_to_disk         opt: $write_each_msp_to_disk";
print "\n#_______________________________________________________________\n";

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sso_files has either xxxxx.sso  or xxxxx.gz
#______________________________________________________
for($i=0; $i < @sso_files; $i++){
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # remove .clu file and save it in $clu
   #________________________________________
   if($sso_files[$i]=~/\.s?clus?/){	   $clu=$sso_files[$i];
      splice(@sso_files, $i, 1);	   $i--;
   }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     If clu file is given, let's open it
#_______________________________________________________________
if(defined($clu)){
    &clu_to_sso_to_msp(\$clu);
}else{
    print "\n# No clu file is given. I open individual fasta or ssearch output\n";
    if($make_each_msp=~/e/){
	 print "\n# You put \"e\" option for each file conversion\n";
	 for($i=0; $i< @sso_files; $i ++){
             &convert_sso_to_msp($sso_files[$i],
                 $single_out_opt,
                 "s=$single_file_name",
		 $add_range, $add_range2,
		 "u=$upper_expect_limit",
                 "t=$eval_threshold",
	         "l=$lower_expect_limit",
	         $new_format);
         }
    }else{
	 print "\n# $0 : You did not put e opt, so I will concatenate all the sso contents to one MSP\n";
	 &convert_sso_to_msp(@sso_files,
	    $single_out_opt, "s=$single_file_name",
	    $add_range,
	    $add_range2,
	    "u=$upper_expect_limit",
	    "l=$lower_expect_limit",
            "t=$eval_threshold",
	    $new_format);
    }
}


print "\n# Sarah! Please remember that $0 used default E value of \"$upper_expect_limit\" for upper limit";
print "\n#        Also,  remember that $0 used default E value of \"$lower_expect_limit\" for lower limit\n";


#___________________________________________________________________
# Title     : convert_sso_to_msp
# Usage     : &convert_sso_to_msp(@ARGV, $single_out_opt);
# Function  : This takes sso file(s) and produces MSP file. It
#             concatenate sso file contents when more than one
#             sso file is given.
# Example   : &convert_sso_to_msp(@ARGV, 'OUT.msp', $single_out_opt);
# Warning   : This capitalize all the input file names when
#              producing xxxxx.msp. xxxxx.sso -> XXXX.sso
# Keywords  :
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
#   $eval_threshold =  by e=
#
# Returns   : the file names created (xxxx.msp, yyyy.msp,,,,)
# Argument  :
# Category  :
# Version   : 2.7
#--------------------------------------------------------------
sub convert_sso_to_msp{
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
   if($vars{'u'}=~/(\.?\d+)/){ $upper_expect_limit = $vars{'u'} };
   if($vars{'l'}=~/(\.?\d+)/){ $lower_expect_limit = $vars{'l'} };
   if($vars{'t'}=~/(\d+)/){ $Score_thresh  = $vars{'t'} };
   if($vars{'e'}=~/(\.?\d+)/){ $upper_expect_limit = $vars{'e'} };
   if($vars{'m'}=~/(\d+)/){ $margin  = $vars{'m'} };
   if($vars{'s'}=~/\S/){ $single_file_name  = $vars{'s'} };
   if($char_opt=~/r2/){  $add_range='r'; $add_range2='r2' }
   if($char_opt=~/r/){   $add_range = 'r' }
   if($char_opt=~/c/){   $create_sso = 'c' }
   if($char_opt=~/s/){   $single_out_opt='s' }
   if($char_opt=~/n/){   $new_format='n' }
   print "\n# (i) File given to convert_sso_to_msp is \"@file\", Normally xxx.sso file names\n";

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
              print "\n# convert_sso_to_msp: File matched  xxxx.sso  format \n";
          }elsif($file[$i]=~/^(\S+)\.out$/){
              $big_out_msp1="\U$1"."\.msp";
              $big_out_msp2="\U$1"."_all".".msp";
              print "\n# convert_sso_to_msp: File matched  xxxx.out  format \n";
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
       print "\n# convert_sso_to_msp: You did not define the big MSP file out format, so $big_out_msp1 \n";
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``~~~~
   #  (1) When File was given(normal situ) to this sub routine
   #_____________________________________________________________
   if(@file == 1){   ## ONE single file input??
	  print "# one file @file is given, OUT will be: $big_out_msp1 \n";
          @sso=@{&open_sso_files(@file, $add_range, $add_range2,
	          "u=$upper_expect_limit",
		  "l=$lower_expect_limit",
		  "m=$margin", $over_write,
		  $new_format,
		  "s=$big_out_msp")};
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # Write MSP and finish off
          #______________________________________________
	  push(@final_out, &write_msp_files(@sso, $big_out_msp1,
	        $single_out_opt, $add_range) );

   }elsif(@file > 1){ ## MOre than 1 file input??
       @sso=@{&open_sso_files(@file, $add_range, $add_range2,
             "l=$lower_expect_limit",
             "u=$upper_expect_limit",
             "m=$margin", $over_write,
             $new_format)};
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # Write MSP and finish off (multiple @file input)
       #______________________________________________
       push(@final_out, @{&write_msp_files(@sso, $big_out_msp2,
             $single_out_opt, $add_range)} ); ## concatenates all the hash ref to one
   }

   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  (2) When NO File but ARRAY is given
   #      Here, you can have SSO files created
   #__________________________________________
   elsif(@array >=1){
       print "\n# In convert_sso_to_msp, \@array is given rather than \@file";
       @sso=@{&open_sso_files(@array, "u=$upper_expect_limit", $add_range2,
                  "l=$lower_expect_limit", $add_range, $create_sso,
                  "m=$margin", $new_format, $over_write)};
       push(@final_out, @{&write_msp_files(@sso, $big_out_msp,
                                               $single_out_opt, $add_range)} );
   }
   return(\@final_out);
}



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
#              convert_clu_to_sso_to_msp
# Options   : USE, convert_clu_to_sso_to_msp, this is obsolute now
# Category  :
# Version   : 1.7
#--------------------------------------------------------------------------------
sub clu_to_sso_to_msp{
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




#________________________________________________________________________
# Title     : read_dir_names_only
# Usage     : @all_dirs_list = @{&read_dir_names_only(\$absolute_path_dir_name, ....)};
# Function  : read any dir names and and then put in array.
# Example   :
# Warning   : This does not report '.', '..'
#             Only file names are reported. Compare with &read_any_dir
# Class     :
# Keywords  :
# Options   :
# Package   :
# Reference :
# Returns   : one ref. of array.
# Tips      :
# Argument  : takes one or more scaler references. ('.', \$path, $path, ... )
# Todo      :
# Author    : A Biomatic
# Version   : 3.1
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub read_dir_names_only{
  my($in_dir, $i,$k, @possible_dirs,
	  @final_files, $full_dir, $pwd, $path,@read_files);
  $pwd=`pwd`; chomp($pwd); $full_dir=1;
  for($k=0; $k < @_; $k++){
	 if   ( ($_[$k] eq '.') || !(defined($_[$k]))){  $in_dir=$pwd;  }
	 elsif(!(ref($_[$k]))){   $in_dir=$_[$k];   }
	 elsif(ref($_[$k])){      $in_dir =${$_[$k]};    }
	 if($in_dir =~ /^([\w\-\.]+)$/){  $in_dir="$pwd\/$in_dir"; $full_dir = 0; }
	 else{ $full_dir =1; }
	 ##########  Main READING PART ##########
	 opendir(DIR1,"$in_dir");
	 @read_files = readdir(DIR1);
	 for($i=0; $i < @read_files; $i ++){
		$read_files[$i]="$in_dir\/$read_files[$i]";
		if( ($read_files[$i] !~ /\/\.\.?$/) && ( -d $read_files[$i]) ){
		  $read_files[$i]=~s/\.\///; ## removing ./ in front of dirs (in bash)
		  push(@final_files, "$read_files[$i]");
		}
	 }
  }
  return([sort @final_files]);
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
#             generalized debug var is added for more verbose printouts.
#--------------------------------------------------------------------
sub set_debug_option{
  my($j, $i, $level );
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
# Version   : 1.7
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
		 $extension=$1;
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
		 if($extension){
			 if($ARGV[$f]=~/\S\.$extension/){
				 push(@input_files, $ARGV[$f] );
			 }else{ next }
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
  \%Final_out;
}               ##<<--- ENd of the sub read_head_box


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
# Keywords  :
# Options   : _  for debugging.
#             #  for debugging.
#             b  for to get just names ($simple_clu_reading)
#             r  for adding ranges in the names
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
# Version   : 2.8
#--------------------------------------------------------------------
sub assign_options_to_variables{
    my($i, %vars, $j, $op, $z, $n, $symb, $value, $var, %val, @val, $ARG_REG,
           $option_table_example, @input_options, $first_border_and_title, $sym, @arg);

    #""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #      Defining small variables for option table reading
    #________________________________________________________________________
    my($g)='gets';                my($if)='if';
    my($is)='is';                 my(@input_files);
    my($o)='or';   my(@arguments) = sort @ARGV;

    #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #  Assigning global arguments(@num_opt, %vars) variables
    #_______________________________________________________________
    for($i=0; $i< @arguments; $i++){
         print "\n# \"$arguments[$i]\"";
         if(($arguments[$i]=~/^(\-?\d+[\.\d+]?)$/)&&   ### it mustn't be a file
             ( !(-f $arguments[$i]) ) ){                ### getting NUM opt
             push(@num_opt, $1);
         }elsif( $arguments[$i]=~/^(\S+)=(\S+)$/){
             $vars{$1}=$2;
             print "\n# (1) assign_options_to_variables: Your variable option for \"$1\" is \"$2\"\n";
         }elsif( $arguments[$i]=~/^(\S+)=[\'\"]([ \S]+)[\'\"]$/){ # just in case Shell fails
             $vars{$1}="$2";
         }elsif( $arguments[$i]=~/^(\S+)=([ \S]+)$/){
             $vars{$1}="$2";
         }
    }


    #""""""""""""""""""""""""""""""""""""""""""""""""""
    #   Some DEFAULT $debug variables for debugging purposes
    #__________________________________________________________
    &set_debug_option;

    #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    #   The main processing of self
    #_______________________________________________________________
    open(SELF, "$0");    ## opens the program you ran to get the options table.
    while(<SELF>){

       if( $first_border_and_title > 5 ){  ## This is to make it read only the first headbox.
            last;                            #  $first_border_and_title is an incremental counter.
       }elsif( /^ *#[_\*\-]{15,}$/ and /^ *# *[Tt][itle]*[ :]*/ ){
            $first_border_and_title++;
            print __LINE__, "# assign_options_to_variables : Title line found\n" if $debug eq 1;
       }elsif(/^ {0,5}#[\t ]{1,50}([\$\%\@].+)$/){
            $op = $1;  ## $op is for the whole input option line which has $xxxx, @xxx, %xxxx format
            $op =~ s/^(\W\w+.+)(\s+\#.*)$/$1/;  ## This is removing any comments in the line.
            #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            ## matching the following line input format.
            ## $av_sc_segment     becomes    a  by  a  # To smooth the SC rates. Gets the averages of
            ## $ARG_REG is for arguments regular expression variable.
            ##  This reg. exp. matches = 'a or A or E or e' part
            ##  which represents alternative prompt arguments possibilities. \=$b$g$is$e$set
            #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # These are optional argument flag like (-a as well as a)
            #__________________________________________________________
            $ARG_REG ='(\S*) *[or=\,]* *(\S*) *[=or\,]* *(\S*) *[=or\,]* *(\S*) *[=or\,]* *(\S*)';

            if(#$op=~/^([\$\@\%])([\w\-]+)[\t ]{0,25}[=|$g|$is]([\t ]+)[bB]y +$ARG_REG/
               #or
               $op=~/^([\$\@\%])([\w\-]+)[\t ]{0,25}[=|$g|$is][\t ]*[\$\@\%]?(\S*) *[bB]y +$ARG_REG/){
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
                      print "\n# Prompt:  \"$arguments[$z]\"\n" if $debug;
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
                               ( $arguments[$z] =~ /^${1}= *([ \S]+)$/) &&  ## NOTE  '[ \S]+' not \S+
                               ( $sym eq '$') ){
                                 ${"$var"}="$1";
                                 print "\n# I am assigning $1 to $var\n";
                                 if($debug eq 1){ print __LINE__,"\$${var} is set to \"$1\"\n";  }
                           }
                      }
                 }
             }
          }
     }
}



#______________________________________________________________________________________
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
         print "\n# (i) \@sso has ", scalar(@sso), " lines. \n"  if $verbose;
         if(@sso > 3000){ # if @sso is very big, I remove the useless contents
             print "\n# (i) open_sso_files: size of \@sso for $file[$i] exceeds 3000 lines, ", scalar(@sso), " !!! \n";
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
                  print "\n# (i) openning text File format xxxx.[fms]sso $file[$i]";
                  open(SSO, "$file[$i]") or die "\n# (E) open_sso_files: Failed to open $file[$i]\n";
                  my @sso=<SSO>;
                  if(@sso < 30){  @sso=`zcat $file[$i]`; }      # if zcat fails to produce output use gunzip -c
                  if(@sso > 3000){ # if @sso is very big, I remove the useless contents
                      print "\n# (i) open_sso_files: size of \@sso is for $file[$i] exceeds 3000 lines, ",
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
    print "\n# (i) \$attach_range_in_names2 is $attach_range_in_names2\n" if $attach_range_in_names2;

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
            print "\n (i) target seq len is  $query_seq_len \n";

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
            print "\n# (i) read_machine_unreadable_sso_lines : Found MATCHed seq $match_seq\n" if $verbose;
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
       print "\n# (i) $query_seq_name\.sso  is created";
       close SSO2;
    }
    push(@out_refs, \%match);
    return(\@out_refs);
}# end of for $i


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
# Version   : 1.4
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
             if($expect=~/^$/){ $expect='0.0'; }
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
                         sprintf("%-5s %-8s %-6s %-4s %-5s %-30s %-4s %-5s %s\n",
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
                    $match{$match_seq}=sprintf("%-5s %-8s %-6s %-4s %-5s %-30s %-4s %-5s %s\n",
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

