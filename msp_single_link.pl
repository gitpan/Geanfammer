#!/usr/bin/perl
#__________________________________________________________________________
# Title     : msp_single_link.pl
# Usage     : msp_single_link.pl <e-val>
# Function  : to single link all the genes in msp files in a given directory and its subdirectories
#             at a given e-value threshhold and make a sorted file called "single_linkage.clus".
# Example   :
# Keywords  :
# Options   :
# Returns   : Sarah A. Teichmann on the 26th August 1997.
# Argument  :
# Version   : 1.1
#----------------------------------------------------------------------------

$|=1;

my $usage="msp_single_link.pl  <e-val> ";

die "\nIncorrect number of arguments.\n\nUsage: $usage\n\ndied" unless ($#ARGV==0);

my $E_cut_main=$ARGV[0];

my @msp_files_main=@{&get_all_msp_files};

my %hash_main=%{&msp_single_link_hash(\@msp_files_main, $E_cut_main)};

&print_clusfile_from_hash(\%hash_main);



#________________Subroutine 1: fetch all msp/msp.gz files and put them into an array________________#

#________________________________________________________________________________
# Title     : get_all_msp_files
# Usage     : @msp_files=@{&get_all_msp_files};
# Function  : puts the names of all msp or msp.gz files in the directory and its subdirectories into an array
# Example   :
# Keywords  : read_msp_files, make_msp_file_array, get_all_msp_files_in_dir
# Date      : 19th September 1997
# Options   :
# Author    : Sarah A. Teichmann, jong@salt2.med.harvard.edu
# Version   : 1.4
#--------------------------------------------------------------------------------
sub get_all_msp_files {
   my (@msp_files_main_final);
   my @msp_files_main=@{&read_file_names_only('.msp','.msp.gz')};
   for(@msp_files_main){
      if($_=~/cluster\.msp/){ next ## xxxx_cluster.msp is the processed file
      }else{ push(@msp_files_main_final, $_) }
   }
   my @dirs=@{&read_dir_names_only('n', '.')};
   my ($i, $j, @msp_files, $file_in_dir, $dir_file);
   print "\n# (i) Trying to read in msp files in subdir\n";
   for ($i=0; $i<@dirs; $i++){
       my $dir=$dirs[$i];
       unless( -d $dir){
          next;
       }
       if( -d $dir){  chdir($dir);  }

       @msp_files=@{&read_file_names_only('.','.msp','.msp.gz')};
       for ($j=0; $j<@msp_files; $j++){
          if($msp_files[$j]=~/cluster\.msp/){ next } ## skipping already existing  20-1_cluster.msp like files
          $file_in_dir=$msp_files[$j];
          $dir_file="$dir"."/"."$file_in_dir";
          push(@msp_files_dirs, $dir_file);
          next;
       }
       chdir('..');
       next;
   }

   @msp_files=sort (@msp_files_main_final, @msp_files_dirs);
   return (\@msp_files);
}



#_________Subroutine 2: parse the msp and msp.gz files and put the genes into a hash with cluster values ________#


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
# Version   : 1.5
#--------------------------------------------------------------------------------
sub msp_single_link_hash {
    my (@msp_files, $i, $j, $k, $e_val, $gene_1, $gene_2,
        @mspcont, $gene_1, $gene_2, $E_cut, %hash, $array);

    $E_cut=0.001; ## Default Eval cut

    if( @_==2 and ref($_[0]) eq 'ARRAY'){
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
            if ($mspcont[$j]=~/^\S+ +(\S+) +\S+ +\d+ +\d+ +(\S+) +\d+ +\d+ +(\S+)/){
                $e_val=$1;
                unless($e_val<=$E_cut){next;}
                $gene_1=$2;
                $gene_2=$3;
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
                # Removing REGION information as Sarah's can not handle regions
                #_______________________________________________________________
                if($gene_1=~/^(\S+)_\d+\-\d+/){ $gene_1=$1 }
                if($gene_2=~/^(\S+)_\d+\-\d+/){ $gene_2=$1 }

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
        next;
    }
    return (\%hash);
} ## msp_single_link_hash



#___________________Subroutine 3: print out the hash with the subroutine values_________________________#

#________________________________________________________________________________
# Title     :print_clusfile_from_hash
# Usage     :&print_clusfile_from_hash(\%hash)
# Function  :To print out a file in cluster file format from an input hash containing the genes as keys and the cluster number as values.
# Example   :
# Keywords  :print_single_linkage_cluster, print_cluster_file
# Options   :
# Version   : 1.0
# Author    : Sarah A. Teichmann
#--------------------------------------------------------------------------------
sub print_clusfile_from_hash {
    my (%hash, $single_linkage_cluster);
    $single_linkage_cluster="single_linkage.clus";
open(SING, ">$single_linkage_cluster");
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
@clus_sizes=sort bynumber (@clus_sizes);

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

}




#______________________________________________________________
# Title     :bynumber
# Usage     :sort bynumber <arrayname>
# Function  :sorts arrays by number
# Example   :sort bynumber (@array)
# Keywords  :sort_by_number, sort_array_by_number
# Options   :
# Returns   :array sorted by number
# Author    : A Scientist
# Version   : 1.0
#--------------------------------------------------------------
sub bynumber {
   $a <=> $b;
}



#______________________________________________________________
# Title     : remove_dup_in_array
# Usage     : @out2 = @{&remove_dup_in_array(\@input1, \@input2,,,,)};
#             @out1 = @{&remove_dup_in_array(\@input1 )};
# Function  : removes duplicate entries in an array.
# Example   : (1,1,1,1,3,3,3,3,4,4,4,3,3,4,4);  --> (1,3,4);
# Warning   :
# Class     :
# Keywords  : merge array elements, remove_repeting_elements,
#             remove_same_array_elements
# Options   :
# Package   :
# Reference :
# Returns   : one or more references.
# Tips      :
# Argument  : one or more refs for arrays or one array.
# Todo      :
# Author    :
# Version   : 1.3
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub remove_dup_in_array{
  my($i, $sort_opt, @out_ref, @nondup,%duplicate, @orig, @out_ref);
  my @in=@_;
  for($i=0; $i<@in; $i++){
		 if($in[$i] eq 's'){
				$sort_opt=1;
				splice(@in, $i, 1);
				$i--;
		 }elsif( (ref($in[$i]) eq 'SCALAR')&&(${$in[$i]} eq 's') ){
				$sort_opt=1;
				splice(@in, $i, 1);
				$i--;
		 }
  }
  for($i=0; $i<@in; $i++){
		  undef(%duplicate);
		  if(ref($in[$i]) eq 'ARRAY'){    @orig = @{$in[$i]};    }
		  else{ @orig=@in }
		  @nondup = grep { ! $duplicate{$_}++ } @orig;    ## NOTE -> $_
		  if($sort_opt==1){ @nondup= sort @nondup }
		  push(@out_ref, \@nondup);  }
  if(@out_ref ==1){ return($out_ref[0]);}
  elsif(@out_ref >1){  return(@out_ref);}
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
# Category  :
# Version   : 2.9
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

