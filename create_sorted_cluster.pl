#!/usr/bin/perl

#Written by Sarah A. Teichmann on December 17th 1996.

#Use: cut-and-pastes all .clu files in a directory into a single sorted file called "sorted_cluster_file".

$|=1;

&create_sorted_cluster;


#______________________________________________________________________________
# Title     : create_sorted_cluster
# Usage     : &create_sorted_cluster
# Function  : to make a "sorted_cluster_file" from the .clu files in a directory
# Example   :
# Keywords  : make_cluster_file, sort_clu_files
# Options   :
# Author    : Sarah A. Teichmann
# Date      : 19th September 1997
# Version   : 1.5
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

    for ($i=0; $i < @clufiles; $i++) {
         open (CLU_FILE, "<$clufiles[$i]") or die "\n# $0: create_sorted_cluster: error opening $clufiles[$i]";
         my $cluster_size;
         while(<CLU_FILE>){
             $line=$_;
             if( /^ *Cluster +size +(\d+)/i){
                 $cluster_size=$1;
                 $hash{$cluster_size} .=$line;
             }elsif (/^ *Cluster +[number]* *\d+/i) {
                 $hash{$cluster_size} .=$line;
             }elsif (/^ *\d+ +\d+ +\S+/) {
                 $hash{$cluster_size} .=$line;
             }
        }
   }

   @sorted_by_size= sort { $a <=> $b } keys %hash;
   my $good_cluster_file="sorted_cluster_file\.gclu";
   open(GOODCLUS, ">$good_cluster_file") || die "\n# $0 create_sorted_cluster: I can not open $good_cluster_file\n";
   for($i=0; $i< @sorted_by_size; $i++){
       print GOODCLUS $hash{$sorted_by_size[$i]};
   }
   return(\$good_cluster_file);
}



#_________________________________________________________________________
# Title     : read_file_names_only
# Usage     : @all_files=@{&read_file_names_only(<dir>, [extension])};
# Function  : read any file names and REMOVES the '.', '..' and dir entries.
#             And then put in array.  This checks if anything is a real file.
#             You can use 'txt' as well as '.txt' as extension
#             You can put multiple file extension (txt, doc, ....)
#               and multiple dir path (/usr/Perl, /usr/local/Perl....)
#               It will fetch all files wanted in all the direc specified
#
# Example   : @all_files=@{&read_file_names_only(\$abs_path_dir_name, ..)};
#             @all_files=@{&read_file_names_only(\$dir1, '.pl', '.txt')};
#             @all_files=@{&read_file_names_only(\$dir1, '.', \$dir2, \$dir3, 'e=pl')};
#
# Warning   : This does not report '.', '..'
#             Only file names are reported. Compare with &read_any_dir
#             extension size should be less than 15 char.
#             It sorts the results!
# Class     :
# Keywords  : filename only, filename_only, read_files_only, read files
#             get_file_names_only, get_files_only, read_files_only
# Options   : "extension name". If you put , 'pl' as an option, it will show
#             files only with '.pl' extension.
#             '-p'  for path also included resulting in '/path/path/file.ext'
#                  rather than 'file.ext' in output @array
#  e='xxx'  for extention xxx
#  '.pl'    for files extended by '.pl'
#  'pl'     for files extended by 'pl', same as above
#
# Package   : File
# Reference :
# Returns   : one ref. of array.
# Tips      :
# Argument  : takes one or more scaler references. ('.', \$path, $path, ... )
# Todo      :
# Author    :  A Biomatic
# Version   : 1.9
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub read_file_names_only{
  my($in_dir, $i,$k, @final_files, @possible_dirs, $ext, @extensions,
          $path_include, @in, @in_dir, $pwd, $extension_given, @read_files);
  $pwd=`pwd`; chomp($pwd);
  $in_dir=$pwd;   @in=@_;
  for($k=0; $k < @in; $k++){
         if   ( $in[$k] eq '.'){ push(@in_dir,$pwd); splice(@in, $k, 1);  $k--; next }
         if( !(ref($in[$k]))){
                if( -d $in[$k]){ push(@in_dir, $in[$k]); splice(@in, $k, 1);    $k--;
                }elsif(!(-f $in[$k]) and $in[$k] =~ /^\-p$/ ){$path_include=1; splice(@in, $k, 1); $k--;
                }elsif(!(-f $in[$k]) and $in[$k] =~ /^\.?(\S+)/){       $extension_given =1; push(@extensions, $1);     splice(@in, $k, 1);       $k--;
                }elsif(!(-f $in[$k]) and $in[$k] =~ /^e=\.?(\S+)/){ $extension_given =1; push(@extensions, $1); splice(@in, $k, 1);$k--;  }
         }elsif(ref($in[$k])){
                if( -d ${$in[$k]}){     push(@in_dir,${$in[$k]});  splice(@in, $k, 1);  $k--;
                }elsif(!(-f ${$in[$k]}) && (${$in[$k]} =~ /^\.?(\S+)/ ) ){
                   $extension_given = 1;  push(@extensions, $1);  splice(@in, $k, 1);  $k--;
                }elsif(!(-f ${$in[$k]}) && (${$in[$k]} =~ /^e=(\S+)/ ) ){
                   $extension_given =1; push(@extensions, $1);  splice(@in, $k, 1);  $k--;
                }
         }
  }
#print "@extensions\n";
  ##########  Main READING PART ##########
  for($k=0; $k< @in_dir; $k++){
         opendir(DIR1,"$in_dir[$k]");
         @read_files = readdir(DIR1);
         for($i=0; $i < @read_files; $i ++){
                if( -f "$in_dir[$k]\/$read_files[$i]" ){
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
                  }else{
                          push(@final_files, $read_files[$i]);
                  }
                }
         }
  }
  sort @final_files;
  return(\@final_files);
}

#________________________________________________________________________
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
# Version   : 1.2
# Used in   :
# Enclosed  :
#--------------------------------------------------------------------
sub remove_dup_in_array{
  my($i, @out_ref, @nondup,%duplicate, @orig, @out_ref);
  for($i=0; $i<@_; $i++){
	  undef(%duplicate);
	  if(ref($_[$i]) eq 'ARRAY'){    @orig = @{$_[$i]};    }
	  else{ @orig=@_ }
	  @nondup = grep { ! $duplicate{$_}++ } @orig;
	  push(@out_ref, \@nondup);  }
  if(@out_ref ==1){ return($out_ref[0]);}
  elsif(@out_ref >1){  return(@out_ref);}
}


sub numerically { $a <=> $b; }
