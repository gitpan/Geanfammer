#!/usr/bin/perl
# Last Update by /gn0/jong/Perl/update_subroutines.pl: Mon Apr 21 12:56:52 BST 1997
#______________________________________________________________
# Title     : delbut.pl
# Usage     : delbut *.zip  (delete files except xxxx.zip)
# Function  :
# Example   :
# Warning   :
# Keywords  :
# Options   : _  for debugging.
~~~~~~~

# Returns   :
# Argument  :
# Category  :
# Version   : 1.2
#--------------------------------------------------------------


&delbut(@ARGV);



#______________________________________________________________
# Title     : delbut
# Usage     : delbut *.zip  (delete files except xxxx.zip)
# Function  :
# Example   :
# Warning   :
# Keywords  :
# Options   : _  for debugging.
#             #  for debugging.
# Returns   :
# Argument  :
# Category  :
# Version   : 1.3
#--------------------------------------------------------------
sub delbut{
		my ($remove_subdir, $i);
		@save_files=@{$_[0]} || @_;
		$remove_subdir=${$_[1]} || $_;
		for(@save_files){
			 if($_=~/^s *$/ and !(-e $_)){ $remove_subdir='s'; next }
			 unless(-e $_){
				 print "\n\n \"$_\" does not exist, so nothing is deleted\n\n";
				 print chr(7);
				 exit;
			 }
		}
		my @files=@{&read_dir_and_file_names_only('.')};
		my @del_files=@{&subtract_array(\@files, \@save_files)};

		for($i=0; $i< @del_files; $i++){
			 if(-d $del_files[$i]){
			    if( $remove_subdir=~/s/i){
						 system("rm -fr $del_files[$i]");
					}else{
					   print "\n# subdir $del_files[$i] has not been deleted\n";
					}
			 }else{
			    unlink($del_files[$i]);
			 }
		}
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
# Title     : subtract_array
# Usage     : @subs = @{&subtract_array(\@array1, \@array2)};
# Function  : removes any occurances of certain elem. of the first
#             input array with second input array.
# Example   : Following will produce (A K C);
#		@array1= qw( A B K B B C);
#  		@array2= qw( B E D);
#  		@subs = @{&subtract_array(\@array1, \@array2)};
# Keywords  : array_subtract, substract_array, ary1_minus_ary2
# Options   :
# Returns   :
# Argument  :
# Version   : 1.4
#--------------------------------------------------------------------
sub subtract_array{
	my(@first)=@{$_[0]};
	my(@second)=@{$_[1]};
	my %counter;
	grep($counter{$_}++, @second );
	return ( [grep(!$counter{$_}, @first)] );
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
# Title     : read_dir_and_file_names_only
# Usage     : @all_dir_and_file=@{&read_dir_and_file_names_only(<dir>, [extension])};
# Function  : read any dir/file names and REMOVES the '.', '..' and dir entries.
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
# Keywords  : filename only, filename_only, read_files_only, read files
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
# Version   : 1.0
#--------------------------------------------------------------------
sub read_dir_and_file_names_only{
		my($in_dir, $i, $j, $x, $k, $dir, @final_files, @possible_dirs, $sort_opt, $ext, @extensions,
			 $path_include, @in, $glob_given, @files_globed, @in_dir, $pwd, $extension_given,
			 %target_file_names, @target_file_names, @read_files);
		$pwd=`pwd`; chomp($pwd);
		$in_dir=$pwd;
		@in=@_;

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#  Directory entry and opts detection
		#_________________________________________
		for($k=0; $k < @in; $k++){
			 if   ( $in[$k] eq '.'){ push(@in_dir,$pwd); splice(@in, $k, 1);  $k--; next }
			 if( !(ref($in[$k]))){
					if($in[$k]=~/D=(\S+)/){
								print "\n# read_dir_and_file_names_only : $1 is used as input dir ";
								push(@in_dir, $1); splice(@in, $k, 1);    $k--; next;  }
					if( -d "$in[$k]" ){
							print "\n# read_dir_and_file_names_only: $in[$k] is a dir";
							if($in[$k]=~/\/\S+$/){
									$path_include=1;  ## If the input dir has '/', I assume path should be added to out file names
									print "\n# read_dir_and_file_names_only: \$path_include is set to 1";
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
					}elsif(!(-f $in[$k]) and $in[$k] =~ /^\-s *$/   ){
								$sort_opt=1; splice(@in, $k, 1); $k--;
					}else{
								print "\n# (W) read_dir_and_file_names_only: $in[$k] not a file, nor dir, a file extnsion?\n";
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
		 print "\n# read_dir_and_file_names_only: Final input directories are : @in_dir";
		 print "\n# read_dir_and_file_names_only: going to \'File name and extension detection\' stage with \@in";
	}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#  File name and extension detection
	#_________________________________________
	for $dir (@in_dir){
			chdir($dir);
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
									 print "\n# read_dir_and_file_names_only: pushing $1 as an extension" if $verbose;
									 $extension_given =1; push(@extensions, $1);
									 splice(@in, $k, 1); $k--;
							 }elsif(!(-f $in[$k]) and $in[$k] =~/^([^\-]{0,8})$/){  ## extension name can not be larger than 8 chars
									 print "\n# read_dir_and_file_names_only: pushing $1 as an extension" if $verbose;
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
		 print "\n# read_dir_and_file_names_only: You used glob for file name, but without extension name\n" if $verbose;
	   return(\@final_files);
	}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	#  Main READING PART
	#_________________________________________________
	print "\n# read_dir_and_file_names_only: \@in_dir is  @in_dir\n";
	for($k=0; $k< @in_dir; $k++){
		 chdir($in_dir[$k]) or die "\n# read_dir_and_file_names_only: could not get into $in_dir[$k]\n";
	   opendir(DIR1, ".");
		 @read_files = readdir(DIR1);
	   if(@read_files < 1){ print "\n# read_dir_and_file_names_only: ERROR??, \@read_files is empty\n\n\n"; }
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
						if( -f "$read_files[$i]" or -d "$read_files[$i]"){ ##
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
											#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
											# The Actual collecting names in HERE, skipping '.' and '..'
											#________________________________________________________
											unless($read_files[$i] eq '.' or $read_files[$i] eq '..'){
											   push(@final_files, $read_files[$i]);
											}
								 }
						}
				}
	 }
	 chdir($pwd);
	 }
	 @final_files=sort @final_files if $sort_opt == 1;
	 return(\@final_files);
}


