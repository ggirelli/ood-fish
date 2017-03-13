#!/usr/bin/env bash

# ------------------------------------------------------------------------------
#
# Date: 170307
# Author: Gabriele Girelli
# Description:	run the recap_blat.R script for every split PSL file in the
# 				specified psl (-i) directory.
#
# ------------------------------------------------------------------------------



# PARAM ========================================================================

# Help string
helps="
 usage:	./recap_split.sh -i psl_dir -o out_dir -s seq_dir
 	[-f seq_file][-t threads][-p path]
 
 Description:
  Run the recap_blat.R script for every split PSL file
  in the specified psl (-i) directory.

 Mandatory arguments:
  -i psl_dir	Directory containing the split PSL files.
  -s seq_dir	Directory containing the seq_file (-f).
  -o out_dir	Output directory.

 Optional arguments:
  -f seq_file	Name of the sequence file in seq_dir (-s).
  -t threads	Number of threads for parallelization.
  -p path	Path to recap_blat.R script. Default './recap_blat.R'
"

# Default values
seq_file="bc.noHeader.fa"
script_path="./recap_blat.R"
threads=1

# Parse options
while getopts hi:o:s:f:t:p: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 1
		;;
		i)
			if [ -d "$OPTARG" ]; then
				psl_dir=$OPTARG
			else
				msg="Invalid -i option, folder not found.\nFolder: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
		o)
			out_dir=$OPTARG
			if [ ! -d "$OPTARG" ]; then
				msg="Output folder not found, creating it."
				mkdir -p $out_dir
			fi
		;;
		s)
			if [ -d "$OPTARG" ]; then
				seq_dir=$OPTARG
			else
				msg="Invalid -s option, folder not found.\nFolder: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
		f)
			seq_file=$OPTARG
		;;
		t)
			if [ 0 -ge "$OPTARG" ]; then
				echo -e "Enforcing a minimum of 1 thread.\n"
			else
				threads=$OPTARG
			fi
		;;
		p)
			if [ -e $OPTARG ]; then
				script_path=$OPTARG
			else
				msg="Invalid -p option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$psl_dir" ]; then
	echo -e "$helps\nMissing mandatory -i option.\n"
	exit 0
fi
if [ -z "$out_dir" ]; then
	echo -e "$helps\nMissing mandatory -o option.\n"
	exit 0
fi
if [ -z "$seq_dir" ]; then
	echo -e "$helps\nMissing mandatory -s option.\n"
	exit 0
fi

# Print options
echo -e "\nPSL input directory: $psl_dir"
echo -e "Sequence input directory: $seq_dir"
echo -e "Sequence file: $seq_file"
echo -e "Threads: $threads"
echo -e "Output directory: $out_dir"
echo -e "Script path: $script_path"
echo -e ""

# RUN ==========================================================================

# Cycle on the files
for file in `ls $psl_dir`
do

	echo -e " · Working on '$file'..."
	suff=`echo $file | rev | cut -d "." -f -1 | rev`
	$script_path $seq_dir/ $psl_dir/  out_dir/ --psl $file \
		--suff $suff --ncores $threads --seq $seq_file & pid=$!
	wait $pid
	
done

# END ==========================================================================

################################################################################
