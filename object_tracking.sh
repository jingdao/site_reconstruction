#!/bin/bash

MATCH=/home/jd/Documents/vslam/improc
SIFT=/home/jd/Downloads/siftDemoV4/sift
ORIGIN=`pwd`
SCRIPT=/home/jd/Documents/site_reconstruction

if [ "$#" -ne "4" ] || ! [ -d $1 ] || ! [ -d $2 ] || ! [ -f $3 ] || ! [ -f $4 ]
then
	echo "./object_tracking.sh image_dir/ output_dir/ scenario.pcd object.pcd"
	exit
fi

generate_site=true
generate_images=true
generate_descriptors=true
generate_target=true
view_result=false
use_sift=false

image_dir=$1
output_dir=$2
scenario=$3
object=$4
width=600
height=400
ref_index=4
match_threshold=0.7
FEATURE=$SCRIPT/cvFeatures.py

if $generate_images
then
	rm $output_dir/*.ppm
	rm $output_dir/*.pgm
	j=0
	for f in `ls $image_dir/*.JPG | sort`
	do
		convert -resize "$width"x"$height"\! $f $output_dir/$j.ppm
		convert $output_dir/$j.ppm $output_dir/$j.pgm
		((j++))
	done
fi

if $generate_site
then
	rm $output_dir/site*
	$SCRIPT/im_synth $scenario $width $height $output_dir/site.ppm $output_dir/depth_buffer.txt $output_dir/cam_settings.txt
	convert $output_dir/site.ppm -blur 0x1 $output_dir/site_blur.ppm
	convert $output_dir/site_blur.ppm $output_dir/site_blur.pgm
fi

if $generate_descriptors
then
	rm $output_dir/*.key
	rm $output_dir/*.match
	j=0
	if $use_sift
	then
		$SIFT < $output_dir/site_blur.pgm > $output_dir/site_blur.key
		$SIFT < $output_dir/$ref_index.pgm > $output_dir/$ref_index.key
	else
		$FEATURE $output_dir/site_blur.pgm $output_dir/site_blur.key
		$FEATURE $output_dir/$ref_index.pgm $output_dir/$ref_index.key
	fi
	while true
	do
		if ! [ -f $output_dir/$j.pgm ]
		then
			break
		fi
		if $use_sift
		then
			$SIFT <	$output_dir/$j.pgm > $output_dir/$j.key
		else
			$FEATURE $output_dir/$j.pgm $output_dir/$j.key
		fi
		$MATCH $match_threshold $output_dir/$j.match $output_dir/$ref_index.key $output_dir/$j.key
		$MATCH $match_threshold $output_dir/$j.site.match $output_dir/$j.key $output_dir/site_blur.key
		((j++))
	done
fi

if $generate_target
then
	cd $output_dir/
	rm target_point.txt
	rm camera_location.txt
	$SCRIPT/match_image target_point.txt $ref_index.pgm
	$SCRIPT/solve_pnp $width $height target_point.txt depth_buffer.txt camera_location.txt
	cd $ORIGIN
fi

if $view_result
then
	$SCRIPT/site_viewer $scenario $object $output_dir/camera_location.txt $output_dir/cam_settings.txt
fi

