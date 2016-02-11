#!/bin/bash

MATCH=/home/jd/Documents/vslam/improc
SIFT=/home/jd/Downloads/siftDemoV4/sift
SCRIPT=`pwd`

if [ "$#" -ne "4" ] || ! [ -d $1 ] || ! [ -d $2 ] || ! [ -f $3 ] || ! [ -f $4 ]
then
	echo "./object_tracking.sh image_dir/ output_dir/ scenario.pcd object.pcd"
	exit
fi

generate_site=true
generate_images=true
generate_descriptors=true
generate_target=true

image_dir=$1
output_dir=$2
scenario=$3
object=$4
width=600
height=400

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
	./im_synth $scenario $width $height $output_dir/site.ppm $output_dir/depth_buffer.txt
	convert $output_dir/site.ppm -blur 0x1 $output_dir/site_blur.ppm
	convert $output_dir/site_blur.ppm $output_dir/site_blur.pgm
fi

if $generate_descriptors
then
	rm $output_dir/*.key
	rm $output_dir/*.match
	$SIFT < $output_dir/site_blur.pgm > $output_dir/site_blur.key
	j=0
	while true
	do
		if ! [ -f $output_dir/$j.pgm ]
		then
			break
		fi
		$SIFT <	$output_dir/$j.pgm > $output_dir/$j.key
		$MATCH $output_dir/$j.match $output_dir/0.key $output_dir/$j.key
		$MATCH $output_dir/$j.site.match $output_dir/$j.key $output_dir/site_blur.key
		((j++))
	done
fi

if $generate_target
then
	cd $output_dir/
	rm target_point.txt
	rm camera_location.txt
	$SCRIPT/match_image target_point.txt 0.pgm
	$SCRIPT/solve_pnp $width $height target_point.txt depth_buffer.txt camera_location.txt
	cd $SCRIPT
fi

./site_viewer $scenario $object $output_dir/camera_location.txt

