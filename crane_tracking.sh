#!/bin/bash

infile=../Videos/Narrawview_MovingBoom_Car\ only.mp4
generate_pcd = false

if $generate_pcd
then
    ../../im_synth downsampled.pcd 1000 1000 site.ppm depth.txt
    convert site.ppm -blu 0x1 site_blur.ppm
    convert site_blur.ppm site_blur.pgm
    ~/Downloads/siftDemoV4/sift < site_blur.pgm > site_blur.key
else
    cp ../car/site.ppm ./
    cp ../car/site_blur.key ./
fi

avconv -i $infile -ss 15 -t 35 -r 1 -vf scale=600:-1 %d.ppm
../../imsegment 1.ppm target_point.txt
for i in `seq 1 35`;do convert $i.ppm $i.pgm;done
for i in `seq 1 35`;do ~/Downloads/siftDemoV4/sift < $i.pgm > $i.key;done
for i in `seq 1 35`;do ~/Documents/vslam/improc 0.8 $i.site.match $i.key site_blur.key;done
../../solve_pnp_2d target_point.txt camera_location_2d.txt site.ppm
../../site_viewer_2d ../filtered.pcd camera_location_2d.txt 


convert -delay 50 {1..35}-seg.ppm anim_seg.gif
convert -delay 50 {1..35}-pcd.ppm anim_pcd.gif
