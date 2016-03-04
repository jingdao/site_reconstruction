#!/usr/bin/python

import cv2
import sys

if len(sys.argv) < 3:
	print "./cvFeatures input.img out.key"
	sys.exit(1)

img = cv2.imread(sys.argv[1],0)
orb = cv2.ORB()
kp = orb.detect(img,None)
kp, des = orb.compute(img, kp)

f = open(sys.argv[2],"w")
descSize = len(des)
descDimension = len(des[0])
f.write(str(descSize)+" "+str(descDimension)+"\n")

for i in range(len(kp)):
	x = kp[i].pt[0]
	y = kp[i].pt[1]
	f.write(str(y)+" "+str(x)+" 0 0\n")
	for j in des[i]:
		f.write(str(j)+" ")
	f.write("\n")

