#!/usr/bin/python

import sys
import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scipy.misc

data_folder = sys.argv[1]
key = 1
target_point = numpy.loadtxt('%s/target_point.txt'%data_folder) 
label_point = numpy.loadtxt('%s/label_point.txt'%data_folder)
target_point = target_point[:,[1,2,4,5,7,8,10,11]].reshape((-1,4,2))
label_point = label_point[:,[1,2,4,5,7,8,10,11]].reshape((-1,4,2))
target_point = target_point[:,[0,1,3,2],:]
label_point = label_point[:,[0,1,3,2],:]

def displayImage(idx):
	print('displayImage %06d.png'%(idx))
	plt.cla()
	I = scipy.misc.imread('%s/%d.ppm'%(data_folder,idx))
	plt.imshow(I)
	ax.add_patch( patches.Polygon(target_point[idx-1],fill=False,edgecolor='#FF0000'))
	ax.add_patch( patches.Polygon(label_point[idx-1],fill=False,edgecolor='#FFFF00'))
	plt.draw()

def press(event):
	global key
	if event.key == 'left' and key > 1:
		key -= 1
		displayImage(key)
	if event.key == 'right' and key < len(label_point):
		key += 1
		displayImage(key)
	if event.key == 'z' and key > 1:
		key -= 10
		displayImage(key)
	if event.key == 'x' and key < len(label_point):
		key += 10
		displayImage(key)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
fig.canvas.mpl_connect('key_press_event', press)
displayImage(key)
plt.show()
