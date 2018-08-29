#!/usr/bin/python

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')

crane_base = numpy.array([26.4,7.7,37.3])
load_width = 2
load_height = 0.2
virtual_f = 600
virtual_cam_pos = numpy.array([2,-1,60])
real_f = 300
cable_length = 9.38
#real_f = 600
#cable_length = 9.79

real_cam_pos = []
f = open('pnp.txt','r')
for l in f:
	if l.startswith('bestParam'):
		k = float(l.split()[1])
		theta = float(l.split()[2])
		Tx = float(l.split()[3])
		Ty = float(l.split()[4])
		x = virtual_cam_pos[0] - Tx * cable_length / k / virtual_f
		y = virtual_cam_pos[1] + Ty * cable_length / k / virtual_f
		z = cable_length * (1 - 1.0/k) + virtual_cam_pos[2]
		real_cam_pos.append([x,y,z])
f.close()

load_pos = numpy.loadtxt('camera_location_2d.txt')
load_pos = load_pos[:,[1,2,4,5,7,8,10,11]]
load_pos = load_pos.reshape((-1,4,2))

load_angle = []
for i in range(len(load_pos)):
	v1 = load_pos[i,1,:] - load_pos[i,0,:]
	v2 = load_pos[i,2,:] - load_pos[i,0,:]
	if numpy.linalg.norm(v1) > numpy.linalg.norm(v2):
		a = numpy.arctan2(-v1[1],v1[0])
	else:
		a = numpy.arctan2(-v2[1],v2[0])
	if a < 0:
		a += numpy.pi
	load_angle.append(a)	
load_pos = (load_pos - [500,500]) * cable_length / virtual_f + virtual_cam_pos[:2]
load_center = load_pos.mean(axis=1)
#print(real_cam_pos)
#print(load_pos)
#print(load_angle)
print('crane_base %f %f %f'%(crane_base[0],crane_base[1],crane_base[2]))

draw = False
for i in range(len(real_cam_pos)) if draw else range(0,len(real_cam_pos),int(numpy.ceil(len(real_cam_pos)/5.0))):
	print('boom_tip %f %f %f'%(real_cam_pos[i][0], real_cam_pos[i][1], real_cam_pos[i][2]))
	print('load %f %f %f %f'%(load_center[i][0],load_center[i][1],real_cam_pos[i][2]-cable_length, load_angle[i] / numpy.pi * 180))
	if draw:
		X = list([crane_base[0]]) + list([real_cam_pos[i][0]]) + list([load_center[i][0]])
		Y = list([crane_base[1]]) + list([real_cam_pos[i][1]]) + list([load_center[i][1]])
		V = [[2,1],[-2,1],[-2,-1],[2,-1],[2,1]]
		R = numpy.array([[numpy.cos(load_angle[i]), -numpy.sin(load_angle[i])], [numpy.sin(load_angle[i]), numpy.cos(load_angle[i])]])
		for v in V:
			v = R.dot(v)
			X.append(load_center[i][0] + v[0])
			Y.append(load_center[i][1] + v[1])
		X = numpy.array(X)
		Y = numpy.array(Y)
		Z = numpy.array(list([crane_base[2]]) + list([real_cam_pos[i][2]]) + list([real_cam_pos[i][2] - cable_length])*6)
		ax.clear()
		ax.scatter(crane_base[0],crane_base[1],crane_base[2],c='r',marker='o')
		max_range = numpy.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
		Xb = 0.5*max_range*numpy.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
		Yb = 0.5*max_range*numpy.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
		Zb = 0.5*max_range*numpy.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
		for xb, yb, zb in zip(Xb, Yb, Zb):
			   ax.plot([xb], [yb], [zb], 'w')
		ax.plot(X,Y,Z)
		plt.pause(0.01)
