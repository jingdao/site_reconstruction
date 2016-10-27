#!/usr/bin/python

import sys
ratio = 0.8
source = None

if len(sys.argv) == 4:
	match = sys.argv[1]
	source = sys.argv[2]
	target = sys.argv[3]
elif len(sys.argv) == 5:
	ratio = float(sys.argv[1])
	target = sys.argv[2]
	id_start = int(sys.argv[3])
	id_end = int(sys.argv[4])
else:
	print 'Usage: kdtree_match match source target'
	print '       kdtree_match ratio target id_start id_end'
	sys.exit(1)


def get2NN(source, target):
#	print len(source),len(source[0])
#	print len(target),len(target[0])
#	import scipy.spatial
#	tree = scipy.spatial.KDTree(target)
#	print 'built tree'
#	dist,index = tree.query(source,k=2,eps=0.1,distance_upper_bound=len(source[0])*10)
	from sklearn.neighbors import NearestNeighbors
	nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(target)
	dist, index = nbrs.kneighbors(source)
	return dist,index

def parseKeyFile(filename):
	f = open(filename,'r')
	points = []
	desc = []
	l = f.readline()
	numPoints = int(l.split()[0])
	dimension = int(l.split()[1])
	l = f.read().split()
	for i in range(numPoints):
		p = (l[i*(dimension+4)] , l[i*(dimension+4)+1])
		points.append(p)
		d = [float(j) for j in l[i*(dimension+4)+4 : i*(dimension+4)+4+dimension]]
		desc.append(d)
	f.close()
	return (points,desc)

if source is not None:
	sourcePoints,sourceDesc = parseKeyFile(source)
	targetPoints,targetDesc = parseKeyFile(target)
	dist, index = get2NN(sourceDesc,targetDesc) 
	f = open(match,"w")
	for i in range(len(dist)):
		if dist[i][0] < dist[i][1] * ratio * ratio:
			x1 = sourcePoints[i][1]
			y1 = sourcePoints[i][0]
			x2 = targetPoints[index[i][0]][1]
			y2 = targetPoints[index[i][0]][0]
			f.write("0 %s %s 1 %s %s\n" % (x1,y1,x2,y2))
	f.close()
else:
	targetPoints,targetDesc = parseKeyFile(target)
	from sklearn.neighbors import NearestNeighbors
	nbrs = NearestNeighbors(n_neighbors=2, algorithm='brute').fit(targetDesc)
	for j in range(id_start,id_end+1):
		source = str(j)+".key"
		try:
			sourcePoints,sourceDesc = parseKeyFile(source)
		except IOError:
			continue
		dist, index = nbrs.kneighbors(sourceDesc)
		match = str(j)+".site.match"
		f = open(match,"w")
		numMatch = 0
		for i in range(len(dist)):
			if dist[i][0] < dist[i][1] * ratio:
#				print i,index[i],dist[i][0]*dist[i][0],dist[i][1]*dist[i][1]
				x1 = sourcePoints[i][1]
				y1 = sourcePoints[i][0]
				x2 = targetPoints[index[i][0]][1]
				y2 = targetPoints[index[i][0]][0]
				f.write("0 %s %s 1 %s %s\n" % (x1,y1,x2,y2))
				numMatch += 1
		f.close()
		print 'Wrote %d matches to %s' % (numMatch,match)

