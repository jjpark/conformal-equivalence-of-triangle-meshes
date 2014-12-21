import numpy as np
from scipy import integrate
from scipy import optimize
from triangle import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

vertices=[]
faces=[]
triangles=[]
boundary_vertices=[]
covered=[]
edge_tri={}
current=[]
corners=[]
corner=[]
target_angl=[]
weirdos=[]

import sys, os
file_name = sys.argv[1]
out_name = file_name[:file_name.index('.')] + '.txt'

sys.setrecursionlimit(10000) # 10000 is an example, try with different values
#egea 238
def main():
	print 'hello'
	lines = [line.strip() for line in open(file_name,'r+')]

	for i in range(len(lines)):
		line = lines[i]
		if ('v'==line[0]):
			line = str.split(line)
			line = line[1:]
			for i in range(3):
				line[i] = float(line[i])

			vertices.append(line)
		elif ('f'==line[0]):
			line = str.split(line)
			line = line[1:]
			for i in range(3):
				line[i] = int(line[i])
			faces.append(line)
		else:
			print 'corrupt data'

	triangle.vertex=vertices
	triangle.face=faces
	for i in range(len(faces)):
		covered.append(0)
		face=faces[i]
		triangles.append(triangle(face[0],face[1],face[2],i))

	edge_tri=triangle.edge_tri
	#t=triangles[0]
	#t.f(0.12,0.2,0.1)
	#print t.lengths
	#print t.new_lengths
	#print t.new_lambda
	#print np.array(t.angles)*180/np.pi
	#print np.array(t.new_angles)*180/np.pi
	"""
	for tri in triangles:
		#pass
		#print tri.vertex_numbers
		tri.f(0.12,0.2,0.1)
		tri.f(0.12,0.2,0.1)
		#print tri.new_lambda
	#print triangle.v_tri
	#	print str(tri.tri_number)+': '+str(tri.angles)
		#print np.array(tri.new_angles)*180/np.pi
	"""
	print '###'
	print len(vertices)
	for key in triangle.edge_tri:
		if len(triangle.edge_tri[key])!=2:
			for v in key:
				if v not in boundary_vertices:
					boundary_vertices.append(v)
				else:
					pass

	####### THIS IS ONLY FOR CONE SINGULARITIES ##########
	#boundary_vertices.append(17)
	#boundary_vertices.append(44)
	######################################################
	## this list contains index of internal vertices
	triangle.vert_wo_bound = range(len(vertices))

	for v in boundary_vertices:
		triangle.vert_wo_bound.remove(v)

	for v in range(len(vertices)):
		if v in corners:
			target_angl.append(np.pi/2)
		elif v in boundary_vertices:
			target_angl.append(np.pi)
		else:
			target_angl.append(2 * np.pi)
	print boundary_vertices


	"""
	a=[0, 2, 3, 6, 7, 237]
	a_sum=0
	for i in a:
		tri=triangles[i]
		print i
		print np.array(tri.angles)*180/np.pi
		a_sum+=tri.angl_at_abs_vert(0)
		print tri.angl_at_abs_vert(0)*180/np.pi
	
	pArray=np.array([])
	for i in range(len(vertices)):
		angle_sum=0
		for t in triangle.v_tri[i]:
			tri=triangles[t]
			angle_sum+=tri.angl_at_abs_vert(i)
		pArray=np.append(pArray, angle_sum)
	print pArray*180/np.pi
	""""""
	pArray=np.array([])

	for i in range(len(vertices)):
		angle_sum=0
		for t in triangle.v_tri[i]:
			tri=triangles[t]
			angle_sum+=tri.angl_at_abs_vert(i)
		pArray=np.append(pArray, 2*np.pi-angle_sum)
	print pArray"""
	
	guess=[]
	guess2=[]
	guess3=[]
	for i in range(len(triangle.vert_wo_bound)):
		guess.append(0.0)
		guess2.append(-10)
		guess3.append(100)
	

	"""
	guess=[]
	for i in range(len(vertices)):
		guess.append(0.0)
	"""
	#energy(guess)
	#energy(guess2)
	#energy(guess3)

	#ePrime(guess)
	#print boundary_vertices
	#print target_angl
	
	
	#energy(u_values)

	start = 0
	if file_name == 'face.obj':
		start = 600

	#tile_nei(triangles[0],(0,1),np.array([0,0,0]),np.array([triangles[0].new_lengths[0],0,0]),-1)
	
	#print triangle.v_tri
	#print lenboundary_vertices)
	
	m=optimize.fmin_l_bfgs_b(energy,guess,fprime=ePrime)
	
	print m[0]
	
	cong=open('temp','w')
	cnt = 0
	for item in m[0]:
		if cnt == 9:
			cong.write(str(item)+'\n')
			cnt = 0
		else:
			cong.write(str(item) + ' ')
			cnt+=1
	cong.close()
	for tri in triangles:
		if tri.new_angles[0]==0 or tri.new_angles[0]==np.pi:
			print 'nono'
	

	f=open("temp","r")
	a=f.readlines()
	u_values=[]
	for tens in a:
		if tens[0] != 'f':
			for ones in tens.split():
				u_values.append(float(ones))
	f.close()
	os.system('rm temp')
	ePrime(u_values)
	tile_not_recursive(triangles[start],triangles[start].edges[0],np.array([0,0,0]),np.array([triangles[start].new_lengths[0],0,0]),-1)
	while(len(current) != 0):
		for tri_num in current:
			tile_neighbors(triangles[tri_num])
	f=open(out_name,'w+')
	for tri in triangles:
		for i in range(3):
			f.write(str(tri.text_coord[i][0] )+'\n')
			f.write(str(tri.text_coord[i][1] )+'\n')
	f.close()
	
	
	#t = triangle([-0.206972, 0.0740737, 0.544664],[-0.398695, -0.0740794, 0.501089],\
	#	[-0.211327, -0.187367, 0.588234])
	#print t.vi

	# starting triangle: 86
	#tile_nei(triangles[0],(0,1),np.array([0,0,0]),np.array([triangles[0].new_lengths[0],0,0]),-1)
	#tile_not_recursive(triangles[1013],(484, 558),np.array([0,0,0]),np.array([triangles[1013].new_lengths[0],0,0]),-1)
	"""
	"""
	"""
	tile_not_recursive(triangles[238],(209, 210),np.array([0,0,0]),np.array([triangles[238].new_lengths[0],0,0]),-1)
	while(len(current) != 0):
		for tri_num in current:
			tile_neighbors(triangles[tri_num])
	
	
	f=open('egea_param.txt','w+')
	for tri in triangles:
		for i in range(3):
			f.write(str(tri.text_coord[i][0] * 150)+'\n')
			f.write(str(tri.text_coord[i][1] * 150)+'\n')
	f.close()
	
	x=[]
	y=[]
	for tri in triangles:
		for i in range(3):
			x.append(tri.text_coord[i][0] * 100)
			y.append(tri.text_coord[i][1] * 100)

	plt.scatter(x,y,marker='+',s=150,linewidths=4 ,c=y, cmap=plt.cm.coolwarm)
	plt.show()"""
	
	"""
	e=open('egea_param.txt','r')
	h=open('egea_sphere.txt','w')
	stereographed = []
	while(True):
		a = e.readline()
		b = e.readline()
		if (not a or not b):
			break

		param = ( float(a.split()[0]), float(b.split()[0]) )
		stereographed.append(stereo_project(param) )


	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	xs,ys,zs=[],[],[]
	for point in stereographed:
		for i in range(3):
			h.write(str(point[i]) + '\n' )
		xs.append(point[0])
		ys.append(point[1])
		zs.append(point[2])
	ax.scatter(xs, ys, zs, c='r', marker='o')
	plt.show()
	"""

def stereo_project(par):
	x = par[0]
	y = par[1]
	result = []
	result.append(2*x/(1+x**2+y**2))
	result.append(2*y/(1+x**2+y**2))
	result.append((-1+x**2+y**2)/(1+x**2+y**2))
	return result

def tile_neighbors(tri):
	assert(covered[tri.tri_number] == 1)
	# this triangle should be in current queue
	assert(tri.tri_number in current)
	# now it's not current
	current.remove(tri.tri_number)
	#get coords
	coords = tri.text_coord
	for i in range(3):
		e=tuple(tri.edges[i])
		if(len(triangle.edge_tri[e])==1):
			continue
		# get the tri_number of neighbor triangle of this edge
		neighbor= triangle.edge_tri[e][(triangle.edge_tri[e].index(tri.tri_number)+1) % 2]

		if(covered[neighbor]==1):
			continue
		
		v1=coords[tri.index_v(e[0])] ##origin vertex
		v2=coords[tri.index_v(e[1])] ##other vertex
		# find the vertex that's not in the current edge
		v3=coords[tri.index_of_remainder_v(e)]
		# calculate the cross product
		cross_z_value=np.cross((v2-v1),(v3-v1))[2]
		
		tile_not_recursive(triangles[neighbor],e,v1,v2,cross_z_value)

	return


def tile_not_recursive(tri, edge, origin, other_v, cross):
	covered[tri.tri_number] = 1
	# add the triangle to the current breadth first search queue
	current.append(tri.tri_number)
	# set text_coord for this triangle
	tri.set_text_coord(edge, origin, other_v, cross)

def tile_nei(tri, edge, origin, other_v, cross):
	covered[tri.tri_number]=1
	#set coordinate for this triangle
	coords=tri.set_text_coord(edge,origin,other_v, cross)
	#iterate through each edge in the tri
	for i in range(3):
		e=tuple(tri.edges[i])
		if(len(triangle.edge_tri[e])==1):
			continue
		# get the tri_number of neighbor triangle of this edge
		neighbor= triangle.edge_tri[e][(triangle.edge_tri[e].index(tri.tri_number)+1) % 2]

		if(covered[neighbor]==1):
			continue
		
		v1=coords[tri.index_v(e[0])] ##origin vertex
		v2=coords[tri.index_v(e[1])] ##other vertex
		# find the vertex that's not in the current edge
		v3=coords[tri.index_of_remainder_v(e)]
		# calculate the cross product
		cross_z_value=np.cross((v2-v1),(v3-v1))[2]
		
		tile_nei(triangles[neighbor],e,v1,v2,cross_z_value)
	return
def energy_rec(u):

	ener=0

	for tri in triangles:
		v_num = tri.vertex_numbers # e.g. 0,20,1001

		tri.set_u(u[v_num[0]],u[v_num[1]],u[v_num[2]])
		ener += - 0.5*np.pi*(u[v_num[0]] + u[v_num[1]] + u[v_num[2]])

	for tri in triangles:
		ener += tri.f()

	for i in range(len(u)):
		angl = target_angl[i]
		ener += 0.5*angl * u[i]

	print 'energy: ' + str(ener)
	return ener

def energy(u):
	ener=0

	for tri in triangles:
		v_num = tri.vertex_numbers
		us = []

		for vert_ind in v_num:
			# if the vertex is boundary
			if vert_ind in boundary_vertices:
				# set the value to 0
				us.append(0.0)
			else:
				# find the index of the vertex in vert_without_bound list
				us.append(u[ triangle.vert_wo_bound.index(vert_ind) ])
		tri.set_u(us[0], us[1], us[2])
		ener += - 0.5*np.pi*(us[0] + us[1] + us[2])

	for tri in triangles:
		ener += tri.f()
		#print tri.uu

	for u_i in u:
		ener += np.pi * u_i
	print 'energy: ' + str(ener)
	return ener

def energy_cone(u):
	ener=0

	for tri in triangles:
		v_num = tri.vertex_numbers
		us = []

		for vert_ind in v_num:
			# if the vertex is boundary
			if vert_ind in boundary_vertices:
				# set the value to 0
				us.append(0.0)
			else:
				# find the index of the vertex in vert_without_bound list
				us.append(u[ triangle.vert_wo_bound.index(vert_ind) ])
		tri.set_u(us[0], us[1], us[2])
		ener += - 0.5*np.pi*(us[0] + us[1] + us[2])

	for tri in triangles:
		ener += tri.f()
		#print tri.uu

	for u_i in u:
		ener += np.pi * u_i
	ener -= np.pi * 0.5 * (u[corner[0]]+u[corner[1]])
	print 'energy: ' + str(ener)
	return ener

def ePrime_rec(u):
	#print u
	print 'u_size: ' + str(np.linalg.norm(u))
	bound=[]
	for tri in triangles:
		v_num = tri.vertex_numbers # e.g. 0,20,1001

		tri.set_u(u[v_num[0]],u[v_num[1]],u[v_num[2]])

	pArray=np.array([])
	for v in range(len(vertices)):
		angle_sum = 0
		angl = target_angl[v]
		for t in triangle.v_tri[v]: #list of triangle indexes that contain a vertex v
			tri = triangles[t]
			angle_sum += tri.angl_at_abs_vert(v)
			if v in boundary_vertices:
				bound.append(angle_sum)
		pArray = np.append(pArray, (angl - angle_sum) * 0.5)
	ang_size=0
	for elements in pArray:
		ang_size+=abs(elements)
	#print pArray
	flat = ang_size/(len(pArray))
	max_ang = max(pArray)
	print 'flat: ' + str(flat)
	print 'max angle: '+str(max_ang)
	print 'angle at boundaries: '
	#print bound
	print str(pArray[corners[0]])+' '+str(pArray[corners[1]])+' '+str(pArray[corners[2]])
	#print 'size :' + str(np.linalg.norm(pArray) )

	return pArray

def ePrime_cone(u):
	#print u
	print 'u_size: ' + str(np.linalg.norm(u))

	for tri in triangles:
		# index of vertices in a triangle
		v_num = tri.vertex_numbers
		us = []
		for vert_ind in v_num:
			# if the vertex is boundary
			if vert_ind in boundary_vertices:
				# set the value to 0
				us.append(0)
			else:
				# find the index of the vertex in vert_without_bound list
				us.append(u[ triangle.vert_wo_bound.index(vert_ind) ])
		tri.set_u(us[0], us[1], us[2])

	pArray=np.array([])
	for v in triangle.vert_wo_bound:
		angle_sum = 0
		for t in triangle.v_tri[v]: #list of triangle indexes that contain a vertex v
			tri = triangles[t]
			angle_sum += tri.angl_at_abs_vert(v)
		pArray = np.append(pArray, (2*np.pi - angle_sum) * 0.5)
	pArray[corner[0]] -= np.pi * 0.5
	pArray[corner[1]] -= np.pi * 0.5
	ang_size=0
	for elements in pArray:
		ang_size+=abs(elements)
	#print pArray
	flat = ang_size/(len(pArray))
	max_ang = max(pArray)
	print 'flat: ' + str(flat)
	print 'max angle: '+str(max_ang)
	print 'size :' + str(np.linalg.norm(pArray) )
	"""
	if flat < 0.0005 and max_ang < 0.0055:
		cong=open('finally','w')
		cong.write('flat: ' + str(flat) + ' max: '+str(max_ang)+'\n')
		cnt = 0
		for item in u:
			if cnt == 9:
				cong.write(str(item)+'\n')
				cnt = 0
			else:
				cong.write(str(item) + ' ')
				cnt+=1
		cong.close()
		"""
	
	#print pArray

	return pArray


def ePrime(u):
	#print u
	print 'u_size: ' + str(np.linalg.norm(u))

	for tri in triangles:
		# index of vertices in a triangle
		v_num = tri.vertex_numbers
		us = []
		for vert_ind in v_num:
			# if the vertex is boundary
			if vert_ind in boundary_vertices:
				# set the value to 0
				us.append(0)
			else:
				# find the index of the vertex in vert_without_bound list
				us.append(u[ triangle.vert_wo_bound.index(vert_ind) ])
		tri.set_u(us[0], us[1], us[2])

	pArray=np.array([])
	for v in triangle.vert_wo_bound:
		angle_sum = 0
		for t in triangle.v_tri[v]: #list of triangle indexes that contain a vertex v
			tri = triangles[t]
			angle_sum += tri.angl_at_abs_vert(v)
		pArray = np.append(pArray, (2*np.pi - angle_sum) * 0.5)
	ang_size=0
	for elements in pArray:
		ang_size+=abs(elements)
	#print pArray
	flat = ang_size/(len(pArray))
	max_ang = max(pArray)
	print 'flat: ' + str(flat)
	print 'max angle: '+str(max_ang)
	print 'size :' + str(np.linalg.norm(pArray) )
	"""
	if flat < 0.0005 and max_ang < 0.0055:
		cong=open('finally','w')
		cong.write('flat: ' + str(flat) + ' max: '+str(max_ang)+'\n')
		cnt = 0
		for item in u:
			if cnt == 9:
				cong.write(str(item)+'\n')
				cnt = 0
			else:
				cong.write(str(item) + ' ')
				cnt+=1
		cong.close()
		"""
	
	#print pArray

	return pArray

def p_hessian(u,p):

	for tri in triangles:
		# index of vertices in a triangle
		v_num = tri.vertex_numbers
		us = []
		for vert_ind in v_num:
			# if the vertex is boundary
			if vert_ind in boundary_vertices:
				# set the value to 0
				us.append(0)
			else:
				# find the index of the vertex in vert_without_bound list
				us.append(u[ triangle.vert_wo_bound.index(vert_ind) ])
		tri.set_u(us[0], us[1], us[2])

	hArray = np.array([])

	for v_index in range(len(triangle.vert_wo_bound)):
		summation = 0
		# set of edges for a vertex
		v_edge = []
		# two triangles for an edge
		tris = []

		for edge in triangle.edge_tri:
			abs_index = triangle.vert_wo_bound[v_index]
			if triangle.vert_wo_bound[v_index] in edge: # real v_index
				if edge[(edge.index[v_index] + 1)%2] not in boundary_vertices:
					v_edge.append(edge)

		for e_ij in v_edge: #for an edge that contains v_i
			v_j = e_ij[ (e_ij.index(v_index) + 1) % 2 ] # get index for v_j
			tris = triangle.edge_tri[e_ij] # two or one triangles
			w_ij = 0
			for t in tris: # boundary edges have one triangles
				t_ijk = triangles[t]
				print 'cotangent'
				w_ij += (np.tan(t_ijk.new_angles[t_ijk.index_of_remainder_v (e_ij) ] ) )**(-1)
			summation += w_ij * (p[v_index] - p[v_j])

	## for each vertex number
	for v_i in range(len(vertices)):
		# set u values
		for t in triangle.v_tri[v_i]:
			tri=triangles[t]
			vi=tri.vertex_numbers[0]
			vj=tri.vertex_numbers[1]
			vk=tri.vertex_numbers[2]
			tri.set_u(u[vi], u[vj], u[vk])
		##################################

		summation = 0
		# set of edges for a vertex
		v_edge = []
		# two triangles for an edge
		tris = []
		### find v_i in the edge
		for key in triangle.edge_tri: # key is edge
			if v_i in key:
				# now we have edges that contain a vertex
				v_edge.append(key)
		for e_ij in v_edge: #for an edge that contains v_i
			v_j = e_ij[ (e_ij.index(v_i) + 1) % 2 ] # get index for v_j
			tris = triangle.edge_tri[e_ij] # two or one triangles
			w_ij = 0
			for t in tris: # boundary edges have one triangles
				t_ijk = triangles[t]
				print 'cotangent'
				w_ij += (np.tan(t_ijk.new_angles[t_ijk.index_of_remainder_v (e_ij) ] ) )**(-1)
			summation += w_ij * (p[v_i] - p[v_j])
		hArray = np.append(hArray, summation)
	print hArray
	return hArray


def edge_order(self, edge):
	if (edge[0] > edge[1]):
		return (edge[1], edge[0])
	return edge

def bound_to_zero(index):
	if index in boundary_vertices:
		return 0
	else:
		return index

if __name__ == "__main__":
    main()

