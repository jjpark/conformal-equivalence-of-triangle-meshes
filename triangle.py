import numpy as np
from scipy import integrate
class triangle:
	tri_number=0

	vertex=[] #global
	face=[] #global
	v_tri={} #global
	z=3
	edge_tri={}
	vert_wo_bound=[]


	def __init__(self, v1, v2, v3, num): #i:v1-1, j:v2-1, k:v3-1

		### vertexes are given as indices appears in .obj raw file
		self.new_lambda=[] # ij, jk, ki order
		self.new_lengths=[] #ij, jk, ki order
		self.angles=[] # order i j k
		self.ai=0 #angles
		self.aj=0
		self.ak=0
		self.li=0 #length opposite to vi -- l_jk
		self.lj=0 #l_ki
		self.lk=0 #l_ij
		self.uu=[]

		self.edges=[]
		self.text_coord=[0,0,0] # i, j, k np.array

		self.tri_number=num
		self.vi = triangle.vertex[v1-1]
		self.vj = triangle.vertex[v2-1]
		self.vk = triangle.vertex[v3-1]
		self.vertex_numbers=[v1-1,v2-1,v3-1]
		self.set_lengths()
		self.lengths=[self.lk, self.li, self.lj] #l_ij, l_jk, l_ki order
		self.new_lengths=[self.lk, self.li, self.lj]
		self.new_lambda=list(2*np.log(np.array(self.new_lengths)))
		self.set_angles()
		self.new_angles=[0,0,0] # angles at vi vj vk order

		a=[v1,v2,v3]
		for elem in a:
			if elem-1 in triangle.v_tri:
				triangle.v_tri[elem-1].append(num)
			else:
				triangle.v_tri[elem-1]=[num]

		self.edges=[(v1-1,v2-1),(v2-1,v3-1),(v3-1,v1-1)]
		for i in range(3):
			self.edges[i]=self.edge_order(self.edges[i])
		for edge in self.edges:
			if edge in triangle.edge_tri:
				triangle.edge_tri[edge].append(num)
			else:
				triangle.edge_tri[edge]=[num]

	def edge_order(self, edge):
		if (edge[0] > edge[1]):
			return (edge[1], edge[0])
		return edge

	def get_angle(self,i,j,k): 
		### get angle for vertex opposite of edge length i
		### i,j,k are edge lengths
		if i > j + k:
			print '!'
			print self.tri_number
			return np.pi
		elif j > i + k:
			return 0
		elif k > i + j:
			return 0
		return 2*np.arctan( ((i+j-k)*(i+k-j)/(j+k-i)/(i+j+k))**0.5 )

	def angl_at_abs_vert(self,ver):
		assert(ver in self.vertex_numbers)
		ind = self.vertex_numbers.index(ver)
		return self.new_angles[ind]

	def set_angles(self):
		self.ai=self.get_angle(self.li, self.lj, self.lk)
		self.aj=self.get_angle(self.lj, self.li, self.lk)
		self.ak=self.get_angle(self.lk, self.lj, self.li)
		self.angles=[self.ai, self.aj, self.ak]

	def set_new_angles(self):
		self.new_angles[0]=(self.get_angle(self.new_lengths[1],self.new_lengths[0],self.new_lengths[2]) )
		self.new_angles[1]=(self.get_angle(self.new_lengths[2],self.new_lengths[0],self.new_lengths[1]) )
		self.new_angles[2]=(self.get_angle(self.new_lengths[0],self.new_lengths[2],self.new_lengths[1]) )
		
	def set_lengths(self):
		self.li=self.get_length(self.vj,self.vk)
		self.lj=self.get_length(self.vi,self.vk)
		self.lk=self.get_length(self.vj,self.vi)

	def get_length(self,a,b):
		return ( (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 )**0.5

	def set_new_lambdas(self, u_i, u_j, u_k):
		### u_x is determined by the index of the face
		### determine tilda_lambdas based on given u values
		self.new_lambda[0]=(self.tilda_lambda(0, u_i, u_j) ) #new_lambda_ij
		self.new_lambda[1]=(self.tilda_lambda(1, u_j, u_k) )#new_lambda_jk
		self.new_lambda[2]=(self.tilda_lambda(2, u_k, u_i) )#new_lambda_ki

	def set_new_lengths(self, u_i, u_j, u_k):
		self.new_lengths[0] = np.e**((u_i+u_j)/2.0) * self.lengths[0] #l_ij

		self.new_lengths[1] = np.e**((u_j+u_k)/2.0) * self.lengths[1]
		self.new_lengths[2] = np.e**((u_k+u_i)/2.0) * self.lengths[2]

	def tilda_lambda (self, i, u1, u2):
		lam = 2*np.log(self.lengths[i])
		return lam + u1 + u2

	def lobachev(self, alpha):
		#assert(alpha>=0 and alpha<=np.pi)
		a=lambda t:np.log(2*np.sin(t) )
		b=-(integrate.quad(a,0,alpha)[0])
		#if(alpha>1.52 and alpha<1.62):
		#print (alpha, b)
		return b

	def set_u(self, u_i, u_j, u_k):
		self.uu=[u_i, u_j, u_k]
		### fuction f returns calculated f including the Lobachevsky function
		### given u[i] for each vertices in the triangle
		self.set_new_lambdas(u_i, u_j, u_k)
		self.set_new_lengths(u_i, u_j, u_k)
		self.set_new_angles()

	def f(self):
		
		first_term=0.5*( self.new_angles[0]*self.new_lambda[1]+\
			self.new_angles[1]*self.new_lambda[2]+self.new_angles[2]*self.new_lambda[0] )
		second_term=0

		for i in range(3):
			if self.new_angles[i] == 0 or self.new_angles[i] == np.pi:
				pass
			else:
				second_term+=self.lobachev(self.new_angles[i])
		return first_term+second_term

	def set_text_coord(self, edge, axis, other, cross):
		### returns the coordinates
		### edge is not list of absolute vertex_numbers of the line segment
		### axis is the coordinate of the first vertex --np.array
		### other is the coord of the second vertex --np.array
		### cross is the cross product value of the caller triangle
		assert((edge[0],edge[1]) in self.edges) #sanity check
		self.text_coord=[0,0,0]
		v1_ind=self.index_v(edge[0]) # index for origin
		v2_ind=self.index_v(edge[1]) # index for second v
		v3_ind=self.index_of_remainder_v(edge) # index for v to set now
		##TODO##
		##Testing function for text_coords and actual new_lengths##

		# set the two given coordinates
		self.text_coord[v1_ind]=axis
		self.text_coord[v2_ind]=other

		#length of edge from axis to other
		l_12=self.new_lengths[self.length_ind_frm_vert_ind(v1_ind,v2_ind)]
		"""
		if self.tri_number==36:
			print self.length_ind_frm_vert_ind(v1_ind,v2_ind)
		"""

		#print str(l_12-self.get_length(axis,other))+' '+str(self.tri_number)
		
		##TODO##
		##test if l_12 is equal to get_length(axis,other)##

		v1_trans=np.array([0,0,0])
		#coord of v2 when v1 is translated to the origin
		#this could be translated to be other again
		v2_trans=other-axis
		#coord of v2 when the edge lies on x-axis
		# this could be rotated to be v2_trans
		v2_trans_hori=np.array([l_12,0,0])
		## angle between x-axis and v2_trans
		rotated_angle=self.angl_rotate(v2_trans)
		#delta is positive when v2 is rotated counterclockwise from x-axis
		delta=self.delta(np.array([1,0,0]), v2_trans)
		# this coord is translated to near 0,0,0 and rotated so that the
		# argument edge is horizontal
		v3_trans_rotated=self.get_third_coord_trans([v1_ind, v3_ind],self.get_direction(cross) )
		#print self.get_length(v3_trans_rotated,v2_trans_hori)-self.new_lengths[self.length_ind_frm_vert_ind(v2_ind,v3_ind)]
		# this coord is rotated about the angle between the x-axis and the original edge argument
		v3_trans=self.rotate_vertex(v3_trans_rotated, delta*rotated_angle)
		# this coord is final
		v3_coord=v3_trans+axis
		#print self.get_length(v3_coord,other)-self.new_lengths[self.length_ind_frm_vert_ind(v2_ind,v3_ind)]
		self.text_coord[v3_ind]=v3_coord

		return self.text_coord

		"""     v3 -- this must be calculated
				/\
	           /  \
			  /	   \
			 /      \
			/        \
		  v1-----------v2   given edge, put this on x-axis and rotate later
		"""

	def angl_rotate(self,v_trans):
		# return the angle between the x axis and the translated v
		return np.arccos(np.dot(v_trans,np.array([1,0,0])) \
			/ np.linalg.norm(v_trans) )

	def index_of_remainder_v(self,edge):
		a=set([0,1,2])
		b=set([self.index_v(edge[0]), self.index_v(edge[1])])
		a-=b
		assert(len(a)==1)
		for i in a:
			return i

	def index_v(self,v):
		return self.vertex_numbers.index(v)

	def delta(self, hori_vec, orig):
		### see if original is left or right or opposite
		### if delta is positive orig is left of hori_vec so
		### the rotation angle should be positive as theta rotates
		### counter clockwise
		### if delta is negative rotation angle should be negative
		### if delta is zero rotation angle should be 180 degree
		delta=hori_vec[0]*orig[1]-hori_vec[1]*orig[0]
		if delta<0:
			return -1
		else:
			return 1

	def rotate_vertex(self, v_trans, angle):
		### get the vertex coord for a coordinate that are moved
		### near the origin (0,0,0)
		### in other words, the third vertex should be rotated to the right angle
		### the value returned must be translated to the original place
		### angle is negative when delta is negative
		return np.array([v_trans[0]*np.cos(angle)-v_trans[1]*np.sin(angle),\
			v_trans[0]*np.sin(angle)+v_trans[1]*np.cos(angle),0])

	def get_third_coord_trans(self, index, direction):
		### index is a list containing indexes of the triangle's vertices
		### the order is 1.origin vertex, 2.third vertex to draw
		### note the coordinate is get assuming a triangle with one
		### vertex at the origin and the second vertex on x axis.
		### if direction is negative the previous triangle was below x-axis and
		## the new coord is above x axis, otherwise below
		theta=direction*self.new_angles[index[0]] # angle at origin vertex
		len_ind=self.length_ind_frm_vert_ind(index[0],index[1]) 
		length=self.new_lengths[len_ind]# length between origin and third vertices
		return np.array([np.cos(theta)*length, np.sin(theta)*length,0])

	def length_ind_frm_vert_ind(self,index1,index2):

		if index1+index2==1:
			len_ind=0
		elif index1+index2==3:
			len_ind=1
		else:
			assert(index1+index2==2)
			len_ind=2
		return len_ind

	def get_direction(self, cross):
		### cross is the cross product of axis vector and another vector
		### from the origin index of the caller triangle
		### if cross is negative, the new triangle must be drawn above x axis
		### otherwise it should be drawn below x axis
		if cross<0:
			return 1
		else:
			return -1