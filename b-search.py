# -*- coding: utf-8 -*-
import math
import random
import numpy as np
import sys
import copy
import time
from gurobipy import *

def lsc(m, n, d, sigma):
	model = Model("set-cover")
	model.setParam('OutputFlag', False )
	y, x = {}, {}

	z = model.addVar(obj=1, vtype="C", name="z")

	for j in range(m):
		y[j] = model.addVar(obj=0, vtype="B", name="y[%s]"%j)

	#d_ord[0] = [2.3, 2.5, 5.6, 6.7]

	#for i in range(m):
	#	coef = [1 for j in range(m)]
	#	var = [y[j] for j in range(m)]
	#	model.setObjective(LinExpr(coef,var), GRB.MINIMIZE)
	model.update()

	model.setObjective(z, GRB.MINIMIZE)

	model.update()

	for i in range(n):
		coef = [0 for k in range(m)]
		for j in range(n):
			if d[i,j] <= sigma:
				coef[j] = 1

		#coef = [1 for j in range(m)]
		var = [y[q] for q in range(m)]
		model.addConstr(LinExpr(coef,var), ">=", 1, name="Assign[%s]"%i)

	for i in range(n):
		coef = [1 for j in range(m)]
		var = [y[j] for j in range(m)]
		model.addConstr(LinExpr(coef,var), "=", z, name="Assign2[%s]"%i)

	model.update()
	model.__data = x,y

	return model

def sort_index(my_list):
	s = [i[0] for i in sorted(enumerate(my_list), key=lambda x:x[1])]
	return s

def get_ub0(d):
	max = np.zeros(len(d))
	
	for j in range(0, len(d)):
		for i in range(0, len(d)):
			if d[i,j] > max[j]:
				max[j] = d[i,j]

	return min(max)

def greedy_sc(d, D, h):
	N = len(d)

	# G pode ser uma lista de adjacencias
	G = np.matrix(N)
	G = d.copy()

	for i in range(0, N):
		for j in range(0, N):
			if d[i][j] > D[h]:
				G[i][j] = 0

	grau = np.zeros(N)

	for i in range(0, N):
		for j in range(0, N):
			if G[i][j] > 0:
				grau[i] = grau[i] + 1

	graus_index = sort_index(grau)

	i_x = N

	covered = np.zeros(N, dtype=np.int)
	y_g = np.zeros(N, dtype=np.int)

	while i_x >= 0:
		i_x = i_x-1
		i = graus_index[i_x]
		num_uncovored_ng = 0

		for k in range(0, N):
			if (G[i][k] > 0 and covered[k] == 0) or (covered[k] == 0 and k == i):
				num_uncovored_ng = num_uncovored_ng + 1

		if num_uncovored_ng > 0: # there are uncovered neighbors
			covered[i] = 1
			y_g[i] = 1

			for j in range(0, N):
				if G[i][j] > 0:
					covered[j] = 1

	centers = np.zeros(N, dtype=np.int)
	centers2 = np.zeros(N, dtype=np.int)
	max_min_dist = 0

	# assign client to facilities
	for j in range(0, N):
		min_dist = float('inf')
		#min_dist2 = float('inf')
		#min_dist = D[len(D)-1]
		#min_dist2 = D[len(D)-1]
		for i in range(0, N):
			if d[i][j] < min_dist and y_g[i]==1:
				#min_dist2 = min_dist
				min_dist = d[i][j]
				#centers2[j] = centers[j]
				centers[j] = i # client j is served by facility i

		if min_dist > max_min_dist:
			max_min_dist = min_dist	

	# close facilities uncovering any client
	y_g = np.zeros(N, dtype=np.int)
	for i in range(0, N):
		y_g[centers[i]] = 1

	for j in range(0, N):
		min_dist2 = float('inf')
		for i in range(0, N):
			if d[i][j] < min_dist2 and y_g[i]==1 and i != centers[j]:
				min_dist2 = d[i][j]
				centers2[j] = i

	# count the number of opened facilities after closing redundant facilities
	num_opened = 0
	for i in range (0, N):
		if y_g[i] == 1:
			num_opened = num_opened + 1

	return y_g, centers, centers2, num_opened, max_min_dist

'''
def find_best_center(centers, j, f): # pass current centers, client j, forbidden center f
	min_dist = float('inf')
	best_c = f
	
	for i in range(0, len(centers)):
		c = centers[i]
		if c != f and d[c][j] < min_dist:
			min_dist = d[c][j]
			best_c = c # find the best center to have client j re-allocated

	return best_c
'''

def close_facility(d, centers, centers2, P):
	# this is step 4!!
	N = len(centers)
	current_centers = dict()

	for j in range(0, N):
		if centers[j] not in current_centers:
			current_centers[centers[j]] = list()
		
		current_centers[centers[j]].append(j)

	max_f = {}
	min_f = float('inf')

	for f in current_centers:
		second_dist = list()
		for client in current_centers[f]:
			second_dist.append(d[centers2[client]][client])

		max_f[f] = max(second_dist)
		
		if max_f[f] < min_f:
			best_fac_to_close = f
	
	ub = 0

	for j in range(0, N):
		if centers[j] == best_fac_to_close:
			centers[j] = centers2[j]
		if d[centers[j]][j] > ub:
			ub = d[centers[j]][j]

	for j in range(0, N):
		min_dist2 = float('inf')
		for i in range(0, N):
			if d[i][j] < min_dist2 and i != centers[j] and i != best_fac_to_close and i in current_centers:
				min_dist2 = d[i][j]
				centers2[j] = i

	return ub, centers, centers2

def bsearch(N, M, P, UB, D, d):
	K = len(D)
	head = -1
	tail = K
	LB = 0
	num_opened = 0
	
	while True:
		# step 0
		if head >= (tail-1):
			LB = D[tail]
			print '********** LB =', LB
			print '********** UB =', UB
			break

		# step 1
		h = int((head + tail)/2)

		# step 2
		y_g, centers, centers2, num_opened, max_min_dist = greedy_sc(d, D, h)

		# step 3
		if num_opened <= P:
			tail = h
			UB = D[h]
		
		else:
			if num_opened == 1:
				break

			# step 4
			while (num_opened > P):
				ub_, centers, centers2 = close_facility(d, centers, centers2, P)
				num_opened = num_opened-1

			#print 'ub_, D[h], UB', ub_, D[h], UB

			if ub_ < UB:
				UB = ub_

			# step 5
			model = lsc(N, M, d, D[h])
			model.optimize()
			x,y = model.__data
			
			# step 6
			if model.ObjVal > P:
				head = h
			else:
				tail = h
				while UB < D[tail]:
					tail = tail - 1

'''
def calc_distance(x, y):
	return math.sqrt(math.pow((x[0]-y[0]),2) + math.pow((x[1]-y[1]),2))

# Transforma em matriz de facility x cliente
# Facility = linha
# Cliente = coluna
def transform_sc(v_fac, v_cli):

	assert(len(v_fac) == len(v_cli))

	d = np.zeros((len(v_fac),len(v_cli)))

	for i in range(0,len(v_fac)):
		for j in range(0,len(v_cli)):
			d[i,j] = calc_distance(v_fac[i],v_cli[j])
	
	return d;

'''

def distance(x1, y1, x2, y2):
	return math.sqrt((x2-x1)**2 + (y2-y1)**2)

def make_data(n):
    # positions of the points in the plane
    fac_x = [random.random() for i in range(n)]
    fac_y = [random.random() for i in range(n)]

    cli_x = [random.random() for i in range(n)]
    cli_y = [random.random() for i in range(n)]

    c = np.zeros((n,n))
    max_c = 0
    for i in range(n):
        for j in range(n):
            c[i][j] = distance(fac_x[i],cli_x[i],fac_y[j],cli_y[j])
            max_c = max(c[i,j],max_c)

    return c, max_c

if __name__ == "__main__":
	random.seed(67)
	n = 200
	d, max_c = make_data(n)
	p = 2
	m = n
	
	start = time.time()
	ub0 = get_ub0(d)
	D = set()
	
    # creating D list    
	for i in range(0, len(d)):
    # update len(d) if M != N
		for j in range(0, len(d)):
			D.add(d[i,j])

	s = sorted(D)
	bsearch(len(d), len(d), p, ub0, s, d)

	print 'time =', time.time() - start