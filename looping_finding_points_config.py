import numpy as np
import matplotlib.pyplot as plt
import math
import csv

# This code was used to find the configuration of N points
# such that the distance to their four closest neighbors is the same
# in a square torus with side ws = 20


def distances_torus(x, y, ws):
	x1, x2 = np.meshgrid(x, x)
	dx = abs(x1 - x2)
	dx = np.where(dx > ws/2, ws - dx, dx)
	y1, y2 = np.meshgrid(y, y)
	dy = abs(y1 - y2)
	dy = np.where(dy > ws/2, ws - dy, dy)
	dists = np.sqrt(dx**2 + dy**2)       
	return dists

# Transform the parameters aplying a rotation, translation and scaling according to parameters
# levelx, levely, centerx, centery, scale
def transform_coordinates(cordinates, parameters):
	levelx, levely, centerx, centery, scal = parameters
	# co/ca, co = levy = levely, ca = levx = levely + levelx
	
	# The rotation angle is that of a triangle with sidey = levely + levelx, 
	# and sidey  = levely. This arrange causes the rotation to move a point at levely, (levely+levelx) 
	# to the x axis, and since there is symmetry, we only search for angles from 0 to pi/4. 

	theta = -math.atan(levely/(levely+levelx))
	
	Rot = np.array([[math.cos(theta), math.sin(theta), 0], 
	                [-math.sin(theta),  math.cos(theta), 0],
	                [0, 0, 1]])

	transit_ida = np.array([[1, 0, centerx],
	                        [0, 1, centery],
	                        [0, 0, 1]])

	transit_vuelta = np.array([[1, 0, -centerx],
	                        	[0, 1, -centery],
	                        	[0, 0, 1]])

	scale = np.array([[scal, 0, 0],
	                  [0, scal, 0],
	                  [0, 0, 1]])

	transformed = np.matmul(transit_vuelta, cordinates) 
	transformed = np.matmul(scale, transformed)
	transformed = np.matmul(Rot, transformed)
	transformed = np.matmul(transit_ida, transformed)
	transformed=transformed.transpose()

	return transformed


ws = 20 #world size 

# Add square numbers to not find trivial cases
configs_found = [int(i**2) for i in range(1, 7) ]

# Search inside a range of conditions looping over scaling vand rotations.
# The points are transformed, and only the ones inside the 2D world are considered.
# If the sum of the variance of the distance to each point to its four closest neighbours
# is low, add that transformation to check later

def search_inside_limits(number_points_per_side, ws, positions_x, positions_y, configs_found):
	x = np.linspace(-30, 30, number_points_per_side + 1)
	y = np.linspace(-30, 30, number_points_per_side + 1)
	X, Y = np.meshgrid(x[:-1], y[:-1])
	Xf = X.flatten(order = 'A')
	Yf = Y.flatten(order = 'A')
	cordinates = np.array([Xf, Yf, [1]*len(Xf) ])
	found = []
	for px in range(positions_x[0], positions_x[1]): 
		for py in range(positions_y[0], positions_y[1]):
			for s in np.linspace(2, 7, 3000):
				# co/ca, co = levy = levely, ca = levx = levely + levelx
				parameters = [px, py, 0, 0,  s]
				transformed = transform_coordinates(cordinates = cordinates, 
													parameters = parameters)
				transformed = transformed[np.all(transformed[:,0:2] >= 0, axis = 1), :]
				transformed = transformed[np.all(transformed[:,0:2] < ws, axis = 1), :]

				dist_mat = distances_torus(transformed[:,0], transformed[:,1], ws)

				if np.shape(dist_mat)[0] >= 4 and np.shape(dist_mat)[0] <= 40 :
					sum_cost = 0		
					for b in dist_mat:
						fila = np.sort(b)
						sum_cost += np.var(fila[1:5])

					if sum_cost < 0.1 and  np.shape(dist_mat)[0]<= 50 and (not np.shape(dist_mat)[0] in configs_found):
						print(f"Found configuration with {np.shape(dist_mat)[0]} points")
						found.append([np.shape(dist_mat)[0], parameters, number_points_per_side ])
						configs_found.append(np.shape(dist_mat)[0])
	return found


search_conditons = [
			[39, [1,13], [1,13]],
			[61, [1,13], [1,13]],
			[60, [1,13], [1,13]],
			[40, [1,9],  [1,9]], 
			[40, [9,18], [1,9]],
			[40, [1,9],  [9,16]],
			[41, [1,9],  [1,9]],
			[41, [1,9],  [9,16]] 
]


found=[]
for case in search_conditons:
	#number_points_per_side, ws, positions_x, positions_y, configs_found 
	esta_vez = search_inside_limits(number_points_per_side = case[0], 
									ws = ws, 
									positions_x = case[1], 
									positions_y = case[2], 
									configs_found = configs_found )
	found = found + esta_vez
	configs_found.append([i[0] for i in esta_vez])


rows = []
for prev_result in found:
	this_scale = prev_result[1][4]
	best_sofar = 20000
	this_cords = []
	x = np.linspace(-30, 30, prev_result[2] + 1)
	y = np.linspace(-30, 30, prev_result[2] + 1)

	X, Y = np.meshgrid(x[:-1], y[:-1])

	Xf = X.flatten(order = 'A')
	Yf = Y.flatten(order = 'A')

	cordinates = np.array([Xf, Yf, [1]*len(Xf) ])
	for i in np.linspace(0.7 * this_scale, 1.3 * this_scale, 10000):
		parameters = prev_result[1][0:4] + [i] # change only the scaling to be more precise
		transformed = transform_coordinates(cordinates = cordinates, parameters = parameters)
		transformed = transformed[np.all(transformed[:,0:2] >= 0, axis = 1) , :]
		transformed = transformed[np.all(transformed[:,0:2] < ws, axis = 1) , :]
		dist_mat = distances_torus(transformed[:,0], transformed[:,1], ws)
		sum_cost = 0	
		for i in dist_mat:
			fila = np.sort(i)
			sum_cost += np.var(fila[1:5]) 

		if sum_cost < best_sofar and np.shape(dist_mat)[0] == prev_result[0] :
			best_sofar = sum_cost
			this_cords = transformed

	# Veryfy 
	dist_mat_final_test = distances_torus(this_cords[:,0], this_cords[:,1], ws)
	si_queda = True
	for j in dist_mat_final_test:
		fila = np.sort(j)
		if any(abs(np.diff(fila[1:5])) > 0.001):
			si_queda = False
	if si_queda:
		rows.append([prev_result[0]]+[str(h[0])+":"+str(h[1]) for h in this_cords])

writer = csv.writer(open("point_configurationsx.csv", mode="w", newline=""))
writer.writerows(rows)
print("Fin")
