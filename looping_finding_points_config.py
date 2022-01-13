
#C:\Users\hugos\Desktop\phd\data\renton_model\experimento
import numpy as np
import matplotlib.pyplot as plt
import math
import csv

# this code was used to fin the configuration of N points
# such that the distance to their four closest neighbors is the same
# in a torus world


def distances_torus(x, y, ws):
	x1, x2 = np.meshgrid(x, x)
	dx = abs(x1 - x2)
	### this wraps the interaction effects
	dx = np.where(dx > ws/2, ws - dx, dx)
	y1, y2 = np.meshgrid(y, y)

	dy = abs(y1 - y2)
	dy = np.where(dy > ws/2, ws - dy, dy)

	dists = np.sqrt(dx**2 + dy**2)       
	return dists

def transform_coordinates(cordinates, parameters, lvl):
	theta_N, centerx, centery, scal = parameters

	theta = -math.atan(lvl/(lvl+theta_N))
	
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


	transformed = np.matmul(transit_vuelta, cordinates) #llevar al inicio

	transformed = np.matmul(scale, transformed) # escalar

	transformed = np.matmul(Rot, transformed) # rotar en contra de sentido de reloj

	transformed = np.matmul(transit_ida, transformed)# regresar a centro original

	transformed=transformed.transpose()

	return transformed

def plot_transformation(cordinates, parameters):
	
	transformed = transform_coordinates(cordinates=cordinates, parameters=parameters)
	
	transformed = transformed[np.all(transformed[:,0:2] >= 0, axis=1) , :]
	transformed = transformed[np.all(transformed[:,0:2] < ws, axis=1 ) , :]

	plt.scatter(transformed[:,0], transformed[:,1])
	for i in range(len(transformed[:,0])):
		plt.annotate(str(i), xy=(transformed[i,0] , transformed[i,1]) )
	plt.xlim([0, ws])
	plt.ylim([0, ws])
	plt.show()


ws = 20
n_row_fst = 41
ya_estan =[ 4 , 9 ,16 ,25 ,36 ,49 , 8 ,18 ,32, 50, 72, 98, 20,5,13,25,40,10,34,17,29,26,37,8,15,24,16,9,36,35]

def search_inside_limits(n_row_fst, ws, level, position, ya_estan ):
	x = np.linspace(-30, 30, n_row_fst + 1)
	y = np.linspace(-30, 30, n_row_fst + 1)

	X, Y = np.meshgrid(x[:-1], y[:-1])

	Xf = X.flatten(order='A')
	Yf = Y.flatten(order='A')

	cordinates = np.array([Xf, Yf, [1]*len(Xf) ])
	#print(cordinates)

	# plt.scatter(cordinates[0,:], cordinates[1,:])
	# for i in range(len(cordinates[0,:])):
	# 	plt.annotate(str(i), xy=(cordinates[0,i] , cordinates[1,i]) )
	# 	plt.xlim([0, ws])
	# 	plt.ylim([0, ws])
	# plt.show()
	#plot_transformation(cordinates = coordinates, parameters= [0.8, 6.666667, 6.666667, 1.5])# para compara con 
	#plot_transformation(cordinates = coordinates, parameters= [math.pi/4, 0, 0, 1.5])# para compara con R

	found = []

	for j in range(position[0], position[1]):
		print(j)
		for lev in range(level[0], level[1]):
			for i in np.linspace(2, 7, 3000):# poner 5000
				parameters = [j, 0, 0 ,  i]
				transformed = transform_coordinates(cordinates=cordinates, parameters=parameters, lvl=lev)
				transformed = transformed[np.all(transformed[:,0:2] >= 0, axis=1) , :]
				transformed = transformed[np.all(transformed[:,0:2] < ws, axis=1 ) , :]
				dist_mat = distances_torus(transformed[:,0], transformed[:,1], ws)
				#plot_transformation(cordinates = cordinates, parameters= parameters)# para compara con R



				if np.shape(dist_mat)[0] >= 4 and np.shape(dist_mat)[0] <= 40 :
					sum_cost = 0		
					for b in dist_mat:
						fila = np.sort(b)
						sum_cost += np.var(fila[1:5])  #sum(abs(fila[1:5] - expected_dif_of_4_closest) )
					#print(sum_cost)

					if sum_cost < 0.1 and  np.shape(dist_mat)[0]<= 50 and (not np.shape(dist_mat)[0] in ya_estan):
						print("encontre una de " + str(np.shape(dist_mat)[0]))
						print(i)
						found.append([lev, np.shape(dist_mat)[0], parameters, n_row_fst ])
						ya_estan.append(np.shape(dist_mat)[0])
	return found


busquedas = [
			[39, [1,13], [1,13]],
			[61, [1,13], [1,13]],
			[60, [1,13], [1,13]],
#			[40, [1,9], [1,9]], # primero n_row_fst, luego level, luego position
#			  [40, [9,18], [1,9]],
#			  [40, [1,9], [9,16]],
#			  [41, [1,9], [1,9]],
#			  [41, [1,9], [9,16]] 
]


found=[]
for caso in busquedas:
	esta_vez = search_inside_limits(caso[0], ws=ws, level= caso[1], position=caso[2], ya_estan=ya_estan )
	found=found+esta_vez
	ya_estan.append([i[0] for i in esta_vez])


rows = []
for n in found:
	this_scale = n[2][3]
	best_sofar = 20000
	this_cords = []
	x = np.linspace(-30, 30, n[3] + 1)
	y = np.linspace(-30, 30,  n[3] + 1)

	X, Y = np.meshgrid(x[:-1], y[:-1])

	Xf = X.flatten(order='A')
	Yf = Y.flatten(order='A')

	cordinates = np.array([Xf, Yf, [1]*len(Xf) ])


	for i in np.linspace(0.9 * this_scale, 1.1 * this_scale, 10000):#poner 8000
		parameters = [n[2][0], 0, 0 ,  i]
		transformed = transform_coordinates(cordinates=cordinates, parameters=parameters, lvl=n[0])
		transformed = transformed[np.all(transformed[:,0:2] >= 0, axis=1) , :]
		transformed = transformed[np.all(transformed[:,0:2] < ws, axis=1 ) , :]
		dist_mat = distances_torus(transformed[:,0], transformed[:,1], ws)
		sum_cost = 0	
		for i in dist_mat:
			fila = np.sort(i)
			sum_cost += np.var(fila[1:5]) 
				#print(sum_cost)
		if sum_cost < best_sofar and np.shape(dist_mat)[0]==n[1] :
			best_sofar = sum_cost
			this_cords = transformed
	
	rows.append([n[1]]+[str(i[0])+":"+str(i[1]) for i in this_cords])

writer = csv.writer(open("configurations/point_configurations3.csv", mode="w", newline=""))
writer.writerows(rows)
print("fin")