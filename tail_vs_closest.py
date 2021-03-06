from numpy import array, sqrt, sin, cos, pi, arctan
import numpy as np
import matplotlib.pyplot as plt
from All_Particles_Calculations import Ring_with_perturbing_galaxy

d_min = np.zeros(100)
tail_lengh = np.zeros(100)
f_array = np.zeros(100)

for a in range (0,100):
	h = 0 #h=0 for point mass @ centre, h=1 for dark matter halo.
	f = 8*a/100 + 22 #'focal lengh' (distance of closst approach in case where M_p is small)
	# keeps f in range 14 to 22
	G = 1.0 #Newton's Constant
	M = 1.0 #Mass of central object
	M_p = 1.0#Mass of 2nd galaxy
	phi = 0 #This is angle of orbit of 2nd galaxy to the plane of the disk

	y_p = 40 #initial y displacement of perturbing galaxy
	x_p = (f - (y_p**2)/(4*f))*cos(phi)
	z_p = (f - (y_p**2)/(4*f))*sin(phi)
	r_p = sqrt(x_p**2 + y_p**2 + z_p**2)

	v_mag_p = sqrt(2*M*G/r_p) #Magnitude of initial velocity of 2nd Galaxy
	theta = arctan(2*f/y_p) #Direction chosen here to match parabolic orbit 
	#Note: theta and phi are NOT conventional polar coordinates

	N = array([12,18,24,30,36]) #Choose number of particles in each ring
	R = array([2.0,3.0,4.0,5.0,6.0]) #Choose Radius of each ring	
	N_tot = np.sum(N) + 10 #Both galactic centres are simulated 5 times each

	ring_4 = Ring_with_perturbing_galaxy(N[4],R[4],f,y_p,phi,v_mag_p,theta,M_p,h)
	# Particles contributing to furtherst point on tail will be from outermost ring
	
	#I want to find closest distance of galactic centres...
	#And the lengh of the tail formed.
	dim = np.shape(ring_4) #first element gives number of time steps
	n = dim[0]
	galac_sep = np.zeros(n)
	dis = np.zeros((n,N[4]))
	for j in range(0,n):
		galac_sep[j] = sqrt((ring_4[j,6*N[4]]-ring_4[j,6*(N[4]+1)])**2 + (ring_4[j,6*N[4]+2]-ring_4[j,6*(N[4]+1)+2])**2
	+(ring_4[j,6*N[4]+4]-ring_4[j,6*(N[4]+1)+4])**2)
		for k in range(0,N[4]):
			i = 6*k
			dis[j,k] = sqrt((ring_4[j,6*(N[4]+1)]-ring_4[j,i])**2 + (ring_4[j,6*(N[4]+1)+2]-ring_4[j,i+2])**2
	+(ring_4[j,6*(N[4]+1)+4]-ring_4[j,i+4])**2)

	d_min[a] = min(galac_sep)
	tail_lengh[a] = np.amax(dis)
	f_array[a] = f
	
	
fig, ax1 = plt.subplots()
ax1.plot(d_min,tail_lengh)
ax1.set_xlabel("Distance of Closest Approach")
ax1.set_ylabel("Tail Lengh")
ax1.set_title("Maximum Tail Lengh vs Distance of Closest Approach")

print('The Focal length is: %s' %f_array)
print('The Distance of closest approach is: %s' %d_min)
print('The Maximum Tail Lengh is: %s' %tail_lengh)

plt.show()