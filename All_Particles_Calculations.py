from numpy import array, sqrt, sin, cos, pi, arctan
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import scipy.integrate as integrate

#creates a time array 
dt = 0.3
T = 700 #upper time (usually 700)
t = np.arange(0.0, T, dt)	

G = 1.0 #Newton's Constant
M = 1.0 #Mass of central object (need to change to 6 in halo case)
R_halo = 7 #Just bigger than outermost ring

def Ring_with_perturbing_galaxy(N,R,f,y_p,phi,v_mag_p,theta,M_p,h):
	
	def derivs(state, t):
	
		dydx = np.zeros_like(state)
		# dydx is the list of derivatives w.r.t. time
		for j in range(0,N):
			i = 6*j
			dydx[i] = state[i+1]
			dydx[i+2] = state[i+3]
			dydx[i+4] = state[i+5]
			q0 = sqrt((state[i]-state[6*(N+1)])**2 + (state[i+2]-state[6*(N+1)+2])**2 + (state[i+4]-state[6*(N+1)+4])**2)
			# q0 is distance from galactic centre of original galaxy
			q1 = sqrt((state[i]-state[6*N])**2 + (state[i+2]-state[6*N+2])**2 + (state[i+4]-state[6*N+4])**2)
			# q1 is distance to centre of perturbing galaxy
			
			if h==0: #point mass at centre of galaxy
				dydx[i+1] = -G*M*(state[i]-state[6*(N+1)])/(q0**3) -G*M_p*(state[i]-state[6*N])/(q1**3)
				dydx[i+3] = -G*M*(state[i+2]-state[6*(N+1)+2])/(q0**3) -G*M_p*(state[i+2]-state[6*N+2])/(q1**3)
				dydx[i+5] = -G*M*(state[i+4]-state[6*(N+1)+4])/(q0**3) -G*M_p*(state[i+4]-state[6*N+4])/(q1**3)
				#The above sorts out gravitational acceleration due to the 2 galactic centres
			
			if h==1: #Dark matter Halo
				if q0 > R_halo and q1 > R_halo:
					dydx[i+1] = -G*M*(state[i]-state[6*(N+1)])/(q0**3) -G*M_p*(state[i]-state[6*N])/(q1**3)
					dydx[i+3] = -G*M*(state[i+2]-state[6*(N+1)+2])/(q0**3) -G*M_p*(state[i+2]-state[6*N+2])/(q1**3)
					dydx[i+5] = -G*M*(state[i+4]-state[6*(N+1)+4])/(q0**3) -G*M_p*(state[i+4]-state[6*N+4])/(q1**3)
				elif q0 <= R_halo and q1 > R_halo:
					dydx[i+1] = -G*M*(state[i]-state[6*(N+1)])/(R_halo**3) -G*M_p*(state[i]-state[6*N])/(q1**3)
					dydx[i+3] = -G*M*(state[i+2]-state[6*(N+1)+2])/(R_halo**3) -G*M_p*(state[i+2]-state[6*N+2])/(q1**3)
					dydx[i+5] = -G*M*(state[i+4]-state[6*(N+1)+4])/(R_halo**3) -G*M_p*(state[i+4]-state[6*N+4])/(q1**3)
				elif q0 > R_halo and q1 <= R_halo:	
					dydx[i+1] = -G*M*(state[i]-state[6*(N+1)])/(q0**3) -G*M_p*(state[i]-state[6*N])/(R_halo**3)
					dydx[i+3] = -G*M*(state[i+2]-state[6*(N+1)+2])/(q0**3) -G*M_p*(state[i+2]-state[6*N+2])/(R_halo**3)
					dydx[i+5] = -G*M*(state[i+4]-state[6*(N+1)+4])/(q0**3) -G*M_p*(state[i+4]-state[6*N+4])/(R_halo**3)
				elif q0 <= R_halo and q1 <= R_halo:
					dydx[i+1] = -G*M*(state[i]-state[6*(N+1)])/(R_halo**3) -G*M_p*(state[i]-state[6*N])/(R_halo**3)
					dydx[i+3] = -G*M*(state[i+2]-state[6*(N+1)+2])/(R_halo**3) -G*M_p*(state[i+2]-state[6*N+2])/(R_halo**3)
					dydx[i+5] = -G*M*(state[i+4]-state[6*(N+1)+4])/(R_halo**3) -G*M_p*(state[i+4]-state[6*N+4])/(R_halo**3)
				#The above sorts out gravitational acceleration due to the 2 galactic centres
		
		d = sqrt((state[6*N]-state[6*(N+1)])**2 + (state[6*N+2]-state[6*(N+1)+2])**2 + (state[6*N+4]-state[6*(N+1)+4])**2)
		#d is distance between galactic centres
		
		#This bit sorts out the 1st galaxy's centre
		dydx[6*(N+1)] = state[6*(N+1)+1]
		dydx[6*(N+1)+2] = state[6*(N+1)+3]
		dydx[6*(N+1)+4] = state[6*(N+1)+5]
		
		dydx[6*(N+1)+1] = -G*M_p*(state[6*(N+1)]-state[6*N])/(d**3)
		dydx[6*(N+1)+3] = -G*M_p*(state[6*(N+1)+2]-state[6*N+2])/(d**3)
		dydx[6*(N+1)+5] = -G*M_p*(state[6*(N+1)+4]-state[6*N+4])/(d**3)
		
		#This bit sorts out the 2nd galaxy's centre
		dydx[6*N] = state[6*N+1]
		dydx[6*N+2] = state[6*N+3]
		dydx[6*N+4] = state[6*N+5]
		
		dydx[6*N+1] = -G*M*(state[6*N]-state[6*(N+1)])/(d**3)
		dydx[6*N+3] = -G*M*(state[6*N+2]-state[6*(N+1)+2])/(d**3)
		dydx[6*N+5] = -G*M*(state[6*N+4]-state[6*(N+1)+4])/(d**3)
		return dydx
	
	if h==0:
		v_mag = sqrt(M*G/R) #This is required speed for circular orbit for point mass centre
	if h ==1:
		v_mag = sqrt(M*G*(R**2)/R_halo**3) #This is required speed for circular orbit in DM Halo
	x=np.zeros(N)
	y=np.zeros(N)
	z=np.zeros(N)
	vx=np.zeros(N)
	vy=np.zeros(N)
	vz=np.zeros(N)
	ran = rand.random(N)*2*pi
	for i in range(0,N):
		x[i] = R*cos(ran[i])
		y[i] = R*sin(ran[i])
		vx[i] = -v_mag*sin(ran[i])
		vy[i] = v_mag*cos(ran[i])
		#keep z and vz arrays 0 to keep disk in z=0 plane

	x_p = (f - (y_p**2)/(4*f))*cos(phi)
	z_p = (f - (y_p**2)/(4*f))*sin(phi)	
	vx_p = v_mag_p*(cos(theta))*(cos(phi))
	vy_p = -v_mag_p*sin(theta)
	vz_p = v_mag_p*(cos(theta))*(sin(phi))

	#This just creates the initial conditions in a convenient format
	state_original_galaxy = np.zeros(6) #original galaxy centre is at rest at origin.
	state_test_particles = np.reshape(np.transpose(np.vstack((x,vx,y,vy,z,vz))),6*N)
	state_pert = np.array([x_p,vx_p,y_p,vy_p,z_p,vz_p]) #2nd galaxy's starting values
	state = np.hstack((state_test_particles,state_pert,state_original_galaxy))
	#Last 6 entries in state are for original galaxy's centre
	#The 6 before them are for the perturbing galaxy
	#This bit does the integration
	r = integrate.odeint(derivs, state, t)
	return r



