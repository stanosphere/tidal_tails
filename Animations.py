from numpy import array, sqrt, sin, cos, pi, arctan
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
from All_Particles_Calculations import Ring_with_perturbing_galaxy
import matplotlib.animation as animation

h = 0 #h=0 for point mass @ centre, h=1 for dark matter halo.
G = 1.0 #Newton's Constant
M = 1.0 #Mass of central object (need to change to 6 in halo case)
M_p = 1.0 #Mass of 2nd galaxy
phi = pi/2 #This is angle of orbit of 2nd galaxy to the plane of the disk

f = 20 #'focal lengh' (distance of closst approach in case where M_p is small)
y_p = 40 #initial y displacement of perturbing galaxy
x_p = (f - (y_p**2)/(4*f))*cos(phi)
z_p = (f - (y_p**2)/(4*f))*sin(phi)
r_p = sqrt(x_p**2 + y_p**2 + z_p**2)

v_mag_p = sqrt(2*M*G/r_p) #Magnitude of initial velocity of 2nd Galaxy
theta = arctan(2*f/y_p) #Direction chosen here to match parabolic orbit 
#Note: theta and phi are NOT conventional polar coordinates

N = array([12,18,24,30,36])*3 #Choose number of particles in each ring
R = array([2.0,3.0,4.0,5.0,6.0]) #Choose Radius of each ring	
N_tot = np.sum(N) + 10 #Both galactic centres are simulated 5 times each

ring_0 = Ring_with_perturbing_galaxy(N[0],R[0],f,y_p,phi,v_mag_p,theta,M_p,h)	
ring_1 = Ring_with_perturbing_galaxy(N[1],R[1],f,y_p,phi,v_mag_p,theta,M_p,h)
ring_2 = Ring_with_perturbing_galaxy(N[2],R[2],f,y_p,phi,v_mag_p,theta,M_p,h)
ring_3 = Ring_with_perturbing_galaxy(N[3],R[3],f,y_p,phi,v_mag_p,theta,M_p,h)
ring_4 = Ring_with_perturbing_galaxy(N[4],R[4],f,y_p,phi,v_mag_p,theta,M_p,h)
all_particles = np.hstack((ring_0,ring_1,ring_2,ring_3,ring_4))

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

d_min = min(galac_sep)
print('The Distance of closest approach is: %s' %d_min)
tail_lengh = np.amax(dis)
print('The Maximum Tail Lengh is: %s' %tail_lengh)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, autoscale_on=False, xlim=(-40, 40), ylim=(-40, 40),
xlabel = 'x displacement',ylabel = 'y displacement')
ax1.grid()
ax1.plot(0,0,marker='o', color='r',markersize = 3)

line1, = ax1.plot([], [], 'o',markersize = 3)
time_template = 'time = %.1f'
time_text1 = ax1.text(0.05, 0.9, '', transform=ax1.transAxes)

def init1():
    line1.set_data([], [])
    time_text1.set_text('')
    return line1, time_text1

def animate1(i):
	x_plot = np.zeros(N_tot)
	y_plot = np.zeros(N_tot)
	for j in range(0,N_tot):
		x_plot[j] = all_particles[i,6*j]
		y_plot[j] = all_particles[i,6*j+2]

	line1.set_data(x_plot, y_plot)
	time_text1.set_text(time_template%(i*0.3))
	return line1, time_text1

ani1 = animation.FuncAnimation(fig1, animate1, np.arange(1, len(all_particles)),
    interval=25, blit=True, init_func=init1)
	
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, autoscale_on=False, xlim=(-35, 35), ylim=(-10, 60),
xlabel = 'y displacement',ylabel = 'z displacement')
ax2.grid()
ax2.plot(0,0,marker='o', color='r',markersize = 3)

line2, = ax2.plot([], [], 'o',markersize = 3)
time_template = 'time = %.1f'
time_text = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)
	
def init2():
    line2.set_data([], [])
    time_text.set_text('')
    return line2, time_text	
	
def animate2(i):
	y_plot = np.zeros(N_tot)
	z_plot = np.zeros(N_tot)
	for j in range(0,N_tot):
		y_plot[j] = all_particles[i,6*j+2]
		z_plot[j] = all_particles[i,6*j+4]

	line2.set_data(y_plot, z_plot)
	time_text.set_text(time_template%(i*0.3))
	return line2, time_text	
	
ani2 = animation.FuncAnimation(fig2, animate2, np.arange(1, len(all_particles)),
    interval=25, blit=True, init_func=init2)
	
plt.show()