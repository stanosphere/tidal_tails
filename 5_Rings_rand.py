from numpy import array, sqrt, sin, cos, pi
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

#This function sets up a single ring of particles
def ring(Ni,Ri):
	
	N = Ni
	R = Ri
	v_mag = sqrt(M*G/R) #This is required speed for circular orbit
	
	def derivs(state, t):

		dydx = np.zeros_like(state)
		# dydx is the list of derivatives w.r.t. time
		for j in range(0,N):
			i = 6*j
			dydx[i] = state[i+1]
			dydx[i+2] = state[i+3]
			dydx[i+4] = state[i+5]
			r = sqrt((state[i])**2 + (state[i+2])**2 + (state[i+4])**2)
			dydx[i+1] = -G*M*(state[i])/(r**3)
			dydx[i+3] = -G*M*(state[i+2])/(r**3)
			dydx[i+5] = -G*M*(state[i+4])/(r**3)
		return dydx

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
		#keep z and vz arrays 0 to keep ploblem in z=0 plane

	#This just creates the initial conditions in a convenient format
	state=np.reshape(np.transpose(np.vstack((x,vx,y,vy,z,vz))),6*N)

	#This bit does the integration
	r = integrate.odeint(derivs, state, t)
	return r

G = 1.0 #Newton's Constant
M = 1.0 #Mass of central object

#creates a time array 
dt = 0.1
t = np.arange(0.0, 200, dt)	
	
N = array([24,36,48,60,72]) #Choose number of particles in each ring
R = array([2.0,3.0,4.0,5.0,6.0]) #Choose Radius of each ring	
N_tot = np.sum(N)

#It didn't like me doing a loop here :(
r0 = ring(N[0],R[0])
r1 = ring(N[1],R[1])
r2 = ring(N[2],R[2])
r3 = ring(N[3],R[3])
r4 = ring(N[4],R[4])
r = np.hstack((r0,r1,r2,r3,r4))

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-6.5, 6.5), 
ylim=(-6.5, 6.5),xlabel = 'x displacement',ylabel = 'y displacement')
ax.grid()
ax.plot(0,0,marker='o', color='r',markersize = 15)

line, = ax.plot([], [], 'o')
time_template = 'time = %.1f'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(i):
	x_plot = np.zeros(N_tot)
	y_plot = np.zeros(N_tot)
	for j in range(0,N_tot):
		x_plot[j] = r[i,6*j]
		y_plot[j] = r[i,6*j+2]

	line.set_data(x_plot, y_plot)
	time_text.set_text(time_template%(i*dt))
	return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(r)),
    interval=25, blit=True, init_func=init)
	
plt.show()
