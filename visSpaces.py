import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

def update(val):
	# Update Resolution from slider
	Res = sRes.val
	# Clear and repopulate axes
	axes.clear()
	draw()	
	
	# Draw new resolution grid
	for r in range(0,int(Res)+1):
		xline = [-R,R]
		yline = [-R + r*2*R/(Res),-R + r*2*R/(Res)]
		axes.plot(xline, yline, color='black', linestyle='--')
		yline = [-R,R]
		xline = [-R + r*2*R/(Res),-R + r*2*R/(Res)]
		axes.plot(xline, yline, color='black', linestyle='--')

	fig.canvas.draw_idle()

def draw():
	# Bounding Surface
	axes.plot(xx,yy, color='black',label='Bounding Surface', linewidth=3)

	# Principle stress axis
	axes.arrow(0,0,0,1.2*R,width=.003*R,head_width=.005,head_length=.005, color='gray',linestyle=':')
	axes.arrow(0,0,1.4*np.sqrt(2/3)*R,-np.sqrt(1/3)*R,width=.003*R,head_width=.005,head_length=.005, color='gray',linestyle='--')
	axes.arrow(0,0,1.4*-np.sqrt(2/3)*R,-np.sqrt(1/3)*R,width=.003*R,head_width=.005,head_length=.005, color='gray',linestyle='--')
	axes.text(0, 1.2*R, '$\sigma_1\'$', color = 'gray')
	axes.text(1.2*np.sqrt(2/3)*R, -np.sqrt(1/3)*R, '$\sigma_2\'$', color = 'gray')
	axes.text(1.2*-np.sqrt(2/3)*R, -np.sqrt(1/3)*R, '$\sigma_3\'$', color = 'gray')

	# Orgin \Pi Axis
	axes.text(0,0,'O$\Pi$',fontsize=18)
	axes.arrow(0,0,R/2,0,width=.005*R,head_width=.05*R,head_length=.05*R, color='black')
	axes.arrow(0,0,0,R/2,width=.005*R,head_width=.05*R,head_length=.05*R, color='black')
	axes.text(0, R/2 + 0.05*R, '$x$', color = 'black')
	axes.text(R/2 + 0.05*R, 0, '$y$', color = 'black')
	
	# Array \Pi Axis
	axes.text(-R,R, 'A$\Pi$', fontsize=18)
	axes.arrow(-R,R,R/2,0,width=.005*R,head_width=.05*R,head_length=.05*R, color='black')
	axes.arrow(-R,R,0,-R/2,width=.005*R,head_width=.05*R,head_length=.05*R, color='black')
	axes.text(-R/2 + 0.05*R, R, '$i$', color = 'black')
	axes.text(-R - 0.05*R, R/2 - 0.1*R, '$j$', color = 'black')

	# Axis settings
	axes.set_xlim(-(5/4)*R,(5/4)*R)
	axes.set_xticks([-R,0,R])
	axes.set_xticklabels(['R','0','R'])
	axes.set_ylim(-(5/4)*R,(5/4)*R)
	axes.set_yticks([-R,0,R])
	axes.set_yticklabels(['R','0','R'])
	axes.grid(False)

	plt.show()


# ---------------------------------------------------------------------

R = 0.1 # Initial bounding surface radius
Res0 = 10 # Initial Resolution

# Set up figure
fig, axes = plt.subplots(1,1)

# Bounding surface coordinates
rr = np.linspace(0,2*np.pi,360)
xx = R*np.sin(rr)
yy = R*np.cos(rr)

# Set up resolution slider
axcolor = 'white'
axRes = plt.axes([0.2, 0.9, 0.65, 0.03], facecolor=axcolor)
sRes = Slider(axRes, 'Resolution', 1, 100, valinit=Res0, valstep=1)
sRes.on_changed(update)

# Draw initial resolution grid
for r in range(0,Res0+1):
	xline = [-R,R]
	yline = [-R + r*2*R/(Res0),-R + r*2*R/(Res0)]
	axes.plot(xline, yline, color='black', linestyle='--')
	yline = [-R,R]
	xline = [-R + r*2*R/(Res0),-R + r*2*R/(Res0)]
	axes.plot(xline, yline, color='black', linestyle='--')

# Draw axes
draw()

plt.show()

