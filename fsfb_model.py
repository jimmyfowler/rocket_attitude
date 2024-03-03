import numpy as np
import matplotlib.pyplot as plt

# turn notebook into functions & main()


# Moment of Inertia
I = 10000

# FSFB GAINS
k1 = 1000
k2 = 0.05

# INITIAL CONDITION
theta_0 = np.pi/6 # rad
omega_0 = np.pi/8 # rad/sec

# define / initialize error
# this is not used to 
err = np.sqrt(theta_0**2 + omega_0**2)

# initialize time and interval tracker
t = np.array([0])
dt = 0.1
i = 0

# initialize the states
theta = np.array([theta_0])
omega = np.array([omega_0])


# main()

while err > 0.001:

    # calculate next state
    next_theta = theta[i] + omega[i]*dt
    next_omega = omega[i] + ( (-k1/I)*theta[i] - (k2*omega[i]) )

    # add next state to array
    theta = np.append(theta, next_theta)
    omega = np.append(omega, next_omega)

    # error = magnitude of current state
    err = np.sqrt(theta[i]**2 + omega[i]**2)
    
    # next interval (i=0 is initial condition)
    t = np.append(t, t[i]+dt)
    i += 1

    if t[i] > 10000:
        print("It took too long!")
        break

print(f"Time to reach goal: {round(t[i], 4)} seconds")
# print(f"Actuator effort spent {effort} Nm")

# convert RAD & RAD/S --> DEG & DEG/S for plots
theta_deg = theta*360/(2*np.pi)
omega_deg_s = omega*360/(2*np.pi)


# plot states over time
fig, ax = plt.subplots(2, 1, figsize=[8,8])

for axis in ax:
    axis.set_ylim([-80, 80])

ax[0].plot(t, theta_deg, 'r')
ax[0].set_title("States over Time", fontsize=18)
ax[0].set_ylabel(r"$\theta$ [deg]", fontsize=15)

ax[1].plot(t, omega_deg_s, 'm')
ax[1].set_title
ax[1].set_ylabel(r"$\omega$ [deg/s]", fontsize=15)
ax[1].set_xlabel("Time [sec]", fontsize=15)

fig.savefig("figs/states_time")
fig.show()

# phase plot
fig2, ax2 = plt.subplots(1)

ax2.set_title("Phase Portrait", fontsize=15)
ax2.plot(theta_deg, omega_deg_s, color='mediumblue')
ax2.plot(theta_deg[0], omega_deg_s[0], 'g.', ms=12) # plot starting point

# Move left y-axis and bottom x-axis to centre, passing through (0,0)
ax2.spines['top'].set_position('center')
ax2.spines['right'].set_position('center')

# # Eliminate upper and right axes
# ax2.spines['top'].set_color('none')
# ax2.spines['right'].set_color('none')

# # Show ticks in the left and lower axes only
# ax2.xaxis.set_ticks_position('bottom')
# ax2.yaxis.set_ticks_position('left')

ax2.set_xlim([-45, 45])
ax2.set_ylim([-45, 45])

ax2.set_xlabel(r"$\theta$ [deg]", fontsize=12)
ax2.set_ylabel(r"$\omega$ [deg/s]", fontsize=12)

fig2.savefig("figs/phase_portrait")
fig2.show()
