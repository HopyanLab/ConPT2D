#!/usr/bin/env /usr/bin/python3
import numpy as np
from matplotlib import pyplot as plt

n = 8
theta = 2*np.pi/n
extra = 0.1
x = np.linspace(np.cos(theta) - extra*(1-np.cos(theta)),
				1 + extra*(1-np.cos(theta)),10)
m = -np.sin(theta)/(1-np.cos(theta))
b = -m
y = m*x + b
fig = plt.figure()
ax = plt.axes()
for i in range(n):
	phi = i*theta
	ax.plot(x*np.cos(phi)-y*np.sin(phi),x*np.sin(phi)+y*np.cos(phi))
ax.set_box_aspect(1)
plt.show()
