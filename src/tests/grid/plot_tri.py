import numpy
import matplotlib.pyplot as plt

data = numpy.loadtxt('surf.check')

u = data[:,0].reshape((100,100))
v = data[:,1].reshape((100,100))

x = data[:,2].reshape((100,100))
y = data[:,3].reshape((100,100))
z = data[:,4].reshape((100,100))

r = numpy.sqrt(x*x + y*y + z*z)
err = 1.-r

plt.figure(1)
plt.contour(u,v,x,40)
plt.colorbar()

plt.figure(2)
plt.contour(u,v,y,40)
plt.colorbar()

plt.figure(3)
plt.contour(u,v,z,40)
plt.colorbar()

plt.figure(4)
plt.contour(u,v,err,40)
plt.colorbar()

plt.show()
