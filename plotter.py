import numpy as np 
import matplotlib.pyplot as plt

def Plot_1():
	a = np.arange(0,8.1,0.1)
	p = 40000
	w = 0.4e4
	l = 8
	E = 2e5
	I = 10**-3
	k = np.sqrt(p/(E*I))
	u = k*l/2
	m_arr = []
	defl_arr = []
	x_arr = []
	for x in a:
		x_arr.append(x)
		moment = ((w*l**2)/(4*u**2))*(np.tan(u)*np.sin(2*u*x/l) + np.cos(2*u*x/l) -1)
		deflection = ((w*l**4)/(16*E*I*u**4))*(np.tan(u)*np.sin(2*u*x/l) + np.cos(2*u*x/l) - 1) - (w*l**2/(8*E*I*u**2))*(x*(l-x))
		m_arr.append(moment)
		defl_arr.append(deflection)
	#plt.plot(x_arr,defl_arr)
	plt.plot(x_arr,defl_arr)
	plt.xlabel('x')
	plt.ylabel('deflection')
	plt.title('deflection vs length()')
	plt.show()

def Plot_2():
	a = np.arange(0,3,0.01)
	p = 78000
	q = 10e4
	l = 3
	E = 2e5
	I = 10**-3
	k = np.sqrt(p/(E*I))
	u = k*l/2
	m_arr = []
	defl_arr = []
	x_arr = []
	for x in a:
		x_arr.append(x)
		#moment = 
		#if(x<l/2):
		deflection = (q/(E*I*k**3))*(np.sin(k*l/2))*np.sin(k*x)/(np.sin(k*l)) - q*k*x/(2*E*I*(k**2))
		#else:
		#	deflection = 
		#m_arr.append(moment)
		defl_arr.append(deflection)
	#plt.plot(x_arr,defl_arr)
	plt.plot(x_arr,defl_arr)
	plt.show()

def Plot_3():
	a = np.arange(0,3,0.1)
	p = 78000
	Ma = -5e3
	Mb = -5e3
	l = 3
	E = 2e5
	I = 10**-3
	k = np.sqrt(p/(E*I))
	u = k*l/2
	m_arr = []
	defl_arr = []
	x_arr = []
	for x in a:
		x_arr.append(x)
		moment =  (Ma*np.cos(k*l) +Mb)*(np.sin(kx))/(E*I*np.sin(k*l)*k**3) + Ma*np.cos(k*x)/(E*I*k**2)
		deflection = -1*(Ma*np.cos(k*l) +Mb)*(np.sin(kx))/(E*I*np.sin(k*l)*k**3) - Ma*np.cos(k*x)/(E*I*k**2)

		#m_arr.append(moment)
		defl_arr.append(deflection)
	#plt.plot(x_arr,defl_arr)
	plt.plot(x_arr,defl_arr)
	plt.xlabel('x from 0 - l')
	plt.xlabel('deflection')
	plt.title('deflection vs length')
	plt.show()

def Plot_4():
	a = np.arange(0,3,0.1)
	p = 78000
	q = 10e4
	l = 3
	E = 2e5
	I = 10**-3
	k = np.sqrt(p/(E*I))
	u = k*l/2
	m_arr = []
	defl_arr = []
	x_arr = []
	for x in a:
		x_arr.append(x)
		#moment =  (Ma*np.cos(k*l) +Mb)*(np.sin(kx))/(E*I*np.sin(k*l)*k**3) + Ma*np.cos(k*x)/(E*I*k**2)
		deflection = (q/(2*E*I*k**3))*(np.sin(k*x) + np.cos(k*l/2 -1)*(np.cos(k*x)/np.sin(k*x)) + ((1- np.cos(k*l/2))/np.sin(k*l/2)) -k*x)
		#m_arr.append(moment)
		defl_arr.append(deflection)
	#plt.plot(x_arr,defl_arr)
	plt.plot(x_arr,defl_arr)
	plt.show()



Plot_2()
Plot_1()