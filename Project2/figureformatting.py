import numpy as np, matplotlib.pyplot as plt
plt.style.use('ggplot')

font = {'family' : 'serif',
        'weight' : 'normal',
        'style'  : 'normal',
        'size'   : 20}
plt.rc('font', **font)

#plt.rc('font',**{'family':'sans-serif', 'sans-serif':['Helvetica']})
#plt.rc('text', usetex = True)
x = np.linspace(0,10,100)

def f(x):
  y = np.cos(x) + np.sin(x)
  return y


plt.title('Results')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.plot(f(x), 'r', linewidth = 3, label = 'f(x)')
plt.legend()
plt.show()
