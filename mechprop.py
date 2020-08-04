import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit, fsolve
import warnings
from scipy.special import erf, erfc
from scipy.stats import t
import time
from tabulate import tabulate

ID       = 1
files    = str(ID)+'DCP.xls'
path     = '/home/fernando/Papers/Inhomogeneities_peroxide/'
filename = path+files
df       = pd.read_excel(filename)
data     = df.values

test = False

def stressfit(l,d,gc,ge):
    denom  = 1 - (d**2)*((l**2) + 2/l - 3)
    factor = l - 1/(l**2)
    term1 = (1 - d**2)*factor/(denom**2)
    term2 = d**2*factor/denom
    term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
    return gc*(term1 - term2) + 2*(ge/1)*term3

def stressfit_C(l,d,gc,ge,C):
    l += C
    denom  = 1 - (d**2)*((l**2) + 2/l - 3)
    factor = l - 1/(l**2)
    term1 = (1 - d**2)*factor/(denom**2)
    term2 = d**2*factor/denom
    term3 = -l**float(-1 - 1) + l**float(1/2 - 1)
    return gc*(term1 - term2) + 2*(ge/1)*term3

def stressfit_b(l,d,b,gc,ge,C):
    l += C
    denom  = 1 - (d**2)*((l**2) + 2/l - 3)
    factor = l - 1/(l**2)
    term1 = (1 - d**2)*factor/(denom**2)
    term2 = d**2*factor/denom
    term3 = -l**float(-b - 1) + l**float(b/2 - 1)
    return gc*(term1 - term2) + 2*(ge/b)*term3

cols     = len(data[0,:])
strain   = {}
stress   = {}
j = 1
k = 1
for i in range(cols):
    if (i+1) % 2 == 0:
        #y axis
        stress["stress{0}".format(k)] = data[:,i]
        k += 1
        print('Added stress ',i)
    else:
        #x axis
        strain["strain{0}".format(j)] = data[:,i]/100+1
        print('Added strain ',i)
        j += 1

l = 1
for i in range(cols//2):
    plt.plot(strain["strain{0}".format(l)],stress["stress{0}".format(l)],label='Sample {0}'.format(l))
    l+=1
plt.legend()
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Stress-Strain curves')
plt.show()

if test:
    lmin = float(input('Select the minimum value for unitary deformation to avoid regions where precycle effect is more important than elastic behaviour: '))
    lmax = float(input('Select the maximum value for unitary deformation to avoid regions where cristalization induced by deformation are more important than elastic behaviour: '))
else:
    lmin = 1.2
    lmax = 3

l = 1
strainselected = {}
strainaux      = {}
stressselected = {}
stressaux      = {}
warnings.filterwarnings('ignore',message='invalid value encountered in less')
for i in range(cols//2):
    strainaux["strain{0}".format(l)]=strain["strain{0}".format(l)][strain["strain{0}".format(l)]<lmax]
    strainselected["strain{0}".format(l)]=strainaux["strain{0}".format(l)][strainaux["strain{0}".format(l)]>lmin]
    stressaux["stress{0}".format(l)]=stress["stress{0}".format(l)][strain["strain{0}".format(l)]<lmax]
    stressselected["stress{0}".format(l)]=stressaux["stress{0}".format(l)][strainaux["strain{0}".format(l)]>lmin]
    l += 1
    #RunTimeWarningError due to NaN comparing to number, proven that this error may be ignored if appears at this point
print('done')
del strain
del stress
del strainaux
del stressaux

l = 1
d  = []
gc = []
ge = []
C  = []
plt.show()
for i in range(cols//2):
    plt.plot(strainselected["strain{0}".format(l)],stressselected["stress{0}".format(l)],label='Sample {0}'.format(l))
    
    popt, pcov = curve_fit(stressfit,xdata=strainselected["strain{0}".format(l)],ydata=stressselected["stress{0}".format(l)],bounds=((-np.inf,0,0),(np.inf,np.inf,np.inf)))
    #popt, pcov = curve_fit(stressfit_C,xdata=strainselected["strain{0}".format(l)],ydata=stressselected["stress{0}".format(l)],bounds=((-np.inf,0,0,-np.inf),(np.inf,np.inf,np.inf,np.inf)))
    plt.plot(strainselected["strain{0}".format(l)],stressfit(l=strainselected["strain{0}".format(l)],d=popt[0],gc=popt[1],ge=popt[2]),label='HSH fitting for sample {0}'.format(l))
    d.append(popt[0])
    gc.append(popt[1])
    ge.append(popt[2])
    #C.append(popt[3])
    l+=1
plt.legend()
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\sigma$'+' (MPa)')
plt.title('Stress-Strain curves')
plt.show()
print('Gc values for the sample ',str(ID),' :')
print(gc)
print('Ge values for the sample ',str(ID),' :')
print(ge)

def Ac(Nc,Ne):
    a = (3*Ne)/(2*Nc)
    return 0.5 + (np.sqrt(a/np.pi)*np.exp(-a))/erf(np.sqrt(a))
rho = 0.92 #g/cm3
T   = 298 #K
R   = 8.3144 #J/(mol K)
Ms  = 131 #g/mol
Na  = 6.022*np.power(10,23) #Avogadro constant

Ne  = []
Nc  = []
nue = []
nuc = []
for i in range(cols//2):
    Nee  = (1/np.sqrt(6))*rho*R*T/(ge[i]*Ms)
    def func(Nc):
        return Nc - rho*R*T*Ac(Nc,Nee)/(gc[i]*Ms)

    nuee = Na*rho/(Nee*Ms)
    Ne.append(Nee)
    nue.append(nuee)

    Ncc  = fsolve(func,x0=1) #Ncc is an array, so we select the component zero which is a float
    Nc.append(Ncc[0])
    nucc = Na*rho/(Ncc[0]*Ms)
    nuc.append(nucc)

Ne  = np.asarray(Ne)
Nc  = np.asarray(Nc)
nue = np.asarray(nue)
nuc = np.asarray(nuc)

n = Nc*(1-Ac(Nc,Ne))/(2*Ac(Nc,Ne))

phi = np.loadtxt(path+'params.dat',dtype=float)[:,2]
phi = phi[0:3]
z   = np.sqrt(0.8*Ne*(0.25*Nc+0.5*n))/(Ne*(1.513-0.763/np.sqrt(Nc)))*(phi**(-2/3))*np.exp(((0.25*Nc+0.5*n)-0.5043*np.sqrt(0.25*Nc+0.5*n))/(0.8*Ne))*erfc(np.sqrt(((0.25*Nc+0.5*n)-0.5043*np.sqrt(0.25*Nc+0.5*n))/(0.8*Ne)))

Sb = 6/(5*Nc)*z*np.sqrt((Nc+2*n)/(0.8*Ne))*np.arctan(Nc/(2*np.sqrt(Nc*n+n**2)))

gcerr  = t.ppf(0.95,df=len(gc)-1)*np.std(gc)/np.sqrt(len(gc))
geerr  = t.ppf(0.95,df=len(ge)-1)*np.std(ge)/np.sqrt(len(ge))
Ncerr  = t.ppf(0.95,df=len(Nc)-1)*np.std(Nc)/np.sqrt(len(Nc))
Neerr  = t.ppf(0.95,df=len(Ne)-1)*np.std(Ne)/np.sqrt(len(Ne))
nucerr = t.ppf(0.95,df=len(nuc)-1)*np.std(nuc)/np.sqrt(len(nuc))
nueerr = t.ppf(0.95,df=len(nue)-1)*np.std(nue)/np.sqrt(len(nue))
nerr   = t.ppf(0.95,df=len(n)-1)*np.std(n)/np.sqrt(len(n))
zerr   = t.ppf(0.95,df=len(z)-1)*np.std(z)/np.sqrt(len(z))
Sberr  = t.ppf(0.95,df=len(Sb)-1)*np.std(Sb)/np.sqrt(len(Sb))

gc  = np.average(gc)
ge  = np.average(ge)
Nc  = np.average(Nc)
Ne  = np.average(Ne)
nuc = np.average(nuc)
nue = np.average(nue)
n   = np.average(n)
z   = np.average(z)
Sb  = np.average(Sb)

if False:
    print(gc)
    print("Gc = {̣̣0} +/- {̣̣1}".format(gc,gcerr))
    time.sleep(1)
    print('Ge = {̣̣0} +/- {̣̣1}'.format(ge,geerr))
    time.sleep(1)
    print('Nc = {̣̣0} +/- {̣̣1}'.format(Nc,Ncerr))
    time.sleep(1)
    print('Ne = {̣̣0} +/- {̣̣1}'.format(Ne,Neerr))
    time.sleep(1)
    print('nuc = {̣̣0} +/- {̣̣1}'.format(nuc,nucerr))
    time.sleep(1)
    print('nue = {̣̣0} +/- {̣̣1}'.format(nue,nueerr))
    time.sleep(1)
    print('n = {̣̣0} +/- {̣̣1}'.format(n,nerr))
    time.sleep(1)
    print('z = {̣̣0} +/- {̣̣1}'.format(z,zerr))
    time.sleep(1)
    print('Sb = {̣̣0} +/- {̣̣1}'.format(Sb,Sberr))

print(tabulate([['Gc','Ge','Nc','Ne','nuc','nue','n','z','Sb'],np.round([gc,ge,Nc,Ne,nuc,nue,n,z,Sb],decimals=2),np.round([gcerr,geerr,Ncerr,Neerr,nucerr,nueerr,nerr,zerr,Sberr],decimals=2)],tablefmt="latex"))