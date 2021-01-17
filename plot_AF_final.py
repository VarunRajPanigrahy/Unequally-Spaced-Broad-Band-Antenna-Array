import math
import random
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


pi = math.pi
N = 20
f = 10 * 10**9  # 10 Ghz X band mid frequency.
c = 3 * 10**8
l = c/f # wavelength lambda.
d = 0.50*l


def cosh_inverse(x):
    return math.log(x + math.sqrt(x**2-1))

def dolph_chebyshev_coefficients(N, RoDb):

    
    Ro = pow(10, RoDb/20)

    window = signal.chebwin(N, RoDb)

    coeffs_normalized = window[N//2:]

    zo = math.cosh(cosh_inverse(Ro)/(N-1))

    normalization_factor = (zo**(N-1))/coeffs_normalized[-1]

    coeffs = [normalization_factor*i for i in coeffs_normalized]

    return coeffs

def C(n,r):
    return math.factorial(n)//(math.factorial(r)*math.factorial(n-r))

def binomial_coefficients(N):
    m = N//2
    coeffs = [C(N-1,i+m) for i in range(m)]
    return coeffs

def array_factor2m(n, theta, coefficients):
    m = n//2
    u = (pi * d/l)*(math.cos(theta))

    AF = 0
    for i in range(1,m+1):
        AF += (coefficients[i-1] * math.cos((2*i-1) * u))
        
    
    return AF

plot_number = 1

def get_AF(n,coefficients):

  global plot_number

  rads = np.arange(0, (np.pi), 0.001)
  AF_actual = [array_factor2m(n, rad, coefficients) for rad in rads]
  normal_factor = max(AF_actual)
  
  thetas = [rad*180/pi for rad in rads]
  for i in range(len(AF_actual)):
      AF_actual[i]/=normal_factor
    


  AF_actual = [20*math.log(abs(i),10) for i in AF_actual]
  plt.grid()
  '''
  f1 = plt.figure(plot_number)
  plt.axes(projection = 'polar')
  plt.polar(rads,AF_actual)
  '''
  f2 = plt.figure(plot_number)
  plt.plot(thetas, AF_actual)

  plot_number += 1

N = 50
RoDb = 20

coefficients_predicted = [ 0.4794866442561804,0.47835599038602616,0.4752483313037993,0.4718802105762074,0.4656589321910359,0.460212184488134,0.45091096204707304,0.44374051851260227,0.4314803444962219,0.4229491825718923,0.4074019468296812,0.3981281414039219,0.3794604321049408,0.37100645503890584,0.3476122907107611,0.34220726455288336,0.309097940454831,0.32356756829127564,0.23783247536598448,0,0.20320776515044808,0.2534872930921506,0.20863275804773024,0.20769372527365113,1.088732657736092 ]

coeffs_dolph = dolph_chebyshev_coefficients(N, RoDb)

get_AF(N, coefficients_predicted)
get_AF(N, coeffs_dolph)

plt.show()