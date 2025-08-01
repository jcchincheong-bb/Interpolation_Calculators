# -*- coding: utf-8 -*-
"""
Script to test general polynomial interpolation of 3 data points

Created on Tue Apr 22 22:58:39 2025

@author: BeanieBabyJ
"""

import numpy as np
import matplotlib.pyplot as plt

# Function
def f(x):
    #return 4 * (1 - 2 * np.cos(x)) * np.exp(-0.36 * x) + 7 + np.sin(2 * x)  # some crazy fucking function
    return (1 / 100)  * (x ** 4 - 32)

# Parameters
a = 5   # Starting value
h = 2  # Step size

# Calculate function values
f_a_minus_h = f(a - h)
f_a = f(a)
f_a_plus_h = f(a + h)

# Coefficients formuale
c0 = ((a**2 + a*h) * f_a_minus_h + 2 * (h**2 - a**2) * f_a + a * (a - h) * f_a_plus_h) / (2 * h**2)
c1 = (-(2*a + h) * f_a_minus_h + 4*a * f_a + (h - 2*a) * f_a_plus_h) / (2 * h**2)
c2 = (f_a_plus_h - 2*f_a + f_a_minus_h) / (2 * h**2)

# Polynomial function
def p(x):
    return c0 + c1 * x + c2 * x**2

# Plotting
x_vals = np.linspace(a - 2*h, a + 2*h, 400)  # We double the step here so we actually get a nice picture 
y_vals = p(x_vals)
f_vals = f(x_vals)

plt.figure(figsize=(8, 5))
plt.plot(x_vals, y_vals, label='Polynomial p(x)', color='blue')
plt.plot(x_vals, f_vals, label='Function f(x)', color='black')
plt.plot([a-h, a, a+h], [f_a_minus_h, f_a, f_a_plus_h], 'ro', label='f(x) sample points')
plt.title('Interpolating Polynomial')
plt.xlabel('x')
plt.ylabel('p(x)')
plt.legend()
plt.grid(True)
plt.show()
