import pygame
import numpy as np
import scipy.integrate as integrate
import cmath

# Constants
h_ = 1
a = 1
m = 1

def f(x):
    return np.exp((-50) * ((x - 0.5) ** 2))

def eigen_function(n_i, x):
    return np.sqrt(2 / a) * np.sin(n_i * np.pi * x / a)

def eigen_value(n_i):
    return (h_ * np.pi * n_i / a) ** 2 / (2 * m)

def normalization_coef(f, min, max):
    N, _ = integrate.quad(f, min, max)
    return 1 / N

def scalar_product(f, n_i, min, max):
    integrand = lambda x: normalization_coef(f, min, max) * f(x) * eigen_function(n_i, x)
    I, _ = integrate.quad(integrand, min, max)
    return I

cn = [scalar_product(f, n, 0, a) for n in range(1, 101)]
A = normalization_coef(f, 0, a)

t = 0
xstep = 0.001
tstep = 0.01
passo = np.arange(0, a, 0.001)
maxf = 0

for x in passo:
    if maxf < f(x) :
        maxf = f(x)
maxf += 1

# pygame setup
pygame.init()
screen = pygame.display.set_mode((800, 800))
pygame.display.set_caption("Wave Function Evolution")

running = True  # Variable to control the main loop

while running:  # Main loop
    
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False  # If the user closes the window, exit the loop
    
    if running:
        screen.fill((0, 0, 0))
        psi_2 = []
        passo_psi_2 = []

        # Axis
        pygame.draw.lines(screen, "white", True, [[100, 150], [700, 150], [700, 650], [100, 650]], 5)

        # Computation of scalar product in order to calculate the new shape of wave function
        for x in passo:
            psi = complex(0., 0.)

            for n in range(len(cn)):
                psi += cmath.rect(cn[n] * eigen_function(n + 1, x), (-1) * eigen_value(n + 1) * t / h_) 

            psi_2.append(abs(psi) ** 2)
            passo_psi_2.append(x)
        
        # Visualisation of the result
        for i in range(len(psi_2) - 1): 
            x1 = int(100 + (600 / a) * passo_psi_2[i])
            y1 = int(650 - (50 / maxf) * psi_2[i])
            x2 = int(100 + (600 / a) * passo_psi_2[i + 1])
            y2 = int(650 - (50 / maxf) * psi_2[i + 1])
            pygame.draw.line(screen, "red", (x1, y1), (x2, y2), 2)

        t += tstep
        # Update the simulation
        pygame.display.flip()

pygame.quit()
