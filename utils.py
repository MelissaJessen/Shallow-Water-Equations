import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
from PIL import Image

def minmod(a, b):
    if a * b <= 0:
        return 0
    else:
        return min(abs(a), abs(b), 1) * np.sign(a)

def superbee(a, b):
    if a * b <= 0:
        return 0
    else:
        return max(0, min(2 * abs(a), abs(b)), min(abs(a), 2 * abs(b)))

def van_leer(a, b):
    if a * b <= 0:
        return 0
    else:
        return (a * abs(b) + b * abs(a)) / (abs(a) + abs(b))
    
def no_limiter(a, b):
    return 0


def roe_solver(hL, huL, hR, huR, g=9.81, lim_func = minmod):
    # Calculate Roe averages
    h_tilde = (hL + hR) / 2
    u_tilde = (np.sqrt(hL) * huL / hL + np.sqrt(hR) * huR / hR) / (np.sqrt(hL) + np.sqrt(hR))

    c = np.sqrt(g * h_tilde)

    # Calculate eigenvalues
    lambda1 = u_tilde - c
    lambda2 = u_tilde + c

    # Compute fluxes
    F_L = np.array([huL, huL**2 / hL + 0.5 * g * hL**2])
    F_R = np.array([huR, huR**2 / hR + 0.5 * g * hR**2])

     # Difference 
    dU = np.array([hR - hL, huR - huL])

    # Compute Roe flux
    if lambda1 >= 0 and lambda2 >= 0:
        return F_L
    elif lambda1 < 0 and lambda2 > 0:
        # Apply flux limiter (minmod)
        dF = F_R - F_L
        flux = 0.5 * (F_L + F_R) - 0.5 * c * dU
        limiter = lim_func(dF[0], c * dU[0])
        flux[0] -= limiter
        return flux
    else:
        return F_R
    
def compute_interface_fluxes(h, hu, N, g, lim_func = minmod):
    fluxes = np.zeros((2, N+1))
    
    # Apply boundary conditions at the walls
    hL, huL = h[0], -hu[0]  # Reflective boundary at x=0
    hR, huR = h[1], hu[1]
    fluxes[:, 0] = roe_solver(hL, huL, hR, huR, g, lim_func)
    
    for i in range(1, N):
        fluxes[:, i] = roe_solver(h[i-1], hu[i-1], h[i], hu[i], g, lim_func)
        
    hL, huL = h[N-2], hu[N-2]
    hR, huR = h[N-1], -hu[N-1]  # Reflective boundary at x=L
    fluxes[:, N] = roe_solver(hL, huL, hR, huR, g, lim_func)
    
    return fluxes

def update_variables(h, hu, dt, fluxes, N, dx):
    h_new = h.copy()
    hu_new = hu.copy()
    
    for i in range(1, N-1):
        h_new[i] = h[i] - dt/dx * (fluxes[0, i+1] - fluxes[0, i])
        hu_new[i] = hu[i] - dt/dx * (fluxes[1, i+1] - fluxes[1, i])
        
    # Apply reflective boundary conditions
    h_new[0] = h_new[1]
    hu_new[0] = -hu_new[1]
    h_new[N-1] = h_new[N-2]
    hu_new[N-1] = -hu_new[N-2]
    
    return h_new, hu_new

# putting it all together
def solve_SWE_1D(t, t_end, h, hu, N, dx, lim_func,CFL=0.9,g=9.81):
    # Store the intermediate states
    h_hist = []
    hu_hist = []
    dt_hist = []

    h_hist.append(h.copy())
    hu_hist.append(hu.copy())

    # Main loop
    while t < t_end:
        # Compute the time step
        u = hu/h
        c = np.sqrt(g*h)
        dt = CFL * dx/np.max(np.abs(u) + c)

        if t + dt > t_end:
            dt = t_end - t
        
        # Compute the fluxes at the interfaces
        fluxes = compute_interface_fluxes(h, hu, N, g, lim_func)

        # Update the variables
        h, hu = update_variables(h, hu, dt, fluxes, N, dx)

        # Store the results
        h_hist.append(h.copy())
        hu_hist.append(hu.copy())
        dt_hist.append(dt)

        # Update the time
        t += dt

    n = len(dt_hist)

    print(f"Number of time steps: {n}")

    return h, hu, h_hist

# Plot function

def subplot_solution(x, h, hu, h_hist, save=False, filename='VSC/tex/plots/solution.png'):
    plt.figure(figsize=(10,5))
    plt.subplot(1,4,1)
    plt.plot(x,h_hist[0])
    plt.title('Initial water height')
    plt.xlabel('Distance (x)')
    plt.ylabel('Water height (h)')

    plt.subplot(1,4,2)  
    plt.plot(x,h)
    plt.title('Water height (h)')
    plt.xlabel('Distance (x)')
    #plt.ylabel('Water height (h)')

    plt.subplot(1,4,3)
    plt.plot(x,hu/h)
    plt.title('Velocity (u)')
    plt.xlabel('Distance (x)')
    #plt.ylabel('Velocity (u)')

    plt.subplot(1,4,4)
    plt.plot(x,hu)
    plt.title('Momentum (hu)')
    plt.xlabel('Distance (x)')
    #plt.ylabel('Momentum (hu)')

    if save:
        plt.savefig(filename)

    plt.show()  

    return 