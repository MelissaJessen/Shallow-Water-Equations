import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from scipy.io import loadmat

import time
import h5py
import os


# Define the global gravity constant
g = 9.81

def FluxX(Q):
    # Initialize F to the same shape as Q
    F = np.zeros_like(Q)

    # Compute the flux in the x-direction    
    F[0, :, :] = Q[1, :, :]
    F[1, :, :] = Q[1, :, :]**2 / Q[0, :, :] + 0.5 * g * Q[0, :, :]**2
    F[2, :, :] = Q[1, :, :] * Q[2, :, :] / Q[0, :, :]

    return F


def FluxY(Q):
    # Initialize G to the same shape as Q
    G = np.zeros_like(Q)

    # Compute the flux in the y-direction
    G[0, :, :] = Q[2, :, :]
    G[1, :, :] = Q[1, :, :] * Q[2, :, :] / Q[0, :, :]
    G[2, :, :] = Q[2, :, :]**2 / Q[0, :, :] + 0.5 * g * Q[0, :, :]**2

    return G

def Lambdax(Q):
    # Initialize L to the same shape as Q
    L = np.zeros_like(Q)

    # Compute velocity and wave speed
    u = Q[1, :, :] / Q[0, :, :]  # u = hu / h
    c = np.sqrt(g * Q[0, :, :])   # c = sqrt(g * h)

    # Characteristic velocities in the x-direction
    L[0, :, :] = u - c
    L[1, :, :] = u
    L[2, :, :] = u + c

    return L

def Lambday(Q):
    # Initialize L to the same shape as Q
    L = np.zeros_like(Q)

    # Compute velocity in y-direction and wave speed
    v = Q[2, :, :] / Q[0, :, :]  # v = hv / h
    c = np.sqrt(g * Q[0, :, :])   # c = sqrt(g * h)

    # Characteristic velocities in the y-direction
    L[0, :, :] = v - c
    L[1, :, :] = v
    L[2, :, :] = v + c

    return L

def minmodarray(a, b):
    nVar, Nx, Ny = a.shape
    c = np.zeros((nVar, Nx, Ny))

    for i in range(nVar):
        # First condition: a * b > 0 and |a| < |b|
        logic1 = (a[i, :, :] * b[i, :, :] > 0) & (np.abs(a[i, :, :]) < np.abs(b[i, :, :]))
        c[i, :, :] = a[i, :, :] * logic1
        
        # Second condition: a * b > 0 and |a| > |b|
        logic2 = (a[i, :, :] * b[i, :, :] > 0) & (np.abs(a[i, :, :]) > np.abs(b[i, :, :]))
        c[i, :, :] += b[i, :, :] * logic2

    return c

def Rusanov(QL, QR, FL, FR, sL, sR):
    smax = np.maximum(np.abs(sL).max(), np.abs(sR).max())
    F = 0.5 * (FR + FL) - 0.5 * smax * (QR - QL)
    return F

def Rusanov(QL, QR, FL, FR, sL, sR):
    smax = np.maximum(np.abs(sL).max(), np.abs(sR).max())
    QL = QL.reshape(3,)
    QR = QR.reshape(3,)
    FL = FL.reshape(3,)
    FR = FR.reshape(3,)
    sL = sL.reshape(3,)
    sR = sR.reshape(3,)
    F = 0.5 * (FR + FL) - 0.5 * smax * (QR - QL)
    return F



def solve_2D_SWE_FVM(Q, Qnew, t, tend, CFL, dx, dy, Nx, Ny):
    Q_all = []
    Q_all.append(Q)
    t_all = []
    t_all.append(t)

     # Time loop
    for n in range(100000):
        sx = Lambdax(Q)  # characteristic velocities in x-direction
        sy = Lambday(Q)  # characteristic velocities in y-direction
        ax = np.max(np.abs(sx))
        ay = np.max(np.abs(sy))
        
        dt = CFL / (ax / dx + ay / dy)
        if t + dt > tend:
            dt = tend - t
        if t >= tend:
            break

        # MUSCL part
        slopeX = np.zeros_like(Q)
        slopeY = np.zeros_like(Q)
        
        slopeX[:, 1:Nx-1, :] = minmodarray(Q[:, 1:Nx-1, :] - Q[:, 0:Nx-2, :], Q[:, 2:Nx, :] - Q[:, 1:Nx-1, :])
        slopeY[:, :, 1:Ny-1] = minmodarray(Q[:, :, 1:Ny-1] - Q[:, :, 0:Ny-2], Q[:, :, 2:Ny] - Q[:, :, 1:Ny-1])

        Qxm = Q - 0.5 * slopeX
        Qxp = Q + 0.5 * slopeX
        Qym = Q - 0.5 * slopeY
        Qyp = Q + 0.5 * slopeY

        Q_t = -(FluxX(Qxp) - FluxX(Qxm)) / dx - (FluxY(Qyp) - FluxY(Qym)) / dy
        Qxm += 0.5 * dt * Q_t
        Qxp += 0.5 * dt * Q_t
        Qym += 0.5 * dt * Q_t
        Qyp += 0.5 * dt * Q_t

        # Recompute physical fluxes and characteristic velocities for Rusanov flux
        fxm, fxp = FluxX(Qxm), FluxX(Qxp)
        gym, gyp = FluxY(Qym), FluxY(Qyp)
        
        sxm, sxp = Lambdax(Qxm), Lambdax(Qxp)
        sym, syp = Lambday(Qym), Lambday(Qyp)

        dtdx = dt / dx
        dtdy = dt / dy

        # Space loop
        for i in range(Nx):
            for j in range(Ny):
                # Numerical fluxes in X
                if i == 0:
                    Qghost = Q[:, i, j].reshape(3, 1, 1).copy()  # Reshape to (3, 1, 1)
                    Qghost[1] = -Qghost[1]
                    Fm = Rusanov(Qghost, Qxm[:, i, j], FluxX(Qghost), fxm[:, i, j], Lambdax(Qghost), sxm[:, i, j])
                    Fp = Rusanov(Qxp[:, i, j], Qxm[:, i+1, j], fxp[:, i, j], fxm[:, i+1, j], sxp[:, i, j], sxm[:, i+1, j])
                elif i == Nx - 1:
                    Qghost = Q[:, i, j].reshape(3, 1, 1).copy()  # Reshape to (3, 1, 1)
                    Qghost[1] = -Qghost[1]
                    Fm = Rusanov(Qxp[:, i-1, j], Qxm[:, i, j], fxp[:, i-1, j], fxm[:, i, j], sxp[:, i-1, j], sxm[:, i, j])
                    Fp = Rusanov(Qxp[:, i, j], Qghost, fxp[:, i, j], FluxX(Qghost), sxp[:, i, j], Lambdax(Qghost))
                else:
                    Fm = Rusanov(Qxp[:, i-1, j], Qxm[:, i, j], fxp[:, i-1, j], fxm[:, i, j], sxp[:, i-1, j], sxm[:, i, j])
                    Fp = Rusanov(Qxp[:, i, j], Qxm[:, i+1, j], fxp[:, i, j], fxm[:, i+1, j], sxp[:, i, j], sxm[:, i+1, j])

                # Numerical fluxes in Y
                if j == 0:
                    Qghost = Q[:, i, j].reshape(3, 1, 1).copy()  # Reshape to (3, 1, 1)
                    Qghost[2] = -Qghost[2]
                    Gm = Rusanov(Qghost, Qym[:, i, j], FluxY(Qghost), gym[:, i, j], Lambday(Qghost), sym[:, i, j])
                    Gp = Rusanov(Qyp[:, i, j], Qym[:, i, j+1], gyp[:, i, j], gym[:, i, j+1], syp[:, i, j], sym[:, i, j+1])
                elif j == Ny - 1:
                    Qghost = Q[:, i, j].reshape(3, 1, 1).copy()  # Reshape to (3, 1, 1)
                    Qghost[2] = -Qghost[2]
                    Gm = Rusanov(Qyp[:, i, j-1], Qym[:, i, j], gyp[:, i, j-1], gym[:, i, j], syp[:, i, j-1], sym[:, i, j])
                    Gp = Rusanov(Qyp[:, i, j], Qghost, gyp[:, i, j], FluxY(Qghost), syp[:, i, j], Lambday(Qghost))
                else:
                    Gm = Rusanov(Qyp[:, i, j-1], Qym[:, i, j], gyp[:, i, j-1], gym[:, i, j], syp[:, i, j-1], sym[:, i, j])
                    Gp = Rusanov(Qyp[:, i, j], Qym[:, i, j+1], gyp[:, i, j], gym[:, i, j+1], syp[:, i, j], sym[:, i, j+1])

                # Finite volume update
                Qnew[:, i, j] = Q[:, i, j] - dtdx * (Fp - Fm) - dtdy * (Gp - Gm)
                

        # Update time and solution
        t += dt
        # Print t with 3 decimals
        print(f't = {t:.3f}')
        t_all.append(t)
        Q = np.copy(Qnew)
        Q_all.append(Q)

    n = len(Q_all)
    print(f'There are time steps: {n}')

    return Q_all, t_all




def solve_2D_SWE_FVM_constant_dt(Q, Qnew, t, tend, CFL, dx, dy, Nx, Ny):
    Q_all = []
    Q_all.append(Q)
    t_all = []
    t_all.append(t)

     # Time loop
    for n in range(100000):
        sx = Lambdax(Q)  # characteristic velocities in x-direction
        sy = Lambday(Q)  # characteristic velocities in y-direction
        ax = np.max(np.abs(sx))
        ay = np.max(np.abs(sy))
        
        #dt = CFL / (ax / dx + ay / dy)
        dt = 0.025
        if t + dt > tend:
            dt = tend - t
        if t >= tend:
            break

        # MUSCL part
        slopeX = np.zeros_like(Q)
        slopeY = np.zeros_like(Q)
        
        slopeX[:, 1:Nx-1, :] = minmodarray(Q[:, 1:Nx-1, :] - Q[:, 0:Nx-2, :], Q[:, 2:Nx, :] - Q[:, 1:Nx-1, :])
        slopeY[:, :, 1:Ny-1] = minmodarray(Q[:, :, 1:Ny-1] - Q[:, :, 0:Ny-2], Q[:, :, 2:Ny] - Q[:, :, 1:Ny-1])

        Qxm = Q - 0.5 * slopeX
        Qxp = Q + 0.5 * slopeX
        Qym = Q - 0.5 * slopeY
        Qyp = Q + 0.5 * slopeY

        Q_t = -(FluxX(Qxp) - FluxX(Qxm)) / dx - (FluxY(Qyp) - FluxY(Qym)) / dy
        Qxm += 0.5 * dt * Q_t
        Qxp += 0.5 * dt * Q_t
        Qym += 0.5 * dt * Q_t
        Qyp += 0.5 * dt * Q_t

        # Recompute physical fluxes and characteristic velocities for Rusanov flux
        fxm, fxp = FluxX(Qxm), FluxX(Qxp)
        gym, gyp = FluxY(Qym), FluxY(Qyp)
        
        sxm, sxp = Lambdax(Qxm), Lambdax(Qxp)
        sym, syp = Lambday(Qym), Lambday(Qyp)

        dtdx = dt / dx
        dtdy = dt / dy

        # Space loop
        for i in range(Nx):
            for j in range(Ny):
                # Numerical fluxes in X
                if i == 0:
                    Qghost = Q[:, i, j].reshape(3, 1, 1).copy()  # Reshape to (3, 1, 1)
                    Qghost[1] = -Qghost[1]
                    Fm = Rusanov(Qghost, Qxm[:, i, j], FluxX(Qghost), fxm[:, i, j], Lambdax(Qghost), sxm[:, i, j])
                    Fp = Rusanov(Qxp[:, i, j], Qxm[:, i+1, j], fxp[:, i, j], fxm[:, i+1, j], sxp[:, i, j], sxm[:, i+1, j])
                elif i == Nx - 1:
                    Qghost = Q[:, i, j].reshape(3, 1, 1).copy()  # Reshape to (3, 1, 1)
                    Qghost[1] = -Qghost[1]
                    Fm = Rusanov(Qxp[:, i-1, j], Qxm[:, i, j], fxp[:, i-1, j], fxm[:, i, j], sxp[:, i-1, j], sxm[:, i, j])
                    Fp = Rusanov(Qxp[:, i, j], Qghost, fxp[:, i, j], FluxX(Qghost), sxp[:, i, j], Lambdax(Qghost))
                else:
                    Fm = Rusanov(Qxp[:, i-1, j], Qxm[:, i, j], fxp[:, i-1, j], fxm[:, i, j], sxp[:, i-1, j], sxm[:, i, j])
                    Fp = Rusanov(Qxp[:, i, j], Qxm[:, i+1, j], fxp[:, i, j], fxm[:, i+1, j], sxp[:, i, j], sxm[:, i+1, j])

                # Numerical fluxes in Y
                if j == 0:
                    Qghost = Q[:, i, j].reshape(3, 1, 1).copy()  # Reshape to (3, 1, 1)
                    Qghost[2] = -Qghost[2]
                    Gm = Rusanov(Qghost, Qym[:, i, j], FluxY(Qghost), gym[:, i, j], Lambday(Qghost), sym[:, i, j])
                    Gp = Rusanov(Qyp[:, i, j], Qym[:, i, j+1], gyp[:, i, j], gym[:, i, j+1], syp[:, i, j], sym[:, i, j+1])
                elif j == Ny - 1:
                    Qghost = Q[:, i, j].reshape(3, 1, 1).copy()  # Reshape to (3, 1, 1)
                    Qghost[2] = -Qghost[2]
                    Gm = Rusanov(Qyp[:, i, j-1], Qym[:, i, j], gyp[:, i, j-1], gym[:, i, j], syp[:, i, j-1], sym[:, i, j])
                    Gp = Rusanov(Qyp[:, i, j], Qghost, gyp[:, i, j], FluxY(Qghost), syp[:, i, j], Lambday(Qghost))
                else:
                    Gm = Rusanov(Qyp[:, i, j-1], Qym[:, i, j], gyp[:, i, j-1], gym[:, i, j], syp[:, i, j-1], sym[:, i, j])
                    Gp = Rusanov(Qyp[:, i, j], Qym[:, i, j+1], gyp[:, i, j], gym[:, i, j+1], syp[:, i, j], sym[:, i, j+1])

                # Finite volume update
                Qnew[:, i, j] = Q[:, i, j] - dtdx * (Fp - Fm) - dtdy * (Gp - Gm)
                

        # Update time and solution
        t += dt
        # Print t with 3 decimals
        print(f't = {t:.3f}')
        t_all.append(t)
        Q = np.copy(Qnew)
        Q_all.append(Q)

    n = len(Q_all)
    print(f'There are time steps: {n}')

    return Q_all, t_all




