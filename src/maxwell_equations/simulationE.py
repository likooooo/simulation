#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from math import pow
from scipy.constants import epsilon_0
from visualizer import display_image

def E(q, r0, r1, k = 1/ (4 * np.pi * epsilon_0)):
    r_diff = r1 - r0
    distance = np.linalg.norm(r_diff, axis=-1, keepdims=True)
    E_magnitude =  k * q / distance**2
    E_direction = r_diff / distance
    return E_magnitude * E_direction
def U(q, r0, r1, k = 1/ (4 * np.pi * epsilon_0)):
    r_diff = r1 - r0
    distance = np.linalg.norm(r_diff, axis=-1, keepdims=True)
    s = distance.shape
    if s[-1] == 1 : distance.shape = s[0: len(s) - 1]
    u = k * q / distance
    return u
def symmetry_charges(nq, center, amp_q = 1, r=1):
    if 1 == nq : r = 0
    center = np.array(center)
    n = np.array(range(nq))
    q = [pow(-1 * amp_q, i) for i in n]
    theta = n * 2 * np.pi / nq
    charges = list(zip(q, r * np.cos(theta) , r * np.sin(theta)))
    return charges
def electric_field(x, y, charges):
    Ex, Ey = np.zeros((len(y), len(x))), np.zeros((len(y), len(x)))
    X, Y = np.meshgrid(x, y)
    r1 = np.stack((X, Y), axis=-1)
    for charge in charges:
        E_field = E(charge[0], charge[1 : len(charge)],r1)
        ex, ey = E_field[..., 0], E_field[..., 1]
        Ex += ex
        Ey += ey
    display_image(np.log(np.hypot(Ex, Ey)))
    return Ex, Ey
def electric_potential(x, y, charges):
    U_array = np.zeros((len(y), len(x)))
    X, Y = np.meshgrid(x, y)
    r1 = np.stack((X, Y), axis=-1)
    for charge in charges:
        delta_U  = U(charge[0], charge[1 : len(charge)],r1)
        U_array += delta_U
    print(U_array.shape)
    display_image(U_array)
    Ey, Ex = np.gradient(-U_array)
    return Ex, Ey
def plot_stream(charges, Ex, Ey, x, y, xrange, yrange):
    fig = plt.figure()
    ax = fig.add_subplot()

    # Plot the streamlines with an appropriate colormap and arrow style
    color = 2 * np.log(np.hypot(Ex, Ey))
    ax.streamplot(x, y, Ex, Ey, color=color, linewidth=1, cmap=plt.cm.inferno,
                density=2, arrowstyle='->', arrowsize=1.5)

    # Add filled circles for the charges themselves
    charge_colors = {True: '#aa0000', False: '#0000aa'}
    for charge in charges:
        q = charge[0]
        pos = charge[1 :len(charge)]
        ax.add_artist(Circle(pos, 0.05, color=charge_colors[q > 0]))

    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_xlim(xrange[0], xrange[1])
    ax.set_ylim(yrange[0], yrange[1])
    ax.set_aspect('equal')
    plt.show()

# Grid of x, y points
nx, ny = 64, 64
lx, ux = -2, 2
ly, uy = -2, 2
nq = 2 # the number of charges
def case0():
    charges = symmetry_charges(nq, [(lx + ly)/2, (ux + uy)/2, 1])
    x = np.linspace(lx, ux, nx)
    y = np.linspace(ly, uy, ny)
    Ex, Ey = electric_field(x, y, charges)
    plot_stream(charges, Ex, Ey, x, y, [lx, ux], [ly, uy])
    return Ex, Ey

def case1():
    charges = symmetry_charges(nq, [(lx + ly)/2, (ux + uy)/2, 1])
    x = np.linspace(lx, ux, nx)
    y = np.linspace(ly, uy, ny)
    Ex, Ey = electric_potential(x, y, charges)
    plot_stream(charges, Ex, Ey, x, y, [lx, ux], [ly, uy])
    return Ex, Ey

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="electric filed simulation")
    parser.add_argument("--start", nargs="+", type=float, default=[lx, ly])
    parser.add_argument("--end", nargs="+", type=float, default=[ux, uy])
    parser.add_argument("--grid", nargs="+", type=int, default=[nx, ny])
    parser.add_argument("--n", type=int, default=nq)
    args = parser.parse_args() 
    nx, ny = args.grid[0:2] 
    lx, ly = args.start[0:2] 
    ux, uy = args.end[0:2] 
    nq = args.n
    case1()