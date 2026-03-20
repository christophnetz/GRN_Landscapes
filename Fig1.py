import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.integrate import solve_ivp
from scipy.optimize import root

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
})

# -----------------------
# Individual fitness peaks
# -----------------------
mu1 = np.array([1,2])
mu2 = np.array([3,2])

# -----------------------
# Mean fitness functions
# -----------------------
def mean_fitness(z, B):
    """Mean population fitness (constant factor removed)."""
    g1 = np.exp(-(z-mu1).T @ B @ (z-mu1))
    g2 = np.exp(-(z-mu2).T @ B @ (z-mu2))
    return g1 + g2


def grad_log_mean_fitness(z, B):
    """Gradient of log mean population fitness."""
    g1 = np.exp(-(z-mu1).T @ B @ (z-mu1))
    g2 = np.exp(-(z-mu2).T @ B @ (z-mu2))

    W = g1 + g2

    grad1 = -2 * B @ (z-mu1) * g1
    grad2 = -2 * B @ (z-mu2) * g2

    gradW = grad1 + grad2

    return gradW / W


# -----------------------
# Evolutionary dynamics
# -----------------------

def dynamics(t, z, G, B):
    beta = grad_log_mean_fitness(z, B)
    return G @ beta


from scipy.optimize import minimize

def find_mean_fitness_maxima(B, initial_guesses):
    """
    Find maxima of mean fitness landscape given B (=inv(I + 2P))
    
    Parameters:
    - B: 2x2 matrix for mean fitness
    - initial_guesses: list of 2D points to start optimization
    
    Returns:
    - list of maxima as 2D numpy arrays
    """
    maxima = []
    for guess in initial_guesses:
        # minimize negative mean fitness
        res = minimize(lambda z: -mean_fitness(z, B), guess, method='BFGS')
        if res.success:
            maxima.append(res.x)
    return maxima


def draw_g_matrix(ax, G, center=(0.85, 1.05), scale=0.25, color='black', lw=2):
    """
    Draw a G-matrix ellipse with a cross through the center along eigenvectors.
    """
    # Eigen decomposition
    eigvals, eigvecs = np.linalg.eigh(G)
    # Sort descending
    order = eigvals.argsort()[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    # Ellipse dimensions in axes coordinates
    width = scale * np.sqrt(eigvals[0])
    height = scale * np.sqrt(eigvals[1])

    # Draw ellipse
    angle = np.degrees(np.arctan2(eigvecs[1,0], eigvecs[0,0]))
    ellipse = Ellipse(
        center,
        width=width,
        height=height,
        angle=angle,
        facecolor='none',
        edgecolor=color,
        lw=lw,
        transform=ax.transAxes,
        clip_on=False
    )
    ax.add_patch(ellipse)

    # Draw cross along eigenvectors
    cx, cy = center

    # Loop over eigenvectors
    for i in range(2):
        # Unit vector in data space
        vec = eigvecs[:, i] / np.linalg.norm(eigvecs[:, i])
        # Corresponding axis length
        if i == 0:
            length = width / 2
        else:
            length = height / 2
        # Rotate the unit vector to match ellipse orientation in axes coordinates
        # For axes coordinates, rotation is already included in vec for ellipse
        # Draw both directions
        for sign in [-1, 1]:
            dx = sign * length * vec[0]
            dy = sign * length * vec[1]
            ax.plot(
                [cx - dx, cx + dx],
                [cy - dy, cy + dy],
                color=color,
                lw=lw,
                transform=ax.transAxes,
                clip_on=False
            )
# -----------------------
# Grid
# -----------------------
x = np.linspace(0,4,300)
y = np.linspace(0,4,300)

X, Y = np.meshgrid(x,y)

# -----------------------
# G matrices
# -----------------------
G_matrices = [
    np.array([[0.1,0],[0,0.1]]),
    np.array([[0.1,0.05],[0.05,0.1]]),
    np.array([[0.1,-0.05],[-0.05,0.1]])
]

# -----------------------
# Figure setup
# -----------------------

n = len(G_matrices)

fig = plt.figure(figsize=(6*n,6))
gs = fig.add_gridspec(1, n+1, width_ratios=[1]*n + [0.05], wspace=0.3)

axes = [fig.add_subplot(gs[0,i]) for i in range(n)]
cax = fig.add_subplot(gs[0,n])


import string
subplot_labels = list(string.ascii_uppercase[:len(G_matrices)])


# -----------------------
# Loop over G matrices
# -----------------------

for index, (ax, G, label) in enumerate(zip(axes, G_matrices, subplot_labels)):
    P = G
    I = np.eye(2)

    # Precompute smoothing matrix
    B = np.linalg.inv(I + 2*P)

    # -----------------------
    # Compute mean fitness landscape
    # -----------------------

    Z = np.zeros_like(X)

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            z = np.array([X[i,j], Y[i,j]])
            Z[i,j] = mean_fitness(z, B)

    # -----------------------
    # Plot fitness landscape
    # -----------------------

    cf = ax.contourf(X, Y, Z, levels=60, cmap='viridis', alpha=0.8)
    ax.contour(X, Y, Z, levels=20, colors='black', linewidths=1)

    # -----------------------
    # Vector field
    # -----------------------

    U = np.zeros_like(X)
    V = np.zeros_like(Y)

    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            z = np.array([X[i,j],Y[i,j]])
            dz = dynamics(0,z,G,B)
            U[i,j] = dz[0]
            V[i,j] = dz[1]

    skip = 15

    ax.quiver(
        X[::skip,::skip],
        Y[::skip,::skip],
        U[::skip,::skip],
        V[::skip,::skip],
        color='black',
        scale=None,
        scale_units='xy',
        angles='xy'
    )

    # -----------------------
    # Find saddle point
    # -----------------------

    sol = root(lambda p: grad_log_mean_fitness(p, B), [2,2])
    saddle = sol.x

    # -----------------------
    # Jacobian at saddle
    # -----------------------

    h = 1e-5
    J = np.zeros((2,2))

    for i in range(2):

        e = np.zeros(2)
        e[i] = h

        f1 = dynamics(0, saddle + e, G, B)
        f2 = dynamics(0, saddle - e, G, B)

        J[:,i] = (f1 - f2)/(2*h)

    
    eigvals, eigvecs = np.linalg.eig(J)

    # Only draw if it's a proper saddle
    if np.any(np.real(eigvals) < 0) and np.any(np.real(eigvals) > 0):
        # valid saddle
        stable_index = np.where(np.real(eigvals) < 0)[0][0]
        stable_vec = np.real(eigvecs[:, stable_index])
        stable_vec /= np.linalg.norm(stable_vec)

        eps = 1e-4
        for sign in [-1,1]:
            z0 = saddle + sign*stable_vec*eps
            sol_forward = solve_ivp(dynamics, [0,-400], z0, args=(G,B), max_step=0.05)
            ax.plot(sol_forward.y[0], sol_forward.y[1], 'r', lw=2)
        

    # -----------------------
    # Separatrix
    # -----------------------
    # -----------------------
    # Fixed points
    # -----------------------
    # --- Fixed points ---
    # Compute mean-fitness maxima
    initial_guesses = [mu1, mu2]  # start near old peaks
    mean_peaks = find_mean_fitness_maxima(B, initial_guesses)

    # Plot them
    for mp in mean_peaks:
        ax.plot(mp[0], mp[1], 'wo', markersize=8)  # dots for mean-fitness peaks
    

    ax.plot(saddle[0],saddle[1],'ro',markersize=6)

    ax.set_xlim(0,4)
    ax.set_ylim(0,4)
    ax.set_aspect('equal')

    ax.set_title(
        r"$G = \left[ \begin{array}{cc} "
        f"{G[0,0]:.2f} & {G[0,1]:.2f} \\\\ "
        f"{G[1,0]:.2f} & {G[1,1]:.2f} "
        r"\end{array} \right]$",
        fontsize=14
    )
    
    draw_g_matrix(ax, G, center=(0.95, 1.12), scale=0.6, color='black', lw=2)     
      
    ax.text(
        -0.1,1.05,label,
        transform=ax.transAxes,
        fontsize=18,
        fontweight='bold',
        va='top',
        ha='right'
    )
    ax.set_xlabel("phenotype x")

    if index == 0:
        ax.set_ylabel("phenotype y")
    else:
        ax.set_ylabel("")
        ax.set_yticklabels([])
        

fig.colorbar(cf, cax=cax, label=r'mean population fitness $\bar{W}(\bar{z})$')
plt.savefig("mean_fitness_landscapes.png", dpi=300, bbox_inches='tight')
plt.show()
