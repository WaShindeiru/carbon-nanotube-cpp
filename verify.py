import numpy as np
import matplotlib.pyplot as plt

def show_me():
    rows = []
    with open('log') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 10:
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    pass  # skips header lines like "N =   500"

    data = np.array(rows)

    data = data[1:, :]

    # data = np.loadtxt('log.old')

    it    = data[:, 0]           # col 1: iteration
    t     = data[:, 1]           # col 2: physical time [s]
    Tinst = data[:, 2]           # col 3: instantaneous temperature [K]
    Tavg  = data[:, 3]           # col 4: average temperature [K]
    Epot  = data[:, 4]           # col 5: potential energy
    Ekin  = data[:, 5]           # col 6: kinetic energy
    Etot  = data[:, 6]           # col 7: total energy

    title  = 'N=500, dt=0.05,  fgrav=0.005'
    xlabel = 'time [s]'

    plt.rcParams.update({'font.size': 20, 'font.family': 'sans-serif'})

    # --- temperature_average_500it.png ---
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(t, Tavg, marker='o', markersize=3, label=r'$T_{\mathrm{average\ 500it}}$')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    # ax.set_xlim(300, None)
    # ax.set_ylim(0, 2000)
    ax.legend()
    fig.savefig('temperature_average_500it.png', dpi=100)
    plt.close(fig)

    # --- energy_temporary.png ---
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(t, Epot, marker='o', markersize=3, label='Epot')
    ax.plot(t, Ekin, marker='o', markersize=3, label='Ekin')
    ax.plot(t, Etot, marker='o', markersize=3, label='Etot')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    # ax.set_xlim(300, None)
    ax.legend()
    fig.savefig('energy_temporary.png', dpi=100)
    plt.close(fig)

    # --- etot_temporary.png ---
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(t, Etot, marker='o', markersize=3, label='Etot')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    # ax.set_xlim(300, None)
    # ax.set_ylim(-1800, -1700)
    ax.legend()
    fig.savefig('etot_temporary.png', dpi=100)
    plt.close(fig)

    ref_idx = 0  # placeholder: index of the reference iteration

    Etot = Etot - Etot[ref_idx]

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(t[ref_idx:], Etot[ref_idx:], marker='o', markersize=3, label=r'$E_{tot} - E_{tot}^{ref}$')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Energy difference')
    ax.legend()
    fig.savefig('etot_diff.png', dpi=100)
    plt.close(fig)

    

if __name__ == "__main__":
    show_me()