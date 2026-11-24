"""Compile et lance les quatre programmes, vérifie les ordres et trace les figures."""
import math
import pathlib
import subprocess
import tempfile

import matplotlib.pyplot as plt
import numpy as np

ICI = pathlib.Path(__file__).resolve().parent
(ICI / "figures").mkdir(exist_ok=True)


def lancer(dossier, parametres, sortie):
    """Lance l'exécutable dans un dossier temporaire et renvoie le contenu du fichier de sortie."""
    subprocess.run(["make", "-s", "-C", str(ICI / dossier)], check=True, stderr=subprocess.DEVNULL)
    with tempfile.TemporaryDirectory() as tmp:
        (pathlib.Path(tmp) / "parametres.dat").write_text(parametres + "\n")
        subprocess.run([str(ICI / dossier / "exe")], cwd=tmp, check=True, stdout=subprocess.DEVNULL)
        return (pathlib.Path(tmp) / sortie).read_text()


def figure(nom, xlabel, ylabel):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(which="both", alpha=0.3)
    plt.savefig(ICI / "figures" / nom, dpi=150, bbox_inches="tight")
    plt.close()


def convergence(dossier, cas, nom):
    """cas : (légende, paramètres, ordre minimal). Ordre mesuré en norme infinie, deux grilles les plus fines."""
    for legende, parametres, mini in cas:
        dx, e = np.array([l.split()[:2] for l in lancer(dossier, parametres, "err.dat").splitlines()], float).T
        ordre = math.log(e[-2] / e[-1]) / math.log(dx[-2] / dx[-1])
        print(f"{dossier:10s} {legende:28s} ordre {ordre:.2f} (min {mini}), erreur {e[-1]:.2e}", flush=True)
        assert ordre > mini - 0.15
        plt.loglog(dx, e, "o-", label=f"{legende} (ordre {ordre:.2f})")
    figure(nom, "dx", "erreur max")


convergence("derivees", [
    ("centré", "1 100 1 1 6", 2),
    ("décentré avant", "1 100 2 1 6", 1),
    ("décentré à 4 points", "1 100 3 1 6", 3),
], "derivees.png")

convergence("poisson_1d", [
    ("3 points, régulier", "1 40 1 1 5 0.1 0.1", 2),
    ("5 points, régulier", "1 40 2 1 5 0.1 0.1", 4),
    ("3 points, raide", "1 40 1 5 5 0 1", 2),
    ("5 points, raide", "1 40 2 5 5 0 1", 3),
], "poisson_1d.png")

convergence("chaleur", [
    ("Euler explicite", "1.0 8 1 1 4 0 0 0.1 0.9", 2),
    ("Euler implicite", "1.0 8 2 1 4 0 0 0.1 0.9", 2),
], "chaleur.png")


# Transport-diffusion : créneau entrant (Ogata-Banks) + gaussienne initiale
A, DIFF = 0.5, 1.0e-3


def ogata_banks(x, t):
    r = 2 * math.sqrt(DIFF * t)
    s = np.array([math.erfc((v - A * t) / r) for v in x])
    for i, v in enumerate(x):  # terme de bord : exp() déborde loin du bord, où erfc() est nul
        if A * v / DIFF < 600:
            s[i] += math.exp(A * v / DIFF) * math.erfc((v + A * t) / r)
    return 0.5 * s


def exacte(x, t):
    var = 0.005 + 2 * DIFF * t
    g = math.sqrt(0.005 / var) * np.exp(-((x - 0.7 - A * t) ** 2) / (2 * var))
    return 30 * (ogata_banks(x, t - 2) - ogata_banks(x, t - 3)) + g


def instantane(schema, nx, t_cible):
    blocs = []
    for bloc in lancer("transport", f"10.0 {nx} {schema} 1 1 0.5 0.5", "solution.dat").split("\n \n \n"):
        lignes = [l for l in bloc.splitlines() if l.strip()]
        if lignes:
            blocs.append(np.array([l.split() for l in lignes], float))
    b = min(blocs, key=lambda b: abs(b[0, 2] - t_cible))
    return b[:, 0], b[:, 1], b[0, 2]


for schema, legende in ((3, "Lax-Wendroff"), (1, "amont")):
    errs = []
    for nx in (200, 400, 800, 1600):
        x, u, t = instantane(schema, nx, 10.0)
        errs.append(np.sum(np.abs(u - exacte(x, t))) * (x[1] - x[0]))
        print(f"transport  {legende:14s} nx={nx:5d} t={t:.2f} erreur L1 {errs[-1]:.3e}", flush=True)
    assert all(b < a for a, b in zip(errs, errs[1:]))

for schema, legende in ((1, "amont"), (3, "Lax-Wendroff")):
    x, u, t = instantane(schema, 800, 10.0)
    plt.plot(x, u, label=f"{legende}, t = {t:.2f}")
plt.plot(x, exacte(x, t), "k--", label=f"exacte, t = {t:.2f}")
plt.xlim(2.5, 6.5)
figure("transport.png", "x", "concentration (nx = 800)")
