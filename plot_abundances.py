
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

stellar_parameters = Table.read("brewer_2016_table8.txt", format="cds")
abundances = Table.read("brewer_2016_table9.txt", format="cds")



similar_to_kronos = \
        (np.abs(stellar_parameters["Teff"] - 5803) < 100) \
    *   (np.abs(stellar_parameters["logg"] - 4.33) < 0.1) 

colors = {
    "HD 240429": "r",
    "HD 240430": "b"
}

# Plot abundances for each sequence.
Lodders_2003 = Table.read("Lodders.txt", format="ascii", delimiter="|")

elements = [el.strip("[]").split("/")[0] for el in abundances.dtype.names[2:]]

x = np.array([Lodders_2003["tc50"][Lodders_2003["el"] == el][0] for el in elements])
ordered = np.argsort(x)


fig, ax = plt.subplots()

for index in np.where(similar_to_kronos)[0]:

    match = abundances["Name"] == stellar_parameters["Name"][index]

    y = np.array([abundances[col][match][0] for col in abundances.dtype.names[2:]]) \
      - stellar_parameters["[M/H]"][index]
    

    kwds = dict(c=colors.get(stellar_parameters["Name"][index], None))

    if kwds["c"] is None:
        difference = np.median(y[x > 1200]) - np.median(y[x < 1200])
        if difference > 0.1:
            kwds.update(c="k", zorder=1, alpha=0.5)
        else:
            kwds.update(c="#CCCCCC", zorder=-10, alpha=0.5)
    else:
        kwds.update(zorder=10, alpha=0.9, label=stellar_parameters["Name"][index])

    ax.plot(x[ordered], y[ordered], **kwds)


plt.legend(loc="upper left")


for xi in x:
    ax.axvline(xi, c="#EEEEEE", linestyle=":", zorder=-100)

twin_y = ax.twiny()
twin_y.set_xlim(ax.get_xlim())
twin_y.set_xticks(x[ordered])
twin_y.set_xticklabels([elements[o] for o in ordered])

ax.set_xlabel(r"50% w/w condensation temperature (Lodders 2003) [K]")
ax.set_ylabel(r"[X/M]")

fig.tight_layout()

fig.savefig("kk_similar.png", dpi=300)


raise a