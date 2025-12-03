from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

INPUT_TXT : Path = Path("bench.txt")
VAL_NAMES : list[str] = ["temps total", "temps init", "temps lecture", "temps matrice", "temps floyd", "temps pam"]
lines : list[str] = [ ]

mean_times : list[tuple[str, list[float]]] = [ ]

current_values : list[list[float]] = []
current_name : str

# open the bench file
with INPUT_TXT.open() as file:
    lines = file.readlines()

# read the values
for line in lines :
    if line.startswith(">") :
        # save the mean values to the list
        if current_values != [] :
            values : list[float] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            # since every bench end w/ _____, there is one more [] at the end of current_values
            for v in current_values[:-1] :
                for i in range(6):
                    values[i] += v[i]
            for i in range(6):
                values[i] = round(values[i] / (len(current_values)-1), 3)

            mean_times.append((current_name, values))
        # next set of benchs
        current_name = line.removeprefix("> ").removesuffix("\n")
        current_values = []
    elif line.startswith("_") :
        current_values.append([])
    elif line.startswith("temps"):
        current_values[-1].append(float(line.split(" : ")[1].removesuffix("ms\n")))


# what do we want to draw ?
TARGET_CHRONOS : list[str] = ["temps total", "temps matrice", "temps floyd", "temps pam"] 

x : list[int] = [1, 4, 16]
ys2k : list[list[float]] = [[0.0 for _ in range(3)] for __ in range(6)] 
ys500 : list[list[float]] = [[0.0 for _ in range(3)] for __ in range(6)] 

for bench in mean_times:
        i = x.index(int(bench[0].split(" ")[2]))
        for j in range(len(bench[1])):
            if bench[0].startswith("2k"):
                ys2k[j][i] = bench[1][j]
            else :
                ys500[j][i] = bench[1][j]

# times for 2k nodes
fig, axis = plt.subplots(2, 2, sharey=True)
fig.set_size_inches(18, 10)

for t in range(len(TARGET_CHRONOS)):
    i = VAL_NAMES.index(TARGET_CHRONOS[t])
    ax = axis[t//2][t%2]
    ax.plot(x, ys2k[i], "b-o")
    for j in range(len(x)):
        xval = x[j]
        yval = ys2k[i][j]
        if j == 0:
            if TARGET_CHRONOS[t] == "temps total" :
                ax.text(xval+.4, yval-700, f"{yval}ms", horizontalalignment="left")
            else :
                ax.text(xval-.5, yval+2000, f"{yval}ms", horizontalalignment="left")
        elif j == 1 :
            ax.text(xval, yval+2000, f"{yval}ms", horizontalalignment="left")
        else :
            ax.text(xval, yval+2000, f"{yval}ms", horizontalalignment="right")
    ax.set_title(f"[{TARGET_CHRONOS[t].split(" ")[1].upper()}]")
    ax.set(xlabel="# de coeurs", ylabel="temps (ms)", xticks=x)

plt.savefig(f"../asset/2k_nodes.svg")

plt.clf()

# times for 500 nodes

fig, axis = plt.subplots(2, 2, sharey=True)
fig.set_size_inches(18, 10)

for t in range(len(TARGET_CHRONOS)):
    i = VAL_NAMES.index(TARGET_CHRONOS[t])
    ax = axis[t//2][t%2]
    ax.plot(x, ys500[i], "b-o")
    for j in range(len(x)):
        xval = x[j]
        yval = ys500[i][j]
        if j == 0:
            if TARGET_CHRONOS[t] == "temps total" or TARGET_CHRONOS[t] == "temps pam" :
                ax.text(xval+.5, yval-20, f"{yval}ms", horizontalalignment="left")
            else :
                ax.text(xval-.5, yval+40, f"{yval}ms", horizontalalignment="left")
        elif j == 1 :
            ax.text(xval, yval+40, f"{yval}ms", horizontalalignment="left")
        else :
            ax.text(xval, yval+40, f"{yval}ms", horizontalalignment="right")
    ax.set_title(f"[{TARGET_CHRONOS[t].split(" ")[1].upper()}]")
    ax.set(xlabel="# de coeurs", ylabel="temps (ms)", xticks=x)

plt.savefig(f"../asset/500_nodes.svg")


# cumulative times (500)

plt.clf()

fig, ax = plt.subplots(1)
fig.set_size_inches(18, 10)

y = [0.0, 0.0, 0.0]
ym1 : list[float]

colors : str = ["hotpink", "gray",  "turquoise", "gold", "limegreen", "royalblue", "rebeccapurple"]

for i in range (1, 6):
    ym1 = y.copy()
    y = np.add(y, ys500[i])
    label = VAL_NAMES[i]
    if (VAL_NAMES[i] == "temps lecture") :
        label += "\n    (negligeable)"
    ax.fill_between(x, ym1, y, label=label, color=colors[i-1])
    if (i == 5):
        ax.plot(x, y, "k-o", label="temps total")
    # else :
        # ax.plot(x, y, "k--")
    for j in range(3):
        ax.text(x[j], y[j]-((((y[j]-ym1[j])/2)+5) if i != 2 else 20), f"{'↕' if j != 2 else ''} {round((ys500[i][j]/ys500[0][j])*100, 2)}% {'↕' if j == 2 else ''}",horizontalalignment=("left" if j != 2 else "right"))

ax.vlines([1, 4, 16], 0, ys500[0], "k", "dotted")

ax.legend()
ax.set(xticks=x, xlabel="# de coeurs",  ylabel="temps cumulé (ms)")

fig.savefig("../asset/temps_cumulés_500.svg")