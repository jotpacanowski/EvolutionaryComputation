# %%
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# import typer

# FIGSIZE = (10, 8)
# FIGSIZE = (6, 4)
FIGSIZE = (6, 3.1)
DPI = 200

# plt.style.use('ggplot')


def read_instance_df_xyv(file_name: str | Path):
    inp = pd.read_csv(file_name, sep=";", header=None, names=["x", "y", "v"])
    node_costs = np.asarray(inp.v)
    xy_points = np.asarray([inp.x, inp.y]).T
    assert xy_points.shape[1] == 2

    # compute distance matrix
    D = xy_points[:, None] - xy_points[None, :]
    D = np.sqrt((D**2).sum(-1))
    # round to nearest integer
    D = np.floor(D + 0.5).astype(np.int64)

    return inp, (node_costs, xy_points, D)


def distances(D: np.ndarray, choices: list | np.ndarray):
    dist = 0
    for a, b in zip(choices, choices[1:]):
        dist += D[a, b]
    if choices[-1] != choices[0]:
        # corner case: the looping pair is already counted
        dist += D[choices[-1], choices[0]]
    return dist


def objective(D: np.ndarray, node_costs, choices: list | np.ndarray):
    Z_cost = node_costs[choices].sum()
    Z_dist = distances(D, choices)
    return int(Z_cost), int(Z_dist)


# %%
instances = {}
for each in ["TSPA", "TSPB"]:
    _fn = f"data/instances/{each}.csv"
    # CSV as DataFrame
    # node_costs, xy_points, D
    instances[each] = read_instance_df_xyv(_fn)


# %%
# results = list(Path("data/results/").rglob("*.txt"))
# results = list(Path("data/results-LS/").rglob("*.txt"))
results = list(Path("data/results-6/").rglob("*.txt"))
# print(*[p.as_posix() for p in results], sep='\n')
# print(*[p.stem.split("_") for p in results], sep="\n")


# %%

# PLOTOUT = Path("data/plots/")
PLOTOUT = Path("data/plots-ls-variations/")
PLOTOUT.mkdir(exist_ok=True)

DISPLAY_TITLE = {
    "iterative_LS": "Iterative Local Search",
    "multiplestart_LS": "Multiple Start Local Search",
    "largescale": "Large Neighbourhood Local Search",
    "largescale2": "Large Neighbourhood Search (no LS)",
}

# %%


def plot_tsp(df, seq=None, title="TSP", show=True):
    fig, ax = plt.subplots(figsize=FIGSIZE)

    scatter = ax.scatter(
        df["x"],
        df["y"],
        c=df.v,
        cmap="viridis",
        s=30,
        # s=80,
        # edgecolor='black',
        alpha=1,  # 0.75,
    )
    fcolors = df.v
    if seq is not None:
        fcolors = [("face" if i in seq else None) for i, c in enumerate(fcolors)]
        b = np.zeros(df.shape[0], dtype=bool)
        for i in seq:
            b[i] = True
        _scatter2 = ax.scatter(
            df["x"][b],
            df["y"][b],
            c=None,
            # cmap="viridis",
            s=30,
            # s=80,
            edgecolor="black",
            alpha=1,
        )

    cbar = plt.colorbar(scatter)
    cbar.set_label("Node Cost")

    if seq is not None:
        path_x = df.loc[seq, "x"].tolist()  # tolist() fixes some weird issues
        path_y = df.loc[seq, "y"].tolist()
        ax.plot(
            path_x,
            path_y,
            color="black",
            linestyle="-",
            linewidth=1,
            # marker="x",
            marker=",",
            alpha=0.5,
        )
        # marker=',', alpha=0.5)

        # close the loop
        ax.plot(
            [path_x[-1], path_x[0]],
            [path_y[-1], path_y[0]],
            color="black",
            linestyle="-",
            linewidth=1,
            # marker="x",
            marker=",",
            alpha=0.5,
        )

    # equal scaling for both axes
    plt.gca().set_aspect("equal", adjustable="box")

    ax.set_xlabel("X Coordinate")
    ax.set_ylabel("Y Coordinate")
    ax.set_title(title)
    # plt.grid(True)
    # plt.grid(False)
    ax.set_facecolor("white")
    if show:
        plt.show()
    # return fig, ax


# %%
def main():
    for each_file in results:
        print(each_file.stem)
        try:
            # inst, method, baw = each_file.stem.split("_")
            inst, steepest, initial, nodesedges, baw = each_file.stem.split("_")
            if steepest == "iterative":
                steepest = "Iterative Local Search"
            if steepest == "multiplestart":
                steepest = "Multiple Start Local Search"
            if steepest == "largescale":
                steepest = "Large Neighbourhood Local Search"
            if steepest == "largescale2":
                steepest = "Large Neighbourhood Search (no LS)"
        except ValueError:
            print(f"\n\n  bad name: {each_file.stem.split('_')}\n\n")
            # break
            continue
        if baw != "best":
            continue
        inst = inst.removesuffix(".csv")

        sol = [int(ln) for ln in each_file.read_text().splitlines()]
        assert len(sol) == len(set(sol))
        df, (node_costs, xy_points, D) = instances[inst]
        print(len(sol), "nodes")
        zc, zd = objective(D, node_costs, sol)
        print(f"score {zd+zc}  ({zc} + {zd})")
        assert len(sol) == (len(node_costs) + 1) // 2, "number of selected nodes"
        if sol[0] != sol[1]:
            assert len(set(sol)) == len(sol), "unique indices"
        else:
            assert len(set(sol[1:])) == len(sol[1:]), "unique indices in cycle"

        # em = DISPLAY_TITLE[method.lower()]
        plot_tsp(
            df,
            sol,
            # title=f"{em} - {baw} {inst} solution\nscore {zd+zc}",
            title=f"{steepest}\nBest solution",
            show=False,
        )
        # plt.savefig(PLOTOUT / f"{inst}_{baw}_{method}.png", dpi=DPI)
        fname = PLOTOUT / each_file.with_suffix(".png").name
        plt.savefig(fname, dpi=DPI)
        # plt.show()

        print("")
        s = r"""
\begin{figure}[H] \centering
    \includegraphics[width=\textwidth]{res/plot}
\end{figure}
"""
        # % \caption{Instance ``A''}
        # % \label{fig:instanceA}
        s = s.replace("res/plot", f"plots-ls-variations/{fname.stem}")
        print(s.strip())

        print("[", end="")
        print(*sol, sep=", ", end="]\n")
        print("")


if __name__ == "__main__":
    main()
