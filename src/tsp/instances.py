from pathlib import Path

import numpy as np
import pandas as pd


def read_instance_xyv_df(file_name: str | Path):
    inp = pd.read_csv(file_name, sep=";", header=None, names=["x", "y", "v"])
    # inp.columns = ["x", "y", "v"]
    # assert inp.shape == (200, 3)
    # print(inp.describe())
    return inp


def read_instance_xyv(
    file_name: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    inp = read_instance_xyv_df(file_name)
    node_costs = np.asarray(inp.c)
    inp_points = np.asarray([inp.x, inp.y]).T
    assert inp_points.shape[1] == 2
    # TODO: check int64

    # compute distance matrix
    D = inp_points[:, None] - inp_points[None, :]
    D = np.sqrt((D**2).sum(-1))
    # round to nearest integer
    D = np.floor(D + 0.5).astype(np.int64)

    return node_costs, inp_points, D
