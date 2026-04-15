import pandas as pd
import numpy as np
import re
from scipy.interpolate import griddata
def parseGhia(path, vel):
    vel=vel.upper()
    dir = {
        "U":"Y", "V":"X"
    }
    df = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None
    )
    # Set column names: first column is 'y', remaining are Re values
    re_values = [100, 400, 1000, 3200, 5000, 7500, 10000]
    df.columns = [dir[vel]] + re_values
    # Set 'y' as the index
    # df.set_index("y", inplace=True)
    return df

def parseResults(filename):
    # Read file
    with open(filename, "r") as f:
        lines = f.readlines()
        # Extract RE, MAXIT, URFU from the first line
    header_line = [line for line in lines if "RE=" in line][0]
    RE = float(re.search(r"RE=\s*([\d.]+)", header_line).group(1))
    MAXIT = int(re.search(r"MAXIT=\s*(\d+)", header_line).group(1))
    URFU = float(re.search(r"URFU=\s*([\d.]+)", header_line).group(1))
    zone_line = [line for line in lines if "ZONE" in line][0]
    I = int(re.search(r"I=\s*(\d+)", zone_line).group(1))
    J = int(re.search(r"J=\s*(\d+)", zone_line).group(1))
    metadata={
        "RE": RE,
        "MAXIT": MAXIT,
        "URFU": URFU,
        "I": I,
        "J": J
    }
    
    data_start = 0
    for i, line in enumerate(lines):
        if re.match(r"\s*[-+0-9.]", line):
            data_start = i
            break
    df = pd.read_csv(
        filename,
        sep=r"\s+",
        skiprows=data_start,
        names=["X", "Y", "U", "V", "P"],
        engine="python"
    )

    return metadata, df
def getValues(metadata, df):
    I = metadata["I"]
    J = metadata["J"]
    X = df["X"].values.reshape(J, I)
    Y = df["Y"].values.reshape(J, I)
    U = df["U"].values.reshape(J, I)
    V = df["V"].values.reshape(J, I)
    velocity_mag = np.sqrt(U**2 + V**2)
    return X, Y, U, V, velocity_mag

def cutline(df, *, x=None, y=None, columns=None, method="linear"):

    if (x is None and y is None) or (x is not None and y is not None):
        raise ValueError("Provide exactly one of x or y.")

    if columns is None:
        raise ValueError("You must provide columns to interpolate.")

    # Vertical cutline (x = constant)
    if x is not None:
        varying = np.sort(df["Y"].unique())
        interp_points = np.column_stack([
            np.full(varying.shape, x),
            varying
        ])
        out = pd.DataFrame({
            "X": np.full(varying.shape, x),
            "Y": varying
        })

    # Horizontal cutline (y = constant)
    else:
        varying = np.sort(df["X"].unique())
        interp_points = np.column_stack([
            varying,
            np.full(varying.shape, y)
        ])
        out = pd.DataFrame({
            "X": varying,
            "Y": np.full(varying.shape, y)
        })

    # Interpolate
    for col in columns:
        out[col] = griddata(
            df[["X", "Y"]].values,
            df[col].values,
            interp_points,
            method=method
        )

    return out

def parseResiduals(filename):
    pattern = re.compile(
        r"IT=\s*(\d+)\s+RESU=\s*([-\d.Ee+]+)\s+RESV=\s*([-\d.Ee+]+)\s+RESM=\s*([-\d.Ee+]+)"
    )
    data = []
    with open(filename, "r") as f:
        for line in f:
            match = pattern.search(line)  # use search instead of match
            if match:
                it, resu, resv, resm = match.groups()
                data.append([int(it), float(resu), float(resv), float(resm)])
                # Optional: debug print
                # print(it, resu, resv, resm)
    df = pd.DataFrame(data, columns=["IT", "RESU", "RESV", "RESM"])
    return df