import numpy as np
import xarray as xr
import re
from scipy.interpolate import griddata


class CFDCase:
    """
    Represents a full CFD dataset loaded from results.dat.
    Provides analysis utilities directly on the dataset.
    """

    def __init__(self, filename):
        self.filename = filename
        self.ds = self._load_results()

    # ------------------------------------------------------------
    # INTERNAL PARSER
    # ------------------------------------------------------------
    def _load_results(self):
        with open(self.filename, "r") as f:
            lines = f.readlines()
        print(lines)
        header_line = next(line for line in lines if "RE=" in line)
        zone_line = next(line for line in lines if "ZONE" in line)

        # RE = float(re.search(r"RE=\s*([\d.]+)", header_line).group(1))
        # MAXIT = int(re.search(r"MAXIT=\s*(\d+)", header_line).group(1))
        # URFU = float(re.search(r"URFU=\s*([\d.]+)", header_line).group(1))
        I = int(re.search(r"I=\s*(\d+)", zone_line).group(1))
        J = int(re.search(r"J=\s*(\d+)", zone_line).group(1))

        # Find start of numeric block
        for i, line in enumerate(lines):
            if re.match(r"\s*[-+0-9.]", line):
                data_start = i
                break

        data = np.loadtxt(self.filename, skiprows=data_start)
        X, Y, U, V, P = (data[:, k].reshape(J, I) for k in range(5))

        ds = xr.Dataset(
            {
                "U": (("y", "x"), U),
                "V": (("y", "x"), V),
                "P": (("y", "x"), P),
            },
            coords={
                "x": X[0, :],
                "y": Y[:, 0],
            },
            # attrs={
            #     "RE": RE,
            #     "MAXIT": MAXIT,
            #     "URFU": URFU,
            #     "I": I,
            #     "J": J,
            # }
        )

        # Derived field
        ds["velocity_mag"] = np.sqrt(ds.U**2 + ds.V**2)

        return ds

    # ------------------------------------------------------------
    # CUTLINE
    # ------------------------------------------------------------
    def cutline(self, *, x=None, y=None, fields=None, method="linear"):
        if (x is None and y is None) or (x is not None and y is not None):
            raise ValueError("Provide exactly one of x or y.")

        if fields is None:
            raise ValueError("You must specify fields to interpolate.")

        X, Y = np.meshgrid(self.ds.x.values, self.ds.y.values)
        points = np.column_stack([X.ravel(), Y.ravel()])

        if x is not None:
            varying = self.ds.y.values
            interp_points = np.column_stack([np.full_like(varying, x), varying])
            coords = {"y": varying}
        else:
            varying = self.ds.x.values
            interp_points = np.column_stack([varying, np.full_like(varying, y)])
            coords = {"x": varying}

        out = {}
        for field in fields:
            values = self.ds[field].values.ravel()
            out[field] = (list(coords.keys()),
                          griddata(points, values, interp_points, method=method))

        return xr.Dataset(out, coords=coords)

    # ------------------------------------------------------------
    # RESIDUAL PARSER
    # ------------------------------------------------------------
    def load_residuals(self, filename):
        pattern = re.compile(
            r"IT=\s*(\d+)\s+RESU=\s*([-\d.Ee+]+)\s+RESV=\s*([-\d.Ee+]+)\s+RESM=\s*([-\d.Ee+]+)"
        )

        it_list, resu_list, resv_list, resm_list = [], [], [], []

        with open(filename, "r") as f:
            for line in f:
                m = pattern.search(line)
                if m:
                    it, resu, resv, resm = m.groups()
                    it_list.append(int(it))
                    resu_list.append(float(resu))
                    resv_list.append(float(resv))
                    resm_list.append(float(resm))

        return xr.Dataset(
            {
                "RESU": ("iteration", resu_list),
                "RESV": ("iteration", resv_list),
                "RESM": ("iteration", resm_list),
            },
            coords={"iteration": it_list}
        )
