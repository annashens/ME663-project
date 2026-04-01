import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import re
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from IPython.display import HTML
from matplotlib.animation import FuncAnimation
import os
from matplotlib.animation import FFMpegWriter
class CavityCase:
    root_folder='../case_data/'
    def __init__(self, case_number):

        self.case_folder = self.root_folder+str(case_number)+'/'
        result_path = self.case_folder+str(case_number)+'_result.dat'
        self.ds = self._load_results(self.case_folder+str(case_number)+'_result.dat')

        self.res_df, self.meta = self._parse_run_log(self.case_folder+str(case_number)+'_run.log')
        self.solution_method = self.meta["solution_method"]
        self.N = self.meta["N"]
        self.RE = self.meta["RE"]
        self.clock_time = self.meta["clock_time"]
        self.IT= self.meta["IT"]
        if self.meta["solution_method"] == "fs":
            self.solution_method='EFSM'
            self.time_scheme = self.meta["time_scheme"]
            self.dt = self.meta["DT"]
            self.snapshots, self.snap_times = self._load_snapshots()
        elif self.meta["solution_method"] == "simple":
            self.solution_method="SIMPLE"
        # Extract core fields (assumes names U, V)
        self.U = self.ds["U"].values
        self.V = self.ds["V"].values

        self.x = self.ds["x"].values
        self.y = self.ds["y"].values
        self.X, self.Y = np.meshgrid(self.x, self.y)
        
        # Derived
        self.velocity_mag = np.sqrt(self.U**2 + self.V**2)
        self._psi = self.compute_streamfunction()
    def _parse_run_log(self, log_path):
        with open(log_path) as f:
            lines = f.readlines()        # --- Residual table ---
        rows = [
            [int(p[0]), int(p[1]), float(p[2]), float(p[3]), float(p[4])]
            for line in lines[1:]
            if re.match(r"\s*\d+", line)
            for p in [line.split()]
            if len(p) >= 5
        ]
        df = pd.DataFrame(rows, columns=["IT", "N", "RESORM", "RESORU", "RESORV"])
        # --- Metadata ---ff
        text = "".join(lines)
        meta = {}

        def grab(pattern, cast=float, key=None):
            m = re.search(pattern, text)
            if m:
                meta[key] = cast(m.group(1))
        grab(r"Solver:\s*(\w+)", str, "solution_method")
        grab(r"scheme:\s*(\w+)", str, "conv_scheme")
        grab(r"RE=\s*([0-9.E+-]+)", float, "RE")
        grab(r"N=\s*(\d+)", int, "N")
        grab(r"Selected Time scheme:\s*(\d+)", int, "time_scheme")
        grab(r"DT=\s*([0-9.E+-]+)", float, "DT")
        grab(r"URFU=\s*([0-9.E+-]+)", float, "URFU")
        grab(r"Wall-clock time.*:\s*([0-9.E+-]+)", float, "clock_time")
        meta["IT"] = df["IT"].max()
        return df, meta
    def _load_results(self, result_path):
        with open(result_path, "r") as f:
            lines = f.readlines()
        # --- Get variable names ---
        variable_line = next(line for line in lines if "VARIABLES" in line)
        columns = re.findall(r'"(.*?)"', variable_line)

        # --- Get zone info ---
        zone_line = [line for line in lines if "ZONE" in line][0]
        I = int(re.search(r"I=\s*(\d+)", zone_line).group(1))
        J = int(re.search(r"J=\s*(\d+)", zone_line).group(1))

        # --- Find data start ---
        data_start = 0
        for i, line in enumerate(lines):
            if re.match(r"\s*[-+0-9.]", line):
                data_start = i
                break

        # --- Load data ---
        data = np.loadtxt(result_path, skiprows=data_start)

        # Reshape to (J, I, num_vars)
        data = data.reshape(J, I, len(columns))

        X = data[:, :, 0]  # assuming first column is X
        Y = data[:, :, 1]  # assuming second column is Y
        # --- Extract other variables ---
        data_vars = {}
        for k, col in enumerate(columns):
            data_vars[col] = (("y", "x"), data[:, :, k])
        return xr.Dataset(
            data_vars,
            coords={
                "x": ("x", X[0, :]),  # x varies along columns
                "y": ("y", Y[:, 0])   # y varies along rows
            }
        )
    def _load_snapshots(self, max_snaps=None):
        snap_dir = self.case_folder + "snaps/"
        if not os.path.exists(snap_dir):
            snapshots, snap_times = [], []
            return
        files = [
            f for f in os.listdir(snap_dir)
            if re.match(r"_snap_([0-9.E+-]+)\.dat", f)
        ]

        # Sort by timestamp
        files = sorted(
            files,
            key=lambda f: float(re.search(r"_snap_([0-9.E+-]+)", f).group(1))
        )

        # Downsample to max_snaps
        if max_snaps is not None and len(files) > max_snaps:
            idx = np.linspace(0, len(files) - 1, max_snaps, dtype=int)
            files = [files[i] for i in idx]
        snapshots = []
        snap_times = []

        for f in files:
            t = float(re.search(r"_snap_([0-9.E+-]+)", f).group(1))
            path = snap_dir + f

            ds = self._load_results(path)
            snapshots.append(ds)
            snap_times.append(t)
        return snapshots, snap_times
    def plot_sl_and_vel(self, step=3, density=1, scale=2, cmap="plasma"):
        fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))

        # --- Contour + quiver ---
        contour = axs[1].contourf(
            self.X, self.Y, self.velocity_mag, levels=10, cmap=cmap
        )

        divider = make_axes_locatable(axs[1])
        cax = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(contour, label="Velocity Magnitude", cax=cax, shrink=0.8)

        axs[1].quiver(
            self.X[::step, ::step],
            self.Y[::step, ::step],
            self.U[::step, ::step],
            self.V[::step, ::step],
            color="white",
            scale=scale
        )

        # --- Streamlines ---
        axs[0].streamplot(
            self.X, self.Y, self.U, self.V,
            color='k',
            density=density,
            linewidth=0.8
        )
        # --- Formatting ---
        axs[1].set_title("Velocity Field")
        axs[0].set_title("Velocity Streamlines")

        for ax in axs:
            ax.set_xlim(self.x.min(), self.x.max())
            ax.set_ylim(self.y.min(), self.y.max())
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_aspect('equal', adjustable='box')
        fig.tight_layout()
        return fig, axs
    def animate(self, interval=300, step=3, density=1, scale=2, cmap="plasma", last=5, saveFig = False, title=None):
        if not hasattr(self, "snapshots") or len(self.snapshots) == 0:
            raise ValueError("No snapshots loaded.")

        fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
        snap_times = self.snap_times[:last]
        def update(i):
            for ax in axs:
                ax.clear()

            ds = self.snapshots[i]
            U = ds["U"].values
            V = ds["V"].values
            vel = np.sqrt(U**2 + V**2)

            # --- Right: contour + quiver ---
            contour = axs[1].contourf(self.X, self.Y, vel, levels=10, cmap=cmap)

            axs[1].quiver(
                self.X[::step, ::step],
                self.Y[::step, ::step],
                U[::step, ::step],
                V[::step, ::step],
                color="white",
                scale=scale
            )

            # --- Left: streamlines ---
            axs[0].streamplot(
                self.X, self.Y, U, V,
                color='k',
                density=density,
                linewidth=0.8
            )

            # Titles
            t = self.snap_times[i]
            axs[0].set_title(f"Streamlines (t={t:.0f})")
            axs[1].set_title(f"Velocity Field (t={t:.0f})")
            if title is not None:
                fig.suptitle(title)
            # Formatting
            for ax in axs:
                ax.set_xlim(self.x.min(), self.x.max())
                ax.set_ylim(self.y.min(), self.y.max())
                ax.set_xlabel("X")
                ax.set_ylabel("Y")
                ax.set_aspect('equal', adjustable='box')

        anim = FuncAnimation(
            fig, update,
            frames=len(snap_times),
            interval=interval
        )

        plt.close(fig)  # prevent duplicate static plot
        if saveFig:
            # anim.save(self.case_folder+"animation.gif", writer="pillow", fps=1000/interval)

            fps = 1000 / interval  # keep your frame rate calculation
            writer = FFMpegWriter(fps=fps)
            anim.save(self.case_folder + "animation.mp4", writer=writer)
        return HTML(anim.to_jshtml())
    def v_horizontal_slice(self, x=0.5):
        return self.ds.V.sel(x=x, method='nearest')
    def u_vertical_slice(self,y=0.5):
        return self.ds.U.sel(y=y, method='nearest')
    # --- Primary vortex (global min of ψ) ---
    # --- Compute streamfunction ψ ---
    # def _compute_streamfunction(self):
    #     """
    #     Compute streamfunction ψ from velocity field (U, V) by simple numerical integration.
    #     Returns a 2D array of shape (ny, nx).
    #     """
    #     ny, nx = self.U.shape
    #     dx = self.x[1] - self.x[0]
    #     dy = self.y[1] - self.y[0]

    #     # Initialize ψ = 0 at bottom-left corner
    #     psi = np.zeros_like(self.U)

    #     # Integrate U along y at first column
    #     for j in range(1, ny):
    #         psi[j, 0] = psi[j-1, 0] + self.U[j-1, 0] * dy

    #     # Integrate V along x for each row
    #     for j in range(ny):
    #         for i in range(1, nx):
    #             psi[j, i] = psi[j, i-1] - self.V[j, i-1] * dx

    #     return psi
    def compute_streamfunction(self, time=None):
        if time is None:
            ds=self.ds
        else:
            ds = self.snapshots[self.snap_times.index(time)]
        ny, nx = ds.U.shape
        dx = ds.x[1] - ds.x[0]
        dy = ds.y[1] - ds.y[0]

        # Initialize ψ = 0 at bottom-left corner
        psi = np.zeros_like(ds.U)

        # Integrate U along y at first column
        for j in range(1, ny):
            psi[j, 0] = psi[j-1, 0] + ds.U[j-1, 0] * dy

        # Integrate V along x for each row
        for j in range(ny):
            for i in range(1, nx):
                psi[j, i] = psi[j, i-1] - ds.V[j, i-1] * dx

        return psi
    @property
    def psi(self): 
        return self._psi
    @property
    def primary_vortex_location(self):
        idx = np.unravel_index(np.argmin(self._psi), self._psi.shape)
        return self.x[idx[1]], self.y[idx[0]]
    @property
    def primary_vortex_strength(self):
        return self._psi.min()

    @property
    def secondary_vortex_location(self):
        idx = np.unravel_index(np.argmax(self._psi), self._psi.shape)
        return self.x[idx[1]], self.y[idx[0]]
    @property
    def secondary_vortex_strength(self):
        return self._psi.max()
    @property
    def vortex_strengths(self):
        primary_df = pd.DataFrame({
            "vortex": ["primary"],
            "x": [self.primary_vortex_location[0]],
            "y": [self.primary_vortex_location[1]],
            "psi": [self.primary_vortex_strength]
        })

        secondary_df = pd.DataFrame({
            "vortex": ["secondary"],
            "x": [self.secondary_vortex_location[0]],
            "y": [self.secondary_vortex_location[1]],
            "psi": [self.secondary_vortex_strength]
        })
        return primary_df, secondary_df
    def plot_psi(self, cmap="viridis", markersize=10):
        fig, ax = plt.subplots(1, 1, figsize=(5, 4.5))
        cf = ax.contourf(self.X, self.Y, self.psi, levels=20, cmap=cmap)
        fig.colorbar(cf, label=r"$\psi$", shrink=0.8)
        ax.set_title("Streamfunction")
        ax.set_ylabel("Y")
        ax.set_xlabel("X")
        ax.set_aspect('equal', adjustable='box')
        primary_x, primary_y = self.primary_vortex_location
        secondary_x, secondary_y = self.secondary_vortex_location
        # Plot primary vortex
        ax.scatter(primary_x, primary_y, color='white', marker='o', s=markersize)  # s=50 for marker size
        ax.annotate(rf'$\psi_{{min}}={self.primary_vortex_strength:.4f}$'+f'\n({primary_x:.2f}, {primary_y:.2f})', 
                    xy=(primary_x, primary_y), 
                    xytext=(primary_x + 0.0, primary_y + 0.05),
                    fontsize=9,
                    ha='center',
                    color='white')
        # Plot secondary vortex
        ax.scatter(secondary_x, secondary_y, color='black', marker='o', s=markersize)
        ax.annotate(rf'$\psi_{{max}}={self.secondary_vortex_strength:.4f}$'+f'\n({secondary_x:.2f}, {secondary_y:.2f})', 
                    xy=(secondary_x, secondary_y), 
                    xytext=(secondary_x - 0.05, secondary_y + 0.05), ha='center',
                    fontsize=9,
                    color='black')
        return fig,ax
class CaseManager:
    ghia_psi_min, ghia_psi_max = -0.117929, 0.00175102
    def __init__(self, case_numbers, root_folder="../case_data/", RE=None, solution_method=None, N=None,conv_scheme=None, time_scheme=None, label=None):
        self.cases = {}
        for num in case_numbers:
            case = CavityCase(num)  # Your existing case class
            self.cases[num] = case
        self.RE = RE
        self.solution_method = solution_method
        self.N = N
        self.conv_scheme = conv_scheme
        self.time_scheme = time_scheme
        self.label = label
    def collect_vortex_strengths(self):
        primary_records = []
        secondary_records = []
        for case_id, case in self.cases.items():
            # Get the two DataFrames from vortex_strengths
            primary_df, secondary_df = case.vortex_strengths
            # Add metadata to each
            for df in [primary_df, secondary_df]:
                df["case"] = case_id
                df["RE"] = case.RE
                df["scheme"] = case.solution_method
                df["N"] = case.N
                # Add any other relevant metadata here
            primary_records.append(primary_df)
            secondary_records.append(secondary_df)

        # Concatenate separately
        primary_all = pd.concat(primary_records, ignore_index=True)
        secondary_all = pd.concat(secondary_records, ignore_index=True)

        return primary_all, secondary_all
    def compare_convergence(self):
        # For example, combine convergence residuals from all cases
        import pandas as pd
        records = []
        for case_id, case in self.cases.items():
            df = case.res_df.copy()
            df["case"] = case_id
            df["RE"] = case.RE
            records.append(df)
        return pd.concat(records, ignore_index=True)

    def compare_meta(self):
        meta_dicts = []
        for case_id, case in self.cases.items():
            meta_dicts.append(case.meta)
        meta_df = pd.DataFrame(meta_dicts, index=self.cases.keys())
        # get convergence dataframe
        conv_df = self.compare_convergence()

        if conv_df.empty:
            print("Convergence dataframe is empty!")
            return meta_df

        conv_df['case'] = conv_df['case'].astype(str)

        # aggregate convergence info per case
        conv_summary = conv_df.groupby('case').agg(
            sum_N=('N', 'sum'),
            total_outer_iterations=('IT', 'max')
        ).reset_index()

        # scale inner iterations by logging interval
        LOG_INTERVAL = 20
        conv_summary['total_inner_iterations'] = conv_summary['sum_N'] * LOG_INTERVAL

        # average inner iterations per outer iteration
        conv_summary['avg_inner_per_outer'] = (
            conv_summary['total_inner_iterations'] / conv_summary['total_outer_iterations']
        )
        conv_summary['case'] = conv_summary['case'].astype(int)
        conv_summary = conv_summary.set_index('case')
        # print(conv_summary)

        # join with meta dataframe

        return pd.concat([meta_df, conv_summary], axis=1)
    def compare_time(self,time=2.5):
        it_lst = 2.5/self.compare_meta()["DT"]
        snapshot_lst=[]
        dict_lst = []
        for case_id, case in self.cases.items():
            target = it_lst[case_id]

            i, closest_val = min(
                enumerate(case.snap_times),
                key=lambda x: abs(x[1] - target)
            )

            snapshot_lst.append(case.snapshots[i])
            case.compute_streamfunction(time=closest_val).min()
            dict_lst.append({
                "case": case_id,
                "DT": case.meta["DT"],
                "IT": int(target),
                "psi_min": case.compute_streamfunction(time=closest_val).min(),
            })
        return pd.DataFrame(dict_lst)