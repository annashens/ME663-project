# ME663-project

## How to Run

```bash
./Allrun
```

You'll be prompted to enter a case index. The script will:
1. **Compile** the Fortran code using the Makefile
2. **Run** the simulation
3. **Organize results** into `case_data/{case_index}/` directory

### Manual Steps:
```bash
make                    # Compile the code → creates 'cavity' executable
./cavity                # Run the simulation (will prompt for parameters)
make clean              # Clean up
```

## Files
| File | Purpose |
|------|---------|
| **main.f** | Main entry point; selects between solvers (SIMPLE or FS) |
| **global_data.f** | Global variables and parameters |
| **grid_init.f** | Grid init |
| **boundary.f** | Boundary condition handling |
| **simple.f** / **simple_calc.f** | SIMPLE algorithm implementation |
| **fs.f** / **fs_calc.f** | Fractional step algorithm implementation |
| **fs_sor.f** | SOR pressure solver for fractional step |
| **vanka.f** | Empty placeholder file |


## Key Parameters
In `main.f`, you can modify:
- **RE** – Reynolds number (default: 1000)
- **MAXIT** – Max iterations
- **SOLVER_NAME** – Choice between `'fs'` (Fractional Step) or `'simple'` (SIMPLE)
- **SCHEME_NAME** – Convection scheme (e.g., `'quick'`)
- **DT** / **NTRANS** – Time step and number of time steps (for FS)
- **URFP, URFU, URFV** – Under-relaxation factors