from symbols import i, j, u, v
import sympy as sp
import re
class FaceFlux:
    _relative_equations = {
        "U": {
            'e':  (1/2)*(u[i, j]+u[i+1, j]),
            'w': (1/2)*(u[i-1, j] + u[i, j]),
            'n': (1/2)*(v[i, j] + v[i+1, j]),
            's': (1/2)*(v[i, j-1] + v[i+1, j-1]),
        },
        "V": {
            'e': (1/2)*(u[i, j] + u[i, j+1]),  
            'w': (1/2)*(u[i-1, j] + u[i-1, j+1]),
            'n': (1/2)*(v[i, j] + v[i, j+1]),
            's': (1/2)*(v[i, j] + v[i, j-1]),
        }
    }

    def __init__(self, kind="U"):
        self.kind = kind
        self.dx, self.dy = sp.symbols('dx dy', positive=True)

        # Symbolic flux variables 
        self.e = [sp.Symbol("Fe")]
        self.w = [sp.Symbol("Fw")]
        self.n = [sp.Symbol("Fn")]
        self.s = [sp.Symbol("Fs")]

        self.faces = {'e': self.e, 'w': self.w, 'n': self.n, 's': self.s}

    def build_equations(self, i_loc, j_loc, id_fort=None):
        subs = {i: i_loc, j: j_loc}
        eqs = self._relative_equations[self.kind]

        self.e_eq = eqs['e'].subs(subs)*self.dy
        self.w_eq = eqs['w'].subs(subs)*self.dy
        self.n_eq = eqs['n'].subs(subs)*self.dx
        self.s_eq = eqs['s'].subs(subs)*self.dx

    def fortran_subs(self, id_fort):
        return {
            self.e[0]: sp.Symbol(f"FE{self.kind}({id_fort})"),
            self.w[0]: sp.Symbol(f"FW{self.kind}({id_fort})"),
            self.n[0]: sp.Symbol(f"FN{self.kind}({id_fort})"),
            self.s[0]: sp.Symbol(f"FS{self.kind}({id_fort})"),
        }

class UDS:
    def build(self, phi, flux_cv, i, j):

        phiP = phi[i,j]
        phiE = phi[i+1,j]
        phiW = phi[i-1,j]
        phiN = phi[i,j+1]
        phiS = phi[i,j-1]

        Fe = flux_cv.e[0]
        Fw = flux_cv.w[0]
        Fn = flux_cv.n[0]
        Fs = flux_cv.s[0]

        return (
            phiP*sp.Max(Fe,0.0) - phiE*sp.Max(-Fe,0.0)
            - phiW*sp.Max(Fw,0.0) + phiP*sp.Max(-Fw,0.0)
            + phiP*sp.Max(Fn,0.0) - phiN*sp.Max(-Fn,0.0)
            - phiS*sp.Max(Fs,0.0) + phiP*sp.Max(-Fs,0.0)
        )
class QUICK:
    def build(self, phi, flux_cv, i, j):
        phiP = phi[i,j]
        phiE = phi[i+1,j]
        phiEE = phi[i+2,j]
        
        phiW = phi[i-1,j]
        phiWW = phi[i-2,j]
        
        phiN = phi[i,j+1]
        phiNN = phi[i,j+2]
        
        phiS = phi[i,j-1]
        phiSS = phi[i,j-2]

        Fe = flux_cv.e[0]
        Fw = flux_cv.w[0]
        Fn = flux_cv.n[0]
        Fs = flux_cv.s[0]

        return (
              (-1/8*phiW  + 3/4*phiP + 3/8*phiE)*sp.Max(Fe,0.0) - (-1/8*phiEE + 3/4*phiE + 3/8*phiP)*sp.Max(-Fe,0.0)
            - (-1/8*phiWW + 3/4*phiW + 3/8*phiP)*sp.Max(Fw,0.0) + (-1/8*phiE + 3/4*phiP + 3/8*phiW )*sp.Max(-Fw,0.0)
            + (-1/8*phiS  + 3/4*phiP + 3/8*phiN)*sp.Max(Fn,0.0) - (-1/8*phiNN + 3/4*phiN + 3/8*phiP)*sp.Max(-Fn,0.0) 
            - (-1/8*phiSS + 3/4*phiS + 3/8*phiP)*sp.Max(Fs,0.0) + (-1/8*phiN + 3/4*phiP + 3/8*phiS)*sp.Max(-Fs,0.0)
        )
class CDS:
    def build(self, phi, flux_cv, i, j):
        phiP = phi[i,j]
        phiE = phi[i+1,j]
        phiW = phi[i-1,j]
        phiN = phi[i,j+1]
        phiS = phi[i,j-1]

        Fe = flux_cv.e[0]
        Fw = flux_cv.w[0]
        Fn = flux_cv.n[0]
        Fs = flux_cv.s[0]

        return (
            (Fe/2)*(phiP + phiE)
            - (Fw/2)*(phiW + phiP)
            + (Fn/2)*(phiP + phiN)
            - (Fs/2)*(phiS + phiP)
        )

class HYBRID:
    def build(self, phi, flux_cv, i, j, diff_scheme=None, nu=None, dx=None, dy=None):
        # Cell values
        phiP = phi[i, j]
        phiE = phi[i+1, j]
        phiW = phi[i-1, j]
        phiN = phi[i, j+1]
        phiS = phi[i, j-1]

        # Fluxes
        Fe = flux_cv.e[0]
        Fw = flux_cv.w[0]
        Fn = flux_cv.n[0]
        Fs = flux_cv.s[0]

        # Diffusion contributions
        diff_expr = diff_scheme.build(phi, nu, dx, dy, i, j)
        Dw = sp.diff(diff_expr, phiW)
        De = sp.diff(diff_expr, phiE)
        Ds = sp.diff(diff_expr, phiS)
        Dn = sp.diff(diff_expr, phiN)

        return (
              phiP*sp.Max(Fe, De + Fe/2,0.0) - phiE*sp.Max(-Fe, De - Fe/2,0.0) 
            - phiW*sp.Max(Fw, Dw + Fw/2,0.0) + phiP*sp.Max(-Fw, Dw -Fw/2,0.0)
            + phiP*sp.Max(Fn, Dn + Fn/2,0.0) - phiN*sp.Max(-Fn, Dn - Fn/2,0.0)
            - phiS*sp.Max(Fs, Ds + Fs/2,0.0) + phiP*sp.Max(-Fs, Dw -Fs/2,0.0)
        )
    
class SecondOrderDiffusion:
    def build(self, phi, nu, dx, dy,i, j):
        xDiff = nu*(phi[i+1,j] - phi[i,j])*dy/dx - nu*(phi[i,j] - phi[i-1,j])*dy/dx
        yDiff = nu*(phi[i,j+1] - phi[i,j])*dx/dy - nu*(phi[i,j] - phi[i,j-1])*dx/dy
        return xDiff + yDiff
    
class FourthOrderDiffusion:
    def build(self,phi,nu,dx,dy,i,j):
        xDiff = nu/12 * dy/dx * (-phi[i+2,j] + 16*phi[i+1,j] - 30*phi[i,j] + 16*phi[i-1,j] - phi[i-2,j])
        yDiff = nu/12 * dx/dy * (-phi[i,j+2] + 16*phi[i,j+1] - 30*phi[i,j] + 16*phi[i,j-1] - phi[i,j-2])
        return xDiff+yDiff
    
class sourceTerm:
    def build(self, kind, p, i,j, dx, dy):
        if kind == "U":
            return  (p[i,j] - p[i+1,j])*dy
        elif kind == "V":
            return (p[i,j] - p[i,j+1])*dx
        else:
            raise ValueError("kind must be 'U' or 'V'")
class MomentumFVM:
    def __init__(self, kind, offset,
                 convection_scheme=None,
                 diffusion_scheme=None,):  
        self.kind = kind
        
        di, dj = (offset,0) if kind=="U" else (0,offset)
        self.i = i + di
        self.j = j + dj
        # self.coord = i + di, j + dj
        self.coord={
            "p": (self.i, self.j),
            "e": (self.i+1, self.j),
            "ee": (self.i+2, self.j),
            "w": (self.i-1, self.j),  
            "ww": (self.i-2, self.j),          
            "n": (self.i, self.j+1),
            "nn": (self.i, self.j+2),
            "s": (self.i, self.j-1),
            "ss": (self.i, self.j-2)
        }
        self.id_fort = offset +1
        self.phi = sp.IndexedBase('phi')
        self.flux_cv = FaceFlux(self.kind)
        self.flux_cv.build_equations(self.i, self.j) 
        self.convection_scheme = convection_scheme
        self.diffusion_scheme = diffusion_scheme
        self.nu, self.dx, self.dy = sp.symbols('NU dx dy', positive=True)
        self.p = sp.IndexedBase('p')
        self.source = sourceTerm().build(self.kind, self.p, self.i, self.j, self.dx, self.dy)
        
        # convection and diffusion
        if isinstance(self.convection_scheme, HYBRID):
            self.conv = self.convection_scheme.build(
                self.phi, self.flux_cv, self.i, self.j, nu=self.nu, dx=self.dx, dy=self.dy, diff_scheme=self.diffusion_scheme
            )
            self.diff = 0
        else:
            self.conv = self.convection_scheme.build(
                self.phi, self.flux_cv, self.i, self.j
            )
            self.diff = self.diffusion_scheme.build(
                self.phi, self.nu, self.dx, self.dy, self.i, self.j
            )
        expr=sp.expand(sp.simplify(self.conv - self.diff + self.source))
        phi_terms = sorted(
            expr.atoms(sp.Indexed),
            key=lambda x: x.indices
        )
        self.eq = sp.collect(expr, phi_terms)
        self.A = self._build_coefficients()
        self.R =self._build_residual()
        # store face flux equations locally for printing
        self.F ={}
        self.F['e'] = self.flux_cv.e_eq
        self.F['w'] = self.flux_cv.w_eq
        self.F['n'] = self.flux_cv.n_eq
        self.F['s'] = self.flux_cv.s_eq
        return

    def _build_coefficients(self):
        phi = self.phi
        i=self.i
        j=self.j
        self.test="test"
        return {
            'p': sp.simplify(sp.diff(self.eq, phi[i,j])),
            'e': sp.simplify(-sp.diff(self.eq, phi[i+1,j])),
            'ee': sp.simplify(-sp.diff(self.eq, phi[i+2,j])),
            'w': sp.simplify(-sp.diff(self.eq, phi[i-1,j])),
            'ww': sp.simplify(-sp.diff(self.eq, phi[i-2,j])),
            'n': sp.simplify(-sp.diff(self.eq, phi[i,j+1])),
            'nn': sp.simplify(-sp.diff(self.eq, phi[i,j+2])),
            's': sp.simplify(-sp.diff(self.eq, phi[i,j-1])),
            'ss': sp.simplify(-sp.diff(self.eq, phi[i,j-2])),
        }
        
    def _build_residual(self): 
        i, j = self.i, self.j
        phi = self.phi
        R = -(self.eq - self.A['p']*phi[i,j])
        # Clean up negatives nicely
        R = sp.factor_terms(R)
        R = sp.simplify(R)
        return R

    def __repr__(self):
        label_str=f"{self.kind}-CV at {self.coord["p"]}\n"
        velocityFlux_str = (f"Face Flux: \n Fe={self.F['e']}, \n Fw={self.F['w']}, \n Fn={self.F['n']}, \n Fs={self.F['s']})\n")
        A_str = (f"A coefficients:\n"
                 +f" AE{self.kind}{self.coord["p"]}={self.A['e']}\n"
                 +f" AEE{self.kind}{self.coord["p"]}={self.A['ee']}\n"
                 +f" AW{self.kind}{self.coord["p"]}={self.A['w']}\n"
                 +f" AWW{self.kind}{self.coord["p"]}={self.A['ww']}\n"
                 +f" AN{self.kind}{self.coord["p"]}={self.A['n']}\n"
                 +f" ANN{self.kind}{self.coord["p"]}={self.A['nn']}\n"
                 +f" AS{self.kind}{self.coord["p"]}={self.A['s']}\n"
                 +f" ASS{self.kind}{self.coord["p"]}={self.A['ss']}\n"
                 +f" AP{self.kind}{self.coord["p"]}={self.A['p']}\n")
        R_str = (f"B{self.kind}{self.coord["p"]}:\n"+
                 f"{self.R}")
        return label_str+velocityFlux_str+A_str+R_str
    def _to_fortran(self, print_code=True):

        id_fort = self.id_fort  # or offset+1 or whatever you want
        subs_flux = self.flux_cv.fortran_subs(id_fort)

        lines_F = [f"      ! Face Flux for {self.kind}-CV at ({self.i},{self.j}):"]
        lines_A = [f"      ! A coefficients for {self.kind}-CV at ({self.i},{self.j}):"]
        lines_B = [f"      ! B residual for {self.kind}-CV at ({self.i},{self.j}):"]
        
        
        lines_B.append(f"   B{self.kind}({self.id_fort})={str(self.source)}")

        # ---- Face flux definitions ----
        flux_eqs = {
            'e': self.flux_cv.e_eq,
            'w': self.flux_cv.w_eq,
            'n': self.flux_cv.n_eq,
            's': self.flux_cv.s_eq,
        }

        for d, expr in flux_eqs.items():
            fort_name = subs_flux[self.flux_cv.faces[d][0]]
            lines_F.append(f"{fort_name} = {expr}")

        lines_Ap=[]
        for d, expr in self.A.items():
            A_var = f"A{d}{self.kind}({self.id_fort})"
            A_expr_fort = expr.subs(subs_flux)            
            A_str = f"      {A_var} = {A_expr_fort}"
            if d!="p":
                lines_A.append(A_str)
                sign = "+" 
            else:
                A_str = "\n     &    + ".join(" + ".join(A_str.split(" + ")[i:i+2]) for i in range(0, len(A_str.split(" + ")), 2))
                lines_Ap.insert(0, A_str)
                sign ="-"
            lines_B.append(f"  &    {sign} {A_var}*{self.kind}{self.coord[d]}")
        lines_A.append(lines_Ap[0])
        lines_B.append("\n")
        output = ("\n      ".join(lines_F) + "\n\n" + "\n".join(lines_A) + "\n\n" + "\n   ".join(lines_B)).upper().replace("[", "(").replace("]", ")")
        output=output
        output=re.sub(r"\bNU\b", "(1/RE)", output)
        if print_code:
            print(output)
        else:
            return output
       