from __future__ import print_function
from math import exp, log, sqrt, pi, cos, sin, ceil
from datetime import datetime
from fenics import *
from dolfin import *
import numpy as np
import time
class setLoad:
    def __init__(self, loadType, fMax, T):
        self.loadType, self.fMax, self.T = loadType, fMax, T
    def tStep(self):
        if self.loadType == "Gaussian":
            freqCut = 0.1
            sigma = sqrt(-log(freqCut) / (2 * pi ** 2 * self.fMax ** 2))
            load_at_t = Expression('t <= 0.05 * T + 10*sigma && t >= 0.05 * T 1/(sqrt(2 * pi) * sigma) * exp(-pow((t-.05*T-5*sigma), 2) /(2 * pow(sigma, 2))) : 0, degree = 2, sigma = sigma, pi = pi, t = 0, T = self.T, f = self.fMax')
            self.sigma = sigma
        elif self.loadType == "Plane":
            load_at_t = Expression('cos(f * t)', degree = 2, f = self.fMax, t = 0)
        elif self.loadType == "sil_plane_sil":
            c1 = ceil((0.05*self.T*self.fMax)/(2*pi) / 2.) * 2
            c2 = ceil((0.75*self.T*self.fMax)/(2*pi) / 2.) * 2
            load_at_t = Expression('t <= c2*2*pi/f && t >= c1*2*pi/f ? sin(f*t) : 0', degree = 2, pi = pi, f = self.fMax, t = 0, c1 = c1, c2 = c2)
        return load_at_t
    def sig(self):
        return self.sigma
class AcousticSolver:
    def __init__(self, T, Z, c, fMax, rho0):
        self.T, self.Z, self.c, self.fMax, self.rho0 = T, Z, c, fMax, rho0
    def setDomain(self, type, speakerSize, lx, ly):
        self.lx, self.ly, self.speakerSize = lx, ly, speakerSize
        # Mesh
        wavelength = self.c / self.fMax
        nElx, nEly = int(20 * lx / wavelength), int(20 * ly / wavelength)
        if type == "Waveguide":
        nEly = 1
        mesh = RectangleMesh(Point(0.0), Point(lx, ly), nElx, nEly)
        # Boundary conditions
        if speakerSize < wavelength/20:
        speakerSize = wavelength/8
        tol = 1E-14
        class Wall(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary
                #return on_boundary and near(x[0], lx, tol)
        class Speaker(SubDomain):
            def inside(self, x, on_boundary):
                #return on_boundary and near(x[0], 0, tol)
                return on_boundary and (near(x[0], 0, tol) and near(x[1], ly/2, speakerSize/2))
        # Mark boundaries
        boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 1)
        bx1, bx2 = Wall(), Speaker()
        bx1.mark(boundary_markers, 1)
        bx2.mark(boundary_markers, 2)
        self.ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)
        self.mesh = mesh
        print("Rectangular Mesh with " + str(nElx) + " x " + str(nEly) + " Elements")
    def solve(self, loadType, micLoc):
        c, Z, T, ds = self.c, self.Z, self.T, self.ds
                # Previous time steps
        V = FunctionSpace(self.mesh, 'P', 1)
        p = TrialFunction(V)
        v = TestFunction(V)
        pf = Function(V)
        p_old = Function(V)
        p_d_old = Function(V)
        p_dd_old = Function(V)
        p_d_new = Function(V)
        p_dd_new = Function(V)
        # Dirichlet
        def boundary(x, on_boundary):
            return (near(x[0], self.lx, 1E-14) and near(x[1], self.ly, 1E-14))
        bc = DirichletBC(V, 0, boundary, method = "pointwise") # Dirichlet pin
        # Newmark parameters
        gamma, beta = 0.5, 0.25
        # Load defition
        load = setLoad(loadType, self.fMax, T)
        a_N = load.tStep()
        sigma = load.sig()
        self.sigma = sigma
        f = self.rho0 * a_N * v * ds(2)
        # Updateable linear system
        class linearSystem():
            def __init__(self, dt):
                self.dt = dt
            def coeffs(self, gamma, beta, Z, ds, c):
                self.Z, self.ds, self.c = Z, ds, c
                self.gamma, self.beta = gamma, beta
                self.coeff1 = 1/(beta*self.dt**2)
                self.coeff2 = (1-2*beta)/(2*beta)
                self.coeff3 = gamma/(beta*self.dt)
                return self.coeff1, self.coeff2, self.coeff3
            def functions(self, p, v, p_old, p_d_old, p_dd_old, p_d_new, p_dd_new):
                self.p = p
                self.v = v
                self.p_old = p_old
                self.p_d_old = p_d_old
                self.p_dd_old = p_dd_old
                self.p_d_new = p_d_new
                self.p_dd_new = p_dd_new
            def leftSide(self):
                m_left = self.coeff1/(self.c**2) * self.p * self.v * dx
                k_left = dot(grad(self.p), grad(self.v)) * dx
                c_left = self.coeff3/self.Z * self.p * self.v * self.ds(1)
                a = m_left + k_left + c_left
                    return a
            def mass_right(self):
                m_right = 1/(self.c**2)*(self.coeff1*(self.p_old + self.p_d_old * self.dt)+ self.coeff2 * self.p_dd_old) * self.v * dx
                    return m_right
            def damp_right(self):
                c_right = 1/self.Z * (self.p_d_old + self.dt * (1-self.gamma) *self.p_dd_old - self.coeff3 * (self.p_old + self.dt * self.p_d_old)- self.coeff2 * self.gamma * self.dt * self.p_dd_old) * self.v *self.ds(1)
                return c_right
        # Time loop
        dt0 = (0.05*T)/10
        dt1 = 10*sigma/100
        dt2 = 1/(2*self.fMax)
        dt_list = [dt0, dt1, dt2]
        t, n, t_update, nodeOut, loadIn, time_vec, sol_list = 0, 0, 0, [], [], [], []
        #vtkfile = File('acoustics/solution.pvd')
        print("Inializing solver with ", str(int((T-10*sigma-0.05*T)/dt2+110)) ," time steps")
        def convert(seconds):
            min, sec = divmod(seconds, 60)
            hour, min = divmod(min, 60)
            return hour, min, sec
        for i in range(len(dt_list)):
            # while loop before, during and after pulse
            t1 = [0, 0.05*T, 0.05*T+10*sigma]
            t2 = [0.05*T, 0.05*T+10*sigma, T]
            # Get dt
            dt = dt_list[i]
            # get a, l based on dt
            linear = linearSystem(dt)
            coeff1, coeff2, coeff3 = linear.coeffs(gamma, beta, Z, ds, c)
            linear.functions(p, v, p_old, p_d_old, p_dd_old, p_d_new, p_dd_new)
            a = linear.leftSide()
            m_right, c_right = linear.mass_right(), linear.damp_right()
            l = f - c_right + m_right # Neumann Load
            # get coeffs
            while t1[i] <= t < t2[i]:
                t = dt * n + t_update
                a_N.t = t
                A, b = assemble_system(a, l, bc) # Dirichlet pin
                solve(A, pf.vector(), b)
                #vtkfile << (pf, t)
                # Update previous solution
                p_dd_new = coeff1 * (pf - p_old - p_d_old * dt) - coeff2 * p_dd_old
                p_d_new = p_d_old + dt * ((1-gamma) * p_dd_old + gamma * p_dd_new)
                p_dd_old.assign(p_dd_new)
                p_d_old.assign(p_d_new)
                p_old.assign(pf)
                n += 1
                nodeOut.append(pf(micLoc))
                if t < .05*T + 10*sigma:
                    loadIn.append(1/(sqrt(2*pi)*sigma) *exp(-(t-.05*T-5*sigma)**2/(2*sigma**2)))
                else:
                    loadIn.append(0)
                time_vec.append(t)
                # Calculation time
                sol = ceil(t/T*100)
                time_count = time.perf_counter()
                if sol % 5 == 0:
                    if sol not in sol_list:
                        hr, min, sec = convert(int(time_count))
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        print("Current Time = " , current_time , " Completed ", "{:.0f}".format(sol), "% of Solution. Elapsed time: " ,
                        "{:.0f}".format(hr), " hr, ", "{:.0f}".format(min), " min," + "{:.0f}".format(sec) + " s")
                        sol_list.append(sol)
                        self.nodeOut, self.time_vec, self.loadIn = nodeOut, time_vec, loadIn

            t_update = t
            n = 0
    def outputSol(self):
        nodeFile = "nodeOut" + str(self.fMax) + ".txt"
        timeFile = "time" + str(self.fMax) + ".txt"
        loadFile = "load" + str(self.fMax) + ".txt"
        sigmaFile = "sigma" + str(self.fMax) + ".txt"
        fMaxFile = "fMax.txt"
        np.savetxt(nodeFile, self.nodeOut)
        np.savetxt(timeFile, self.time_vec)
        np.savetxt(loadFile, self.loadIn)
        np.savetxt(sigmaFile, np.atleast_1d(self.sigma))
        np.savetxt(fMaxFile, np.atleast_1d(self.fMax))
# Initiate Solver
a_solver = AcousticSolver(T = 1, Z = 1E6, c = 343, fMax = 1000, rho0 = 1.225)
a_solver.setDomain("N/A", 0.1, lx = 5, ly = 5)
a_solver.solve("Gaussian", micLoc = Point(.5, 2.5))
a_solver.outputSol()
