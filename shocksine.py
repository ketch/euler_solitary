#!/usr/bin/env python
# encoding: utf-8
r"""Shu-Osher problem.
   1D compressible inviscid flow (Euler equations)."""
import numpy as np
gamma = 1.4
gamma1 = gamma - 1.

def switch_to_periodic(solver,state):
    #Change to periodic BCs after initial pulse 
    from clawpack import pyclaw
    if state.t>50 and solver.bc_lower[0]==pyclaw.BC.wall:
        solver.bc_lower[0] = pyclaw.BC.periodic
        solver.bc_upper[0] = pyclaw.BC.periodic


def setup(use_petsc=False,iplot=False,htmlplot=False,outdir='./_output',solver_type='sharpclaw',kernel_language='Fortran',char_decomp=2):
    """
    Solve the Euler equations of compressible fluid dynamics.
    This example involves a shock wave impacting a sinusoidal density field.
    """
    from clawpack import riemann
    rs = riemann.euler_with_efix_1D

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)
        #solver.time_integrator = 'SSP33'
        #solver.a, solver.b, solver.c = a, b, c
        #solver.cfl_desired = 0.4
        #solver.cfl_max = 0.5
        solver.char_decomp = char_decomp
        import sharpclaw1
        solver.fmod = sharpclaw1
    elif solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)

    #solver.before_step = switch_to_periodic

    solver.bc_lower[0]=pyclaw.BC.wall
    solver.bc_upper[0]=pyclaw.BC.wall

    # Initialize domain
    xl = -10.; xr = 2000.
    mx=300*(xr-xl)
    x = pyclaw.Dimension('x',xl,xr,mx)
    domain = pyclaw.Domain([x])
    state = pyclaw.State(domain,solver.num_eqn)

    state.problem_data['gamma']= gamma
    state.problem_data['gamma1']= gamma1
    state.problem_data['efix'] = True
    #state.problem_data['tfluct_solver'] = True

    xc =state.grid.x.centers
    epsilon=0.8; w = 8.
    #state.q[0,:] = (xc<0.)*1.25*np.exp(-((xc-xl)/w)**2) + (xc>=0)*(1+epsilon*np.sin(5*xc))
    state.q[0,:] = 1.0 + epsilon*np.sin(5*xc)
    velocity = xc*0. #(xc<0.)*2.629369*np.exp(-((xc-xl)/w)**2) 
    state.q[1,:] = velocity * state.q[0,:]
    pressure = 1 + 0.15*np.exp(-((xc-xl)/w)**2)
    state.q[2,:] = pressure/gamma1 + 0.5 * state.q[0,:] * velocity**2

    claw = pyclaw.Controller()
    claw.tfinal = 1600.
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = claw.tfinal
    claw.outdir = outdir
    claw.setplot = setplot

    state.grid.add_gauges([[100.],[500.],[1000.],[1500.],[1800.]])

    return claw

#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Density'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    
    # Figure for q[1]
    plotfigure = plotdata.new_plotfigure(name='Energy', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Density'
    plotaxes.xlimits=(-10.,2000.)

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 2
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    
    return plotdata

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
