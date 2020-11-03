from Beam_column import * 
from sympy import Symbol, symbols
#E, I, P = symbols('E, I, P', positive=True)
c = Beam_column(3, 200000, 10**-3, 78000,{'uniform':8000,'point':[10000,1.5]},[0,0],top='hinged',bottom='hinged')
c.boundary_conditions
c.solve_slope_deflection()
c.deflection()
c.slope()
c.Plot_moment()

c = Beam_column(3, 200000, 10**-3, 15000,{'uniform':8000,'point':[10000,1.5]},[0,0],top='hinged',bottom='hinged')
c.boundary_conditions
{'deflection': [(0, 0), (3, 0)], 'slope': []}
c.solve_slope_deflection()
c.deflection()
0.266666666666667*x**2 - 1.13333333333333*x + 0.0534258456896503*sqrt(6)*sin(3.53553390593274*sqrt(6)*x) + 0.00711111111111111*cos(3.53553390593274*sqrt(6)*x) - 0.00711111111111111
c.slope()
0.533333333333333*x - 0.0251415744421884*sqrt(6)*sin(3.53553390593274*sqrt(6)*x) + 1.13333333333333*cos(3.53553390593274*sqrt(6)*x) - 1.13333333333333
c.Plot_moment()



c = Beam_column(3, 200000, 10**-3, 15000,{'uniform':0,'point':[10000,1.5]},[0,0],top='hinged',bottom='hinged')
c.boundary_conditions
c.solve_slope_deflection()
c.deflection()
c.slope()
c.Plot_moment()
c = Beam_column(3, 200000, 10**-3, 15000,{'uniform':0,'point':[10000,1.5]},[0,0],top='hinged',bottom='hinged')
c.boundary_conditions
{'deflection': [(0, 0), (3, 0)], 'slope': []}
c.solve_slope_deflection()
c.deflection()
-0.333333333333333*x + 0.0157134840263677*sqrt(6)*sin(3.53553390593274*sqrt(6)*x)
c.slope()
0.333333333333333*cos(3.53553390593274*sqrt(6)*x) - 0.333333333333333
c.Plot_moment()
<function _lambdifygenerated at 0x7f3300342950>
def _lambdifygenerated(x):
    return (5000.0*x + 15000*y(x))


from Beam_column import * 
from sympy import Symbol, symbols
c = Beam_column(3, 200000, 10**-3, 15000,{'uniform':10000,'point':[0,1.5]},[0,0],top='fixed',bottom='fixed')
c.boundary_conditions
c.solve_slope_deflection()
c.deflection()
c.slope()
c.Plot_moment()

c.boundary_conditions
{'deflection': [(0, 0), (3, 0)], 'slope': [(0, 0)]}
c.solve_slope_deflection()
c.deflection()
0.333333333333333*x**2 - 1.0*x + 0.0471404520791032*sqrt(6)*sin(3.53553390593274*sqrt(6)*x) - 0.491111111111111*cos(3.53553390593274*sqrt(6)*x) + 0.491111111111111
c.slope()
0.666666666666667*x + 1.73633998491363*sqrt(6)*sin(3.53553390593274*sqrt(6)*x) + 1.0*cos(3.53553390593274*sqrt(6)*x) - 1.0
 c.Plot_moment()
<function _lambdifygenerated at 0x7f3300342b90>
def _lambdifygenerated(x):
    return (-5000*x**2 + 15000*x + 15000*y(x) - 7500)