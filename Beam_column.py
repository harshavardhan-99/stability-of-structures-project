from sympy.solvers import linsolve, solve
from sympy.core import Symbol, diff, symbols
from sympy import dsolve, Function, Derivative, Eq, cos, sin, sqrt, tan
from sympy.core.symbol import Dummy
from sympy.printing import sstr
from sympy.plotting import plot, PlotGrid
from sympy import sin, cos, symbols, lambdify
import numpy as np 

class Beam_column(object):
    """
    A Beam_column is a structural member designed to undertake both axial
    compressive and Flexural loads. A column is characterized by its
    cross-sectional profile(second moment of area), its length and
    its material.

    Examples
    ========
    There is a solid round bar 3 m long with second-moment I is used as a
    column with both the ends pinned. Young's modulus of the Column is E.
    The buckling load applied is 78KN

    >>> from sympy.physics.continuum_mechanics.column import Column
    >>> from sympy import Symbol, symbols
    >>> E, I, P = symbols('E, I, P', positive=True)
    >>> c = Beam_column(3, E, I, 15000,{'uniform':8000,'point':[10000,5]},[0,0],top='hinged',bottom='hinged') 
    top="pinned", bottom="pinned")
    >>> c.end_conditions
    {'bottom': 'pinned', 'top': 'pinned'}

    >>> c.loading_condition
    {'uniform':w = 10000,'point':[5000,1.5]}
    >>> c.end_moments('Ma'=M1,'Mb' =M2)
    >>> c.boundary_conditions
    {'deflection': [(0, 0), (3, 0)], 'slope': [(0, 0)]}
    >>> c.moment()
    78000*y(x)
    >>> c.solve_slope_deflection()
    >>> c.deflection()
    C1*sin(20*sqrt(195)*x/(sqrt(E)*sqrt(I)))
    >>> c.slope()
    20*sqrt(195)*C1*cos(20*sqrt(195)*x/(sqrt(E)*sqrt(I)))/(sqrt(E)*sqrt(I))
    >>> c.critical_load()
    pi**2*E*I/9
    """
    def __init__(self, height, elastic_modulus, second_moment, load,loading_conditions={},end_moments=[],top="pinned", bottom="pinned",eccentricity=0, bc_slope=None, bc_deflection=None):
        """
        Parameters
        ==========

        height: Sympifyable
            A symbol or a value representing column's height

        elastic_modulus: Sympifyable
            A symbol or a value representing the Column's modulus of
            elasticity. It is a measure of the stiffness of the Column
            material.

        second_moment: Sympifyable
            A symbol or a value representing Column's second-moment of area
            It is a geometrical property of an area which reflects how its
            points are distributed with respect to its neutral axis.

        load_p: Sympifyable
            A symbol or a value representing the load applied on the Column.
        
        load_f: Sympifyable
            A symbol or a value representing the load applicaton and magnitude of the load.

        eccentricity: Sympifyable (default=None)
            A symbol or a value representing the eccentricity of the load
            applied. Eccentricity is the distance of the point of application
            of load from the neutral axis.

        top: string (default="pinned")
            A string representing the top-end condition of the column.
            It can be: pinned
                       fixed
                       free
        bottom: string (default="pinned")
            A string representing the bottom-end condition of the column.
            It can be: pinned
                       fixed
        bc_slope: list of tuples
            A list of tuples representing the boundary conditions of slope.
            The tuple takes two elements `location` and `value`.
        bc_deflection: list of tuples
            A list of tuples representing the boundary conditions of deflection
            The tuple consists of two elements `location` and `value`.
        """
        self._height = height
        self._elastic_modulus = elastic_modulus
        self._second_moment = second_moment
        self._load = load
        self._eccentricity = eccentricity
        self._moment = 0
        self._end_conditions = {'top':top, 'bottom': bottom}
        self._boundary_conditions = {'deflection': [], 'slope': []}
        if bc_deflection:
            self._boundary_conditions['deflection'] = bc_deflection
        if bc_slope:
            self._boundary_conditions['slope'] = bc_slope
        #self._loading_conditions = {'uniform':w,'point':[point_load-value,distance]}
        #self._end_moments ={'Ma':Ma,'Mb':Mb}
        self._loading_conditions = loading_conditions
        self._end_moments = end_moments
        self._variable = Symbol('x')
        #self._y = Function('y')
        self._deflection = None
        self._slope = None
        self._critical_load = None
        self._apply_load_conditions()

    def __str__(self):
        str_sol = 'Column({}, {}, {})'.format(sstr(self._height), sstr(self._elastic_modulus), sstr(self._second_moment))
        return str_sol

    @property
    def height(self):
        """Height of the column"""
        return self._height

    @property
    def elastic_modulus(self):
        """Elastic modulus of the column"""
        return self._elastic_modulus

    @property
    def second_moment(self):
        """Second moment of the column"""
        return self._second_moment

    @property
    def load(self):
        """Load applied on the column"""
        return self._load

    @property
    def eccentricity(self):
        """Eccentricity of the load applied on the column"""
        return self._eccentricity

    @property
    def end_conditions(self):
        """End-conditions in the form of a dictionary"""
        return self._end_conditions

    @property
    def boundary_conditions(self):
        """Boundary conditions in the form of a dictionary"""
        return self._boundary_conditions

    @property
    def loading_conditions(self):
        """loading condition  """
        return self._loading_conditions

    @property
    def end_moments(self):
        """loading condition  """
        return self._end_moments

    def _apply_load_conditions(self):
        y = Function('y')
        x = self._variable
        P = Symbol('P', positive=True)
        Ma = Symbol('Ma')
        Mb = Symbol('Mb')
        Q = Symbol('Q')
        w = Symbol('w')
        a = self._loading_conditions['point'][1]

        self._moment += P*y(x)
        if self.eccentricity:
            self._moment += P*eccentricity

        # Initial boundary conditions, considering slope and deflection
        # the bottom always zero
        self._boundary_conditions['deflection'].append((0, 0))
        self._boundary_conditions['deflection'].append((self._height, 0))

        if self._end_conditions['top'] == "pinned" and self._end_conditions['bottom'] == "pinned":
            self._boundary_conditions['deflection'].append((self._height, 0))

        elif self._loading_conditions['uniform'] == 0 and self._loading_conditions['point'][0] == 0 and self._end_conditions['top'] == "hinged" and self._end_conditions['bottom'] == "hinged":
            # `Ma` is the end reaction moment
            # `Mb` is the end reaction moment
            #self._boundary_conditions['deflection'].append((self._height, 0))
            # moment  = P*y - ((Ma + Mb)*(x/self.height) + Ma)
            self._moment -= ((Ma + Mb)*(x/self.height) + Ma) 

        elif self._loading_conditions['uniform'] == 0 and self._loading_conditions['point'][0] != 0 and self._end_conditions['top'] == "hinged" and self._end_conditions['bottom'] == "hinged":
            # `Ma` is the end reaction moment
            # `Mb` is the end reaction moment
            #self._boundary_conditions['deflection'].append((self._height, 0))
            # `Q` is the horizontal force at a length form fixed end
            # moment = P*y + Q*(self.height - a)/self.height - ((Ma + Mb)*(x/self.height) + Ma)
            self._moment += Q/2*x* - ((Ma + Mb)*(x/self.height) + Ma)

        elif self._loading_conditions['uniform'] == 0 and self._loading_conditions['point'][0] != 0 and self._end_conditions['top'] == "fixed" and self._end_conditions['bottom'] == "fixed":
            # `Ma` is the end reaction moment
            # `Mb` is the end reaction moment
            self._boundary_conditions['slope'].append((0, 0))
            #self._boundary_conditions['deflection'].append((self._height, 0))
            # `0` is the deflection at the free end
            # moment = P*y + (w*x**2/2) + (w*self.height*x)/2 - ((Ma + Mb)*(x/self.height) + Ma)
            self._moment += (Q*x)/2 -(Q*self.height)/8

        elif self._loading_conditions['uniform'] != 0 and self._loading_conditions['point'][0] == 0 and self._end_conditions['top'] == "hinged" and self._end_conditions['bottom'] == "hinged":
            # `Ma` is the end reaction moment
            # `Mb` is the end reaction moment
            #self._boundary_conditions['deflection'].append((self._height, 0))
            # `0` is the deflection at the free end
            # moment = P*y + (w*x**2/2) + (w*self.height*x)/2 - ((Ma + Mb)*(x/self.height) + Ma)
            self._moment += -(w*x**2/2) + (w*self.height*x)/2 - ((Ma + Mb)*(x/self.height) + Ma)

        elif self._loading_conditions['uniform'] != 0 and self._loading_conditions['point'][0] != 0 and self._end_conditions['top'] == "hinged" and self._end_conditions['bottom'] == "hinged":
            # `Ma` is the end reaction moment
            # `Mb` is the end reaction moment
            #self._boundary_conditions['deflection'].append((self._height, 0))
            # `0` is the deflection at the free end
            #w = unform load 
            #Q = Point load 
            # moment = P*y - combing above above 2 differential equartions to obtain this
            #if(x<a):
            print('Harsha')
            print(Ma , Mb)
            self._moment += -(w*x**2/2) + (w*self.height*x)/2 + Q/2*x - ((Ma + Mb)*(x/self.height) + Ma)
            #else:
            #    self._moment += 

        elif self._loading_conditions['uniform'] != 0 and self._loading_conditions['point'][0] == 0 and self._end_conditions['top'] == "fixed" and self._end_conditions['bottom'] == "fixed":
            # `Ma` is the end reaction moment
            # `Mb` is the end reaction moment
            #self._boundary_conditions['deflection'].append((self._height, 0))
            self._boundary_conditions['slope'].append((0, 0))
            self._boundary_conditions['slope'].append((self._height, 0))

            # `0` is the deflection at the free end
            # moment = P*y + (w*x**2/2) + (w*self.height*x)/2 - ((Ma + Mb)*(x/self.height) + Ma)
            self._moment += -(w*x**2/2) + (w*self.height*x)/2 -(w*self.height**2)/12

        else:
            raise ValueError("{} {} end-condition is not supported".format(sstr(self._end_conditions['top']), sstr(self._end_conditions['bottom'])))

    def solve_slope_deflection(self):
        """
        Solves the differnetial equation of buckling to determine the
        deflection and slope equations.

        Examples
        ========

        A column of timber section is 8m long both ends being fixed.
        Young's modulus of the timber is E and second-moment I.
        The applied load on the column is 360KN.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy import Symbol, symbols
        >>> E, I, P = symbols('E, I, P', positive=True)
        >>> c = Column(8, E, I, 360000, top="fixed", bottom="fixed")
        >>> c.end_conditions
        {'bottom': 'fixed', 'top': 'fixed'}
        >>> c.boundary_conditions
        {'deflection': [(0, 0), (8, 0)], 'slope': [(0, 0)]}
        >>> c.moment()
        -M + 360000*y(x)
        >>> c.solve_slope_deflection()
        >>> c.deflection()
        -M*cos(600*x/(sqrt(E)*sqrt(I)))/360000 + M/360000
        >>> c.slope()
        M*sin(600*x/(sqrt(E)*sqrt(I)))/(600*sqrt(E)*sqrt(I))
        """
        y = Function('y')
        x = self._variable
        P = Symbol('P', positive=True)

        C1, C2 = symbols('C1, C2')
        E = self._elastic_modulus
        I = self._second_moment

        # differnetial equation of Beam_Column buckling
        eq = E*I*y(x).diff(x, 2) + self._moment

        defl_sol = dsolve(Eq(eq, 0), y(x)).args[1]
        slope_sol = defl_sol.diff(x)
        #print(defl_sol)
        #print(slope_sol)
        constants = list(linsolve([defl_sol.subs(x, 0), defl_sol.subs(x, self._height)], C1, C2).args[0])
        #print(constants)
        self._deflection = defl_sol.subs({C1: constants[0], C2:constants[1]})
        self._slope = slope_sol.subs({C1: constants[0], C2:constants[1]})
        #self._moment  = E*I*y(x).diff(x, 2)
        # if deflection is zero, no buckling occurs, which is not the case,
        # so trying to solve for the constants differently
        print(self._deflection)
        print(self._slope)
        if self._deflection == 0:
            self._deflection = defl_sol
            self._slope = slope_sol

            defl_eqs = []
            # taking last two bounndary conditions which are actually
            # the initial boundary conditions.
            for point, value in self._boundary_conditions['deflection'][-2:]:
                defl_eqs.append(self._deflection.subs(x, point) - value)

            # solve for C1, C2 along with P
            solns = solve(defl_eqs, (P, C1, C2), dict=True)
            for sol in solns:
                if self._deflection.subs(sol) == 0:
                    # removing trivial solutions
                    solns.remove(sol)

            # checking if the constants are solved, and subtituting them in
            # the deflection and slope equation
            #print(solns[0].keys())
            if C1 in solns[0].keys():
                self._deflection = self._deflection.subs(C1, solns[0][C1])
                self._slope = self._slope.subs(C1, solns[0][C1])
            if C2 in solns[0].keys():
                self._deflection = self._deflection.subs(C2, solns[0][C2])
                self._slope = self._slope.subs(C2, solns[0][C2])
            if P in solns[0].keys():
                self._critical_load = solns[0][P]




    def critical_load(self):
        """
        Detrmines the critical load (for single bow buckling condition) of
        the given column under the given conditions.

        Examples
        ========
        from Beam_column import * 
        from sympy import Symbol, symbols
        E, I, P = symbols('E, I, P', positive=True)
        c = Beam_column(3, E, I, 15000,{'uniform':8000,'point':[10000,5]},[0,0],top='hinged',bottom='hinged')
        c.boundary_conditions
        c.solve_slope_deflection()
        c.deflection()
        c.slope()





        A solid round bar 10 long is used as a column. The young's modulus
        is E and the second moment is I. One end of the column is fixed
        and other end is free. A load of 15KN is applied on it.

        >>> from sympy.physics.continuum_mechanics.column import Column
        >>> from sympy import Symbol, symbols
        >>> E, I, P = symbols('E, I, P', positive=True)
        >>> c = Column(10, E, I, 15000, top="free", bottom="fixed")
        >>> c.end_conditions
        {'bottom': 'fixed', 'top': 'free'}
        >>> c.boundary_conditions
        {'deflection': [(0, 0), (10, d)], 'slope': [(0, 0)]}
        >>> c.moment()
        -15000*d + 15000*y(x)
        >>> c.solve_slope_deflection()
        >>> c.deflection()
        -d*cos(50*sqrt(6)*x/(sqrt(E)*sqrt(I))) + d
        >>> c.slope()
        50*sqrt(6)*d*sin(50*sqrt(6)*x/(sqrt(E)*sqrt(I)))/(sqrt(E)*sqrt(I))
        >>> c.critical_load()
        pi**2*E*I/400
        """
        y = Function('y')
        x = self._variable
        P = Symbol('P', positive=True)
        Ma = Symbol('Ma')
        Mb = Symbol('Mb')
        Q = Symbol('Q')
        w = Symbol('w')

        if self._critical_load is None:
            defl_eqs = []
            # taking last two bounndary conditions which are actually
            # the initial boundary conditions.
            #print(self._boundary_conditions['deflection'])
            for point, value in self._boundary_conditions['deflection'][-2:]:
                defl_eqs.append(self._deflection.subs(x, point) - value)

            # C1, C2 already solved, solve for P
            self._critical_load = solve(defl_eqs, P, dict=True)[0][P]

        return self._critical_load


    def moment(self):
        """Returns the moment equation in terms of any arbitrary point ``x``
        on the column and the deflection at that point.
        """
        P = Symbol('P', positive=True)
        Ma = Symbol('Ma')
        Mb = Symbol('Mb')
        Q = Symbol('Q')
        w = Symbol('w')
        #eq = E*I*y(x).diff(x, 2)
        return self._moment.subs([(P,self._load),
                                (Ma,self.end_moments[0]),
                                (Mb,self.end_moments[1]),
                                (Q,self.loading_conditions['point'][0]),
                                (w,self.loading_conditions['uniform'])])
        #return self._moment.subs(P, self._load)


    def slope(self):
        """Returns the slope equation in terms of any arbitrary point ``x``
        on the column.
        """
        P = Symbol('P', positive=True)
        Ma = Symbol('Ma')
        Mb = Symbol('Mb')
        Q = Symbol('Q')
        w = Symbol('w')
        #return self._slope.subs(P, self._load)
        #if(self.loading_conditions)
        return self._slope.subs([(P,self._load),
                                (Ma,self.end_moments[0]),
                                (Mb,self.end_moments[1]),
                                (Q,self.loading_conditions['point'][0]),
                                (w,self.loading_conditions['uniform'])])


    def deflection(self):
        """Returns the deflection equation in terms of any arbitrary point ``x``
        on the column.
        """
        P = Symbol('P', positive=True)
        Ma = Symbol('Ma')
        Mb = Symbol('Mb')
        Q = Symbol('Q')
        w = Symbol('w')
        #return self._deflection.subs(P, self._load)
        return self._deflection.subs([(P,self._load),
                                (Ma,self.end_moments[0]),
                                (Mb,self.end_moments[1]),
                                (Q,self.loading_conditions['point'][0]),
                                (w,self.loading_conditions['uniform'])])

    def Plot_moment(self,subs=None):
        import inspect
        P = Symbol('P', positive=True)
        Ma = Symbol('Ma')
        Mb = Symbol('Mb')
        Q = Symbol('Q')
        w = Symbol('w')

        #moment = c.moment()

        height = self._height
        variable = self._variable
        if subs is None:
            subs = {}
        for sym in self.deflection().atoms(Symbol):
            if sym == self._variable:
                continue
            if sym not in subs:
                raise ValueError('Value of %s was not passed.' %sym)
        if self.height in subs:
            height = subs[self._height]
        else:
            height = self._height
        f = lambdify(self._variable, self.deflection(), 'numpy')
        print(f)
        a = np.arange(0,height+0.1,0.1)
        print(inspect.getsource(f))
        diff = (f(a))
        for i in diff:
            if i<=0:
                print(i)
        #plot(self.moment().subs(subs), (self._variable, 0, height))