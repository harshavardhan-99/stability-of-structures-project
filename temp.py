from sympy.solvers import linsolve, solve
from sympy.core import Symbol, diff, symbols
from sympy import dsolve, Function, Derivative, Eq, cos, sin, sqrt, tan
from sympy.core.symbol import Dummy
from sympy.printing import sstr

class CB():
	def __init__(self,loading_conditions,end_moments):
		self._loading_conditions = loading_conditions
		self._end_moments = end_moments
	@property
	def loading_conditions(self):
		return self._loading_conditions
	@property
	def end_moments(self):
		return self._end_moments
