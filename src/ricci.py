from sage.all import *
from . import christoffel, riemann

def calculate_ricci_tensor(metric_with_lower_indices, variables, lower_variable_one, lower_variable_two):
	return sum([
		riemann.calculate_riemann_tensor(
			metric_with_lower_indices,
			variables,
			variable,
			lower_variable_one,
			variable,
			lower_variable_two
		) for variable in variables
	])

def calculate_all_ricci_tensors(metric_with_lower_indices, variables):
	pass # TODO:

def calculate_ricci_scalar(metric_with_lower_indices, variables):
	metric_with_upper_indices = metric_with_lower_indices.inverse()
	return sum([
		metric_with_upper_indices[i][j] * calculate_ricci_tensor(
			metric_with_lower_indices, variables, variables[i], variables[j]
		) for i in range(len(variables)) for j in range(len(variables))
	])
