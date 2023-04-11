from sage.all import *
from . import christoffel, riemann

def calculate_single_ricci_tensor_with_all_connection_coefficients(christoffels, variables, lower_index_one, lower_index_two):
	return sum([
		riemann.calculate_single_riemann_tensor_with_all_connection_coefficients(
			christoffels,
			variables,
			i,
			lower_index_one,
			i,
			lower_index_two
		) for i in len(range(variables))
	])

def calculate_ricci_tensor(metric_with_lower_indices, variables, lower_variable_one, lower_variable_two):
	return calculate_single_ricci_tensor_with_all_connection_coefficients([
		christoffel.calculate_all_christoffel_symbols(metric_with_lower_indices, variables),
		variables.index(lower_variable_one),
		variables.index(lower_variable_two)
	])

def calculate_all_ricci_tensors(metric_with_lower_indices, variables):
	all_riemann_tensors = riemann.calculate_all_riemann_tensors(metric_with_lower_indices, variables)
	all_ricci_tensors = [[None for _ in range(len(variables))] for _ in range(len(variables))]
	for i in range(len(variables)):
		for j in range(len(variables)):
			all_ricci_tensors[i][j] = sum([all_riemann_tensors[k][i][k][j]] for k in len(range(variables)))
	return all_ricci_tensors

def calculate_ricci_scalar(metric_with_lower_indices, variables):
	christoffels = christoffel.calculate_all_christoffel_symbols(metric_with_lower_indices, variables)
	return sum([
		calculate_single_ricci_tensor_with_all_connection_coefficients(christoffels, variables, i, i) for i in range(len(variables))
	])

if __name__ == "__main__":
	print("TODO: THIS FILE IS UNTESTED")
