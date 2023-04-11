from sage.all import *
from . import christoffel

def calculate_single_riemann_tensor_with_all_connection_coefficients(christoffels, variables, upper_index, lower_index_one, lower_index_two, lower_index_three):
	return derivative(christoffels[upper_index][lower_index_three][lower_index_one], variables[lower_index_two]) - \
		derivative(christoffels[upper_index][lower_index_two][lower_index_one], variables[lower_index_three]) + \
		sum(christoffels[upper_index][lower_index_two][i] * christoffels[i][lower_index_three][lower_index_one] for i in range(len(variables))) - \
		sum(christoffels[upper_index][lower_index_three][i] * christoffels[i][lower_index_two][lower_index_one] for i in range(len(variables)))

def calculate_riemann_tensor(metric_with_lower_indices, variables, upper_variable, lower_variable_one, lower_variable_two, lower_variable_three):
	return calculate_single_riemann_tensor_with_all_connection_coefficients(
		christoffel.calculate_all_christoffel_symbols(metric_with_lower_indices, variables),
		variables,
		variables.index(upper_variable),
		variables.index(lower_variable_one),
		variables.index(lower_variable_two),
		variables.index(lower_variable_three)
	)

def calculate_all_riemann_tensors(metric_with_lower_indices, variables):
	christoffels = christoffel.calculate_all_christoffel_symbols(metric_with_lower_indices, variables)
	all_riemann_tensors = [[[[None for _ in range(len(variables))] for _ in range(len(variables))] for _ in range(len(variables))] for _ in range(len(variables))]
	for i in range(len(variables)):
		for j in range(len(variables)):
			for k in range(len(variables)):
				for l in range(len(variables)):
					all_riemann_tensors[i][j][k][l] = calculate_single_riemann_tensor_with_all_connection_coefficients(
						christoffels,
						variables,
						i,
						j,
						k,
						l
					)
	return all_riemann_tensors

if __name__ == "__main__":
	print("TODO: THIS FILE IS UNTESTED")
