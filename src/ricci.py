import itertools
from sage.all import *
import christoffel
import riemann

def calculate_single_ricci_tensor_value_with_all_connection_coefficients(christoffels, variables, lower_index_one, lower_index_two):
	return sum([
		riemann.calculate_single_riemann_tensor_value_with_all_connection_coefficients(
			christoffels,
			variables,
			i,
			lower_index_one,
			i,
			lower_index_two
		) for i in len(range(variables))
	])

def calculate_ricci_tensor_value(metric_with_lower_indices, variables, lower_variable_one, lower_variable_two):
	return calculate_single_ricci_tensor_value_with_all_connection_coefficients([
		christoffel.calculate_all_christoffel_symbols(metric_with_lower_indices, variables),
		variables.index(lower_variable_one),
		variables.index(lower_variable_two)
	])

def calculate_all_ricci_tensor_values(metric_with_lower_indices, variables):
	riemann_tensor = riemann.calculate_all_riemann_tensor_values(metric_with_lower_indices, variables)
	ricci_tensor = [[None for _ in range(len(variables))] for _ in range(len(variables))]
	for lower_variable_one, lower_variable_two in itertools.combinations_with_replacement(variables, 2):
		i = variables.index(lower_variable_one)
		j = variables.index(lower_variable_two)
		ricci_tensor[i][j] = sum([riemann_tensor[k][i][k][j] for k in range(len(variables))])
		ricci_tensor[j][i] = ricci_tensor[i][j]
	return ricci_tensor

def calculate_ricci_scalar(metric_with_lower_indices, variables):
	metric_with_upper_indices = metric_with_lower_indices.inverse()
	ricci_tensor = calculate_all_ricci_tensor_values(metric_with_lower_indices, variables)
	return sum([metric_with_upper_indices[i][j] * ricci_tensor[i][j] for i in range(len(variables)) for j in range(len(variables))])

if __name__ == "__main__":
	r, phi, theta = var("r"), var("phi", latex_name=r"\phi"), var("theta", latex_name=r"\theta")
	metric = Matrix([[r ** 2, 0], [0, r ** 2 * sin(theta) ** 2]])
	print(calculate_all_ricci_tensor_values(metric, [theta, phi]))
