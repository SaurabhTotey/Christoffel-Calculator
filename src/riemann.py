from sage.all import *
from . import christoffel

def calculate_single_riemann_tensor_with_both_metrics(metric_with_lower_indices, metric_with_upper_indices, variables, upper_index, lower_index_one, lower_index_two, lower_index_three):
	previous_cc_values = {}
	def cc(i, j, k):
		if (i, j, k) in previous_cc_values:
			return previous_cc_values[(i, j, k)]
		previous_cc_values[(i, j, k)] = christoffel.calculate_single_christoffel_symbol_with_both_metrics(metric_with_lower_indices, metric_with_upper_indices, variables, i, j, k)
		return previous_cc_values[(i, j, k)]
	return derivative(cc(upper_index, lower_index_three, lower_index_one), variables[lower_index_two]) - \
		derivative(cc(upper_index, lower_index_two, lower_index_one), variables[lower_index_three]) + \
		sum(cc(upper_index, lower_index_two, i) * cc(i, lower_index_three, lower_index_one) for i in range(len(variables))) - \
		sum(cc(upper_index, lower_index_three, i) * cc(i, lower_index_two, lower_index_one) for i in range(len(variables)))

# TODO: this only saves the computation of the inverse of the metric, but it's really the connection coefficients that should be being saved
#  then the calculate_all_riemann_tensors function wouldn't need a separate inner function
def calculate_riemann_tensor(metric_with_lower_indices, variables, upper_variable, lower_variable_one, lower_variable_two, lower_variable_three):
	metric_with_upper_indices = metric_with_lower_indices.inverse()
	return calculate_single_riemann_tensor_with_both_metrics(
		metric_with_lower_indices,
		metric_with_upper_indices,
		variables,
		variables.index(upper_variable),
		variables.index(lower_variable_one),
		variables.index(lower_variable_two),
		variables.index(lower_variable_three)
	)

def calculate_all_riemann_tensors(metric_with_lower_indices, variables):
	christoffels = christoffel.calculate_all_christoffel_symbols(metric_with_lower_indices, variables)
	def calculate_single_riemann_tensor(upper_index, lower_index_one, lower_index_two, lower_index_three):
		lower_variable_two = variables[lower_index_two]
		lower_variable_three = variables[lower_index_three]
		return derivative(christoffels[upper_index][lower_index_three][lower_index_one], lower_variable_two) - \
			derivative(christoffels[upper_index][lower_index_two][lower_index_one], lower_variable_three) + \
			sum(christoffels[upper_index][lower_index_two][i] * christoffels[i][lower_index_three][lower_index_one] for i in range(len(variables))) - \
			sum(christoffels[upper_index][lower_index_three][i] * christoffels[i][lower_index_two][lower_index_one] for i in range(len(variables)))
	all_riemann_tensors = [[[[None for _ in range(len(variables))] for _ in range(len(variables))] for _ in range(len(variables))] for _ in range(len(variables))]
	for i in range(len(variables)):
		for j in range(len(variables)):
			for k in range(len(variables)):
				for l in range(len(variables)):
					all_riemann_tensors[i][j][k][l] = calculate_single_riemann_tensor(i, j, k, l)
	return all_riemann_tensors

if __name__ == "__main__":
	print("TODO: THIS FILE IS UNTESTED")
