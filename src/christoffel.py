import itertools
from sage.all import *

def calculate_single_christoffel_symbol_with_both_metrics(metric_with_lower_indices, metric_with_upper_indices, variables, upper_index, lower_index_one, lower_index_two):
	lower_variable_one = variables[lower_index_one]
	lower_variable_two = variables[lower_index_two]
	return 1 / 2 * sum([
		metric_with_upper_indices[upper_index][i] * (
			derivative(metric_with_lower_indices[i][lower_index_two], lower_variable_one)
			+ derivative(metric_with_lower_indices[lower_index_one][i], lower_variable_two)
			- derivative(metric_with_lower_indices[lower_index_one][lower_index_two], variables[i])
		) for i in range(len(variables))
	])

def calculate_christoffel_symbol(metric_with_lower_indices, variables, upper_variable, lower_variable_one, lower_variable_two):
	metric_with_upper_indices = metric_with_lower_indices.inverse()
	return calculate_single_christoffel_symbol_with_both_metrics(
		metric_with_lower_indices,
		metric_with_upper_indices,
		variables,
		variables.index(upper_variable),
		variables.index(lower_variable_one),
		variables.index(lower_variable_two)
	)

def calculate_all_christoffel_symbols(metric_with_lower_indices, variables):
	metric_with_upper_indices = metric_with_lower_indices.inverse()
	all_christoffel_symbols = [[[None for _ in variables] for _ in variables] for _ in variables]
	for i, _ in enumerate(variables):
		for lower_variable_one, lower_variable_two in itertools.combinations_with_replacement(variables, 2):
			j = variables.index(lower_variable_one)
			k = variables.index(lower_variable_two)
			all_christoffel_symbols[i][j][k] = calculate_single_christoffel_symbol_with_both_metrics(
				metric_with_lower_indices,
				metric_with_upper_indices,
				variables,
				i,
				j,
				k
			)
			all_christoffel_symbols[i][k][j] = all_christoffel_symbols[i][j][k]
	return all_christoffel_symbols

def print_all_christoffel_symbols(metric_with_lower_indices, variables):
	all_christoffel_symbols = calculate_all_christoffel_symbols(metric_with_lower_indices, variables)
	for i, j, k in itertools.product(range(len(variables)), repeat=3):
		symbol = all_christoffel_symbols[i][j][k]
		if symbol != 0:
			pretty_print(
				LatexExpr(r"\Gamma^{")
				+ latex(variables[i])
				+ LatexExpr(r"}_{")
				+ latex(variables[j])
				+ latex(variables[k])
				+ LatexExpr(r"} = ")
				+ latex(symbol.simplify_full().expand())
			)
			print()

if __name__ == "__main__":
	r, phi = var("r"), var("phi", latex_name=r"\phi")
	metric_with_lower_indices = Matrix([[1, 0], [0, r ** 2]])
	print_all_christoffel_symbols(metric_with_lower_indices, [r, phi])
