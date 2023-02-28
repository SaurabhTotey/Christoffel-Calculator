from sage.all import *

def calculate_christoffel_symbol(metric_with_lower_indices, variables, lower_variable_one, lower_variable_two, upper_variable):
	lower_index_one = variables.index(lower_variable_one)
	lower_index_two = variables.index(lower_variable_two)
	upper_index = variables.index(upper_variable)
	metric_with_upper_indices = metric_with_lower_indices.inverse()
	return 1 / 2 * sum([
		metric_with_upper_indices[upper_index][i] * (
			derivative(metric_with_lower_indices[i][lower_index_two], lower_variable_one)
			+ derivative(metric_with_lower_indices[lower_index_one][i], lower_variable_two)
			- derivative(metric_with_lower_indices[lower_index_one][lower_index_two], variables[i])
		) for i in range(len(variables))
	])

def print_all_christoffel_symbols(metric_tensor_with_lower_indices, variables):
	# TODO: this could be sped up
	#  the metric with upper indices doesn't need to be recomputed every iteration
	#  swapping the lower variables doesn't change the Christoffel symbol
	for upper_variable in variables:
		for lower_variable_one in variables:
			for lower_variable_two in variables:
				result = calculate_christoffel_symbol(
					metric_tensor_with_lower_indices,
					variables,
					lower_variable_one,
					lower_variable_two,
					upper_variable
				)
				if result != 0:
					pretty_print(
						LatexExpr(r"\Gamma^{")
						+ latex(upper_variable)
						+ LatexExpr(r"}_{")
						+ latex(lower_variable_one)
						+ latex(lower_variable_two)
						+ LatexExpr(r"} = ")
						+ latex(result.simplify_full().expand())
					)
					print()

if __name__ == "__main__":
	r, phi = var("r"), var("phi", latex_name=r"\phi")
	metric_tensor = Matrix([[1, 0], [0, r ** 2]])
	print_all_christoffel_symbols(metric_tensor, [r, phi])
