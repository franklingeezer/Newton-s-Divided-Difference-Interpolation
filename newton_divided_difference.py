"""
Newton's Divided Difference Interpolation
------------------------------------------
Numerical Methods Presentation - Topic: Interpolation (Unequal Intervals)

This method builds an interpolating polynomial for a set of (x, y) data
points that are NOT necessarily equally spaced. It works by constructing
a divided difference table and using it as the coefficients of Newton's
interpolating polynomial:

    P(x) = f[x0] + f[x0,x1](x - x0) + f[x0,x1,x2](x - x0)(x - x1) + ...

where f[x0, x1, ..., xk] denotes the k-th order divided difference.
"""


def compute_divided_difference_table(x_values, y_values):
    """
    Builds the full divided difference table.

    Returns a 2D list `table` where table[i][j] is the j-th order
    divided difference starting at point i.
    table[i][0] = y_i  (the 0th order "difference" is just the y value)
    """
    n = len(x_values)
    table = [[0.0] * n for _ in range(n)]

    # 0th column is just the y values
    for i in range(n):
        table[i][0] = y_values[i]

    # Fill each subsequent column using the divided difference formula:
    # f[x_i,...,x_{i+j}] = (f[x_{i+1},...,x_{i+j}] - f[x_i,...,x_{i+j-1}])
    #                       / (x_{i+j} - x_i)
    for j in range(1, n):
        for i in range(n - j):
            numerator = table[i + 1][j - 1] - table[i][j - 1]
            denominator = x_values[i + j] - x_values[i]
            table[i][j] = numerator / denominator

    return table


def print_divided_difference_table(x_values, y_values, table):
    """Pretty-prints the divided difference table."""
    n = len(x_values)
    print("\nDivided Difference Table")
    print("-" * 60)
    header = f"{'x':>10} {'f[x]':>12}"
    for order in range(1, n):
        header += f"{'Order ' + str(order):>14}"
    print(header)
    print("-" * 60)

    for i in range(n):
        row = f"{x_values[i]:>10.4f} {table[i][0]:>12.4f}"
        for j in range(1, n - i):
            row += f"{table[i][j]:>14.4f}"
        print(row)
    print("-" * 60)


def newton_interpolate(x_values, y_values, x_target):
    """
    Evaluates the Newton divided difference interpolating polynomial
    at x_target and returns (result, table, coefficients).
    """
    n = len(x_values)
    table = compute_divided_difference_table(x_values, y_values)

    # The coefficients of the Newton polynomial are the top row
    # of the divided difference table: f[x0], f[x0,x1], f[x0,x1,x2], ...
    coefficients = [table[0][j] for j in range(n)]

    result = coefficients[0]
    product_term = 1.0
    for k in range(1, n):
        product_term *= (x_target - x_values[k - 1])
        result += coefficients[k] * product_term

    return result, table, coefficients


def build_polynomial_string(x_values, coefficients):
    """Builds a human-readable string of the Newton polynomial."""
    terms = [f"{coefficients[0]:.4f}"]
    for k in range(1, len(coefficients)):
        factors = "".join(f"(x - {x_values[i]:.4f})" for i in range(k))
        terms.append(f"{coefficients[k]:.4f}{factors}")
    return " + ".join(terms)


def get_user_input():
    """Reads any number of (x, y) points and a target x value from the user."""
    print("=== Newton's Divided Difference Interpolation ===")
    n = int(input("Enter the number of data points: "))

    x_values = []
    y_values = []
    print("Enter the data points (x values do not need equal spacing):")
    for i in range(n):
        x = float(input(f"  x[{i}] = "))
        y = float(input(f"  y[{i}] = "))
        x_values.append(x)
        y_values.append(y)

    x_target = float(input("\nEnter the x value to interpolate at: "))
    return x_values, y_values, x_target


def main():
    x_values, y_values, x_target = get_user_input()

    result, table, coefficients = newton_interpolate(x_values, y_values, x_target)

    print_divided_difference_table(x_values, y_values, table)

    print(f"\nNewton's Interpolating Polynomial:")
    print(f"P(x) = {build_polynomial_string(x_values, coefficients)}")

    print(f"\nP({x_target}) = {result:.6f}")


if __name__ == "__main__":
    main()
