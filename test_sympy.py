import sympy as sp

# give me a short tutorial on how to use sympy

# 1. define a symbol
x = sp.Symbol("x")
y = sp.Symbol("y")
z = sp.Symbol("z")

# 2. define an expression
expr = x**2 + 2 * x + 1

# 3. print the expression
print(expr)

# 4. expand the expression
print(sp.expand(expr))

# 5. factor the expression
print(sp.factor(expr))

# 6. simplify the expression
print(sp.simplify(expr))

# 7. substitute a value for a symbol
print(expr.subs(x, 2))

# 8. substitute a value for a symbol and evaluate the expression
print(expr.subs(x, 2).evalf())

# %%
# calculate the future value of a series of cash flows (deposits) in the context of compound interest

# 1. define the symbols
a0 = sp.Symbol("a0")
ry = sp.Symbol("r_y")  # yearly interest rate
m = sp.Symbol("m")
rm = (1 + ry) ** sp.Rational(1, 12) - 1  # monthly interest rate


# 2. define the expression
fv = a0 * (((1 + rm) ** m - 1) / rm)
td = m * a0
ti = fv - td
tipy = ti / m * 12

# calculate the derivative of tipy
dm_tipy = sp.diff(tipy, m)
dr_tipy = sp.diff(tipy, ry)


# %%
