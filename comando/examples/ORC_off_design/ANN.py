"""Code to create functions corresponding to artificial neural networks."""
# Copyright © 2020 Institute of Energy and Climate Research
# Energy Systems Engineering (IEK-10)
# Forschungszentrum Jülich GmbH
# Tel.: +49 2461 61-96307
# http://www.fz-juelich.de/iek/iek-10/EN
# 52425 Jülich, Germany
#
# AUTHORS: Marco Langiu
import numpy as np

import comando
from comando import Variable, Parameter, Problem

BOUNDS = -10, 10


class Neuron:
    """A simple model for the transfer function of a neuron.

    An input vector is mapped to one output via the expression:
        activation_function(weights.dot(inputs) + bias)
    """

    def __init__(self, weights, bias, activation_function):
        self.n_in = len(weights)
        self.weights = np.array(weights)
        self.bias = bias
        self.activation_function = activation_function

    def __call__(self, inputs):
        """Get the neuron's activation with the given inputs."""
        return self.activation_function(self.weights.dot(inputs) + self.bias)


class FFNN:
    """A simple model for a Feed Forward Neural Network (FFNN).

    A FFNN is modeled via a list of layers, each of which contains a list of
    Neurons. Input weights, bias and activation function for each neuron is
    thus given within two levels of lists.

    For each layer l an input vector is mapped to an output vector:
        for i, neuron in enumerate(layer):
            activation_function(weights[l][i].dot(inputs[l]) + bias[l][i])
    The input to the first layer is thus mapped to the output of the last one.

    Arguments
    ---------
    weight_arrays : list of list of arrays with input weights
    bias_lists : list of list of bias values
    af_lists : list of list of activation functions
    in_bounds : bounds for input variables assumed during training
    out_bounds : bounds for output variables assumed during training
    """

    def __init__(self, weight_arrays, bias_lists, af_lists,
                 in_bounds=None, out_bounds=None):
        assert len(weight_arrays) == len(bias_lists) == len(af_lists), \
            "the length of `weight_arrays`, `bias_lists` and `af_lists` must" \
            " be the same!"
        self.layers = []
        params = zip(weight_arrays, bias_lists, af_lists)
        for weights_list, bias_list, af_list in params:
            layer = []
            neuron_params = zip(weights_list, bias_list, af_list)
            for weights, bias, activation_function in neuron_params:
                layer.append(Neuron(weights, bias, activation_function))
            self.layers.append(layer)
        self.n_in = len(self.layers[0])
        self.n_out = len(self.layers[-1])
        if in_bounds is None:
            in_bounds = [(-1, 1) for i in range(self.n_in)]
        if out_bounds is None:
            out_bounds = [(-1, 1) for i in range(self.n_out)]
        self.in_lbs = np.array([b[0] for b in in_bounds])
        self.in_delta = np.array([b[1] - b[0] for b in in_bounds])
        self.out_lbs = np.array([b[0] for b in out_bounds])
        self.out_delta = np.array([b[1] - b[0] for b in out_bounds])

    def input_transform(self, data):
        """Transform from input data to the interval [-1, 1]."""
        # return (np.array(data) - self.in_lbs) / self.in_delta * 2 - 1
        return np.array(data) * 2 / self.in_delta \
            - self.in_lbs * 2 / self.in_delta - 1

    def output_transform(self, data):
        """Transform from the interval [-1, 1] to output data."""
        # return (np.array(data) + 1) * self.out_delta / 2 + self.out_lbs
        return np.array(data) * self.out_delta / 2 \
            + self.out_delta / 2 + self.out_lbs

    def __call__(self, *inputs):
        """Compute the outputs of the FFNN with given inputs."""
        inputs = self.input_transform(inputs)
        for layer in self.layers:
            inputs = np.array([neuron(inputs) for neuron in layer])
        return self.output_transform(inputs)

    def dump(self, filename):
        """Dump the FFNN's current parameter values to json."""
        import json

        weights = []
        biases = []
        for layer in self.layers:
            layer_weights = []
            layer_biases = []
            for neuron in layer:
                layer_weights.append(neuron.weights.tolist())
                layer_biases.append(neuron.bias)
            weights.append(layer_weights)
            biases.append(layer_biases)

        data = {
            'in_lbs': self.in_lbs.tolist(),
            'in_delta': self.in_delta.tolist(),
            'weights': [l_weights for l_weights in weights],
            'biases': [l_biases for l_biases in biases],
            'out_lbs': self.out_lbs.tolist(),
            'out_delta': self.out_delta.tolist()
        }
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)


class ScalarFFNN(FFNN):
    """An ANN with only a single output."""

    def __call__(self, *inputs):
        """Compute the output of the ScalarFFNN with given inputs."""
        return FFNN.__call__(self, *inputs)[0]


def ann_from_json(filename: str):
    """Create a function from a json encoded ANN."""
    import json
    with open(filename, 'r') as f:
        data = json.load(f)
    in_lbs = np.asarray(data['in_lbs'])
    in_delta = np.asarray(data['in_delta'])
    weights = [np.asarray(l_weights) for l_weights in data['weights']]
    biases = [np.asarray(l_biases) for l_biases in data['biases']]
    out_lbs = np.asarray(data['out_lbs'])
    out_delta = np.asarray(data['out_delta'])
    afuncs_num = [*[np.tanh] * (len(weights) - 1), lambda x: np.atleast_1d(x)]
    afuncs_sym = [*[comando.tanh] * (len(weights) - 1), lambda x: x]

    def ann_func(*args):
        # NOTE: The seemingly strange use of transposition ensures that we can
        #       also use this function with multidimensional arrays, e.g., for
        #       convenient plotting!

        # Scale from input range to [-1, 1]
        input = ((np.array(args).T - in_lbs) / in_delta * 2 - 1).T

        # Propagate scaled input through layers
        if any(isinstance(arg, comando.Expr) for arg in args):
            for l_weights, l_biases, afunc in zip(weights, biases, afuncs_sym):
                input = np.array([afunc(sum(i * w
                                            for i, w in zip(input.T, weight))
                                        + bias)
                                  for weight, bias
                                  in zip(l_weights, l_biases)])
        else:
            for l_weights, l_biases, afunc in zip(weights, biases, afuncs_num):
                input = afunc([(input.T @ weight).T + bias
                               for weight, bias
                               in zip(l_weights, l_biases)])
        # Revert scaled output to original range
        return (((input + 1) / 2) * out_delta + out_lbs)[0]

    ann_func.__name__ = str(filename).split('.')[0]

    return ann_func


MODELS = {'h_liq': ('liq', 'ps', 'h'),
          'h__sat_liq': ('sat', 'p', 'h_sat_liq'),
          'h__sat_vap': ('sat', 'p', 'h_sat_vap'),
          'h_vap': ('vap', 'ps', 'h'),
          'T_vap': ('vap', 'ph', 'T'),
          's_liq': ('liq', 'ph', 's'),
          's__sat_liq': ('sat', 'p', 's_sat_liq'),
          's__sat_vap': ('sat', 'p', 's_sat_vap'),
          's_vap': ('vap', 'ph', 's'),
          'T_liq': ('liq', 'ph', 'T'),
          'T__sat': ('sat', 'p', 'T_sat'),
          'Pr_liq': ('liq', 'ph', 'Pr'),
          'Pr_vap': ('vap', 'ph', 'Pr'),
          'rho_vap': ('vap', 'ps', 'rho'),
          'v_liq': ('liq', 'ph', 'rho_inv'),
          'v_vap': ('vap', 'ph', 'rho_inv')}


def make_ann(quantity=None, **spec):
    """Create an ANN from the 'Network Databank' for the ORC working fluid.

    Arguments
    ---------
    state: either 'liq', 'sat', or 'vap'
    model:
        - if state == 'liq': either 'ph' or 'ps'
        - if state == sat: 'p'
        - if state == vap: either 'ph', 'ps' or 'pT'
    quantity:
        - if state == 'liq': either 'eta_inv', 'lambda', 'Pr', 'rho',
                                    'rho_inv', T, and
                                    - if model == 'ph': s
                                    - if model == 'ps': h
        - if state == sat: either 'h_sat_liq', 'h_sat_vap', 'rho_sat_liq',
                                  'rho_sat_vap', 's_sat_liq', 's_sat_vap' or
                                  'T_sat'
        - if state == vap: either 'eta_inv', 'lambda', 'Pr', 'rho',
                                    'rho_inv', and
                                    - if model == 'ph': s and T
                                    - if model == 'ps': h and T
                                    - if model == 'pT': h and s
    """
    import os
    import csv
    import comando
    tanh = comando.tanh

    if spec:
        try:
            state, model, qty = spec['state'], spec['model'], spec['qty']
            quantity = qty
        except KeyError as e:
            print('When specifying the desired ANN via keyword arguments you '
                  "provide 'state', 'model' and 'qty'!")
            raise e
    elif quantity in MODELS:
        state, model, qty = MODELS[quantity]
    else:
        raise NotImplementedError(f'No model for quantity {quantity} '
                                  'is available!')

    base_name = f'{os.path.dirname(os.path.realpath(__file__))}/' \
        f'Network Databank/{state}/{model}/{qty}/{state}_{model}_{qty}'

    def vals(key, lines):
        with open(f'{base_name}_{key}.csv', 'r') as f:
            reader = csv.reader(f)
            result = [[float(s) for s in row] for i, row in enumerate(reader)
                      if i in lines]
        return result

    weight_arrays = []
    weight_arrays.append(vals('IW', range(1, 7)))   # input layer
    weight_arrays.append(vals('LW', range(4, 10)))  # first hidden layer
    weight_arrays.append(vals('LW', [14]))          # second hidden layer
    af_lists = [(tanh, ) * 6, (tanh, ) * 6, [lambda x: x]]

    in_bounds = vals('bounds', [0, 1])
    out_bounds = vals('bounds', [2, 3])
    print('Input bounds:')
    for v, lb, ub in zip(model, *in_bounds):
        print(f'  {v}: [{lb}, {ub}]')
    print('Output bounds:')
    print(f'  {quantity}: [{out_bounds[0][0]}, {out_bounds[1][0]}]')

    return ScalarFFNN(weight_arrays, vals('BW', [1, 3, 5]), af_lists,
                      [*zip(*in_bounds)], [*zip(*out_bounds)])


def CachedVar(name, cache=None):
    """If a Variable with that name was created before, return that."""
    if cache is None:
        return Variable(name, bounds=BOUNDS)
    if name in cache:
        return cache[name]
    else:
        var = Variable(name, bounds=BOUNDS)
        cache[name] = var
        return var


def make_fit_ann(neurons_per_layer, in_bounds=None, out_bounds=None):
    """Make an ANN for fitting, i.e., with variables for weights and biases.

    # Naming
    (l0)  l1  l2  ... ln
    --------------------
    i1    n11 n21     o1
    i2    n12 n22     o2
    i3    n13
    """
    try:
        cache = make_fit_ann.cache
    except AttributeError:
        cache = make_fit_ann.cache = {}

    weight_arrays = []
    bias_lists = []
    n_in = neurons_per_layer[0]
    # NOTE: By taking the number of hidden layers as part of the variable name,
    #       we can re-use values of previous optimizations with the same number
    #       of hidden layers and thus hopefully obtain a good initial guess
    #       when broadening the network.
    n_h = len(neurons_per_layer) - 2
    tanh = comando.tanh
    hidden_layer_neurons = neurons_per_layer[1:-1]
    af_lists = [*((tanh, ) * n for n in hidden_layer_neurons), [lambda x: x]]
    for l, n in enumerate(neurons_per_layer[1:], 1):  # For each layer
        weight_array = np.empty((n, n_in), 'object')
        bias_list = np.empty(n, 'object')
        for ni in range(n):  # For each neron per layer
            # For each input of the previous layer
            weight_array[ni] = [CachedVar(f'w{n_h}_{l}_{ni + 1}_{i}', cache)
                                for i in range(1, n_in + 1)]
            bias_list[ni] = CachedVar(f'b{n_h}_{l}_{ni}', cache)
        weight_arrays.append(weight_array)
        bias_lists.append(bias_list)
        n_in = n

    fit_ann = \
        ScalarFFNN(weight_arrays, bias_lists, af_lists, in_bounds, out_bounds)

    def dump(filename):
        weights = []
        biases = []
        for layer in fit_ann.layers:
            layer_weights = []
            layer_biases = []
            for neuron in layer:
                layer_weights.append([weight.value
                                      for weight in neuron.weights])
                layer_biases.append(neuron.bias.value)
            weights.append(layer_weights)
            biases.append(layer_biases)
        tanh = comando.tanh
        af_lists = [*((tanh, ) * n for n in hidden_layer_neurons),
                    [lambda x: x]]
        ScalarFFNN(weights, biases, af_lists, in_bounds, out_bounds). \
            dump(filename)

    fit_ann.dump = dump
    return fit_ann


def make_data(expr, points=None):
    """Make data for a regression.

    Arguments
    ---------
    expr : Expression
        The mathematical expression for which data is to be generated
    points : Mapping
        A mapping from the variables in `expr` to an iterable of points.
        A subset of variables can also be specified, for all other variables
        a 20-point linspace between lower and upper bound is used.

    Returns
    -------
    data : dict
        A dictionary with values for `expr` at the given `points`
    """
    # We lambdify the expression for faster evaluation, to do this we need to
    # replace our custom Variable symbols with SymPy Dummy symbols
    from itertools import product
    from comando.utility import lambdify, get_vars

    # from comando.utility import get_vars, get_pars
    # pars = {p: p.value for p in get_pars(expr)}
    func = lambdify(expr)

    if points is None:
        points = dict()
    vars = get_vars(expr)
    for v in vars:
        if v not in points:
            from numpy import linspace
            points[v] = linspace(*v.bounds)

    # We then evaluate the function on a grid and store the results
    return {val: func(*val) for val in product(*(points[v] for v in vars))}


def make_max_err_fit(vars, data_or_func, hidden_layer_neurons, tol=1e-3,
                     T_max=60, points_per_variable=100, initial_points=None,
                     iter_time=5, enforced_error_bound=float('inf'),
                     solver='baron', filename=None, **solver_options):
    """Fit ANN, iteratively adding points that give largest error."""
    from itertools import product
    from time import time

    import numpy as np

    from comando.utility import lambdify, evaluate  # , get_vars

    # DATA GENERATION
    # func = lambdify(expr)
    # vars = get_vars(expr)

    # Generating validation points and data
    validation_points = {}
    for v in vars:
        validation_points[v] = np.linspace(*v.bounds, points_per_variable)
    # evaluate the function on a grid and store the results
    n_in = len(vars)

    if callable(data_or_func):
        func = data_or_func
        data = {vals: func(*vals)
                for vals in product(*(validation_points[v] for v in vars))}
    else:
        data = data_or_func

        def func(*vals):
            return data[vals]
    lb, ub = min(data.values()), max(data.values())

    # Actual ANN in domain [lb, ub]
    ann = make_fit_ann([n_in, *hidden_layer_neurons, 1],
                       [v.bounds for v in vars], [[lb, ub]])

    # Scaling from [lb, ub] to [-1, 1]
    delta = ub - lb
    stretch = 2 / delta
    offset = lb * stretch + 1

    # NOTE: scaling and unscaling cannot be used directly for errors as those
    #       present differences, so the offset is eliminated and we only need
    #       to multiply with / divide by stretch!
    def scale(val):
        return val * stretch - offset

    def unscale(scaled_val):
        return (scaled_val + offset) / stretch

    # Sacled ANN for training, domain [-1, 1]
    scaled_ann = make_fit_ann([n_in, *hidden_layer_neurons, 1],
                              [v.bounds for v in vars])

    # For fitting original variables become parameters
    x_data = [Parameter(v.name) for v in vars]
    f_data = Parameter('f_data')
    model = scaled_ann(*x_data)
    # model = ann(*x_data)
    args = sorted(model.free_symbols, key=lambda s: s.name)
    print(f'# Neural network fit with {n_in} input/s and '
          f'{len(hidden_layer_neurons)} hidden layers, with\n'
          f'# the following nubers of neurons {hidden_layer_neurons}\n'
          f'def f({", ".join(str(s) for s in args)}):\n'
          f'    return {ann(*x_data)}')
    max_scaled_err = comando.Variable('max_scaled_err', bounds=(0, None))

    ERR = abs(model - f_data)
    # ERR = (model - f_data) ** 2
    # ERR = (1 - model / f_data - f_data / model) ** 2
    # Initial problem formulation with scaled data from validation points
    n_data = len(data)
    P = comando.Problem(max_scaled_err, timesteps=[range(n_data), n_data],
                        constraints={'max_err_bound':
                                     0 >= (ERR - max_scaled_err)})
    # SETTING DATA VALUES
    scaled_data = {k: scale(v) for k, v in data.items()}
    # if n_in == 1:
    #     P[x_data[0].name] = scaled_data.keys()
    # else:
    for i, v in enumerate(vars):
        P[v.name] = [k[i] for k in scaled_data.keys()]
    P['f_data'] = scaled_data.values()

    # max_err.value = max(comando.utility.evaluate(abs(model - f_data)))

    # from comando.interfaces.maingo_api import MaingoProblem
    # mp = MaingoProblem(P)
    # mp.solve(epsilonA=tol, maxTime=T_max)

    # from comando.interfaces.maingo_ale import solve
    # options = dict(epsilonA=tol, maxTime=T_max)
    # solve(P, 'ANN_TEST', reuse=False, **options)

    if solver.lower() == 'baron':
        from comando.interfaces.baron import solve as b_solve, baron_str_map
        baron_str_map['Abs'] = lambda arg: f'(({arg}) ^ 2) ^ 0.5'
        # NumLoc=-1
        # MaxIteroption=0
        # TDo, MDo, LBTTDo, andOBTTDo = 0
        MAXTIME = 'MaxTime'
        options = dict(epsa=tol, epsr=tol, MaxTime=iter_time, **solver_options)

        def solve(**options):
            b_solve(P, 'ANN_TEST', silent=True, cse=True, reuse=False,
                    **options)
    elif solver.lower() == 'maingo':
        from comando.interfaces.maingo_api import MaingoProblem
        mp = MaingoProblem(P)
        MAXTIME = 'maxTime'
        options = dict(epsilonA=tol, epsilonR=tol, maxTime=iter_time,
                       **solver_options)

        def solve(**options):
            mp.solve(**options)
    else:
        raise NotImplementedError('No routine is implemented for solver '
                                  + str(solver) + '!')

    # Initialize training data
    training_data = {} if initial_points is None \
        else {k: scale(v) for k, v in initial_points.items()}

    def scaled_func(*args):
        return scale(func(*args))

    if n_in == 1:  # Plot if 1D
        import matplotlib.pyplot as plt
        plt.figure()
        ax = plt.gca()
        x = next(iter(validation_points.values()))
        y = scaled_func(x)
        ax.plot(x, y, label='function')
        ax.scatter(x, y, marker='+', zorder=3, c='k', lw=0.5,
                   label='validation points')
        if training_data:
            ax.scatter(training_data.keys(), training_data.values(), lw=0.5,
                       zorder=3, c='k', marker='o', label='initial points')
    # elif n_in == 2:
    #     from mpl_toolkits.mplot3d import Axes3D
    #     import matplotlib.pyplot as plt
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection='3d')
    #     x, y = np.array([*scaled_data.keys()]).T
    #     z = np.array([*scaled_data.values()])
    #     X = x.reshape(points_per_variable, -1)
    #     Y = y.reshape(points_per_variable, -1)
    #     Z = z.reshape(points_per_variable, -1)
    #     ax.scatter(X, Y, Z, zorder=2, label='function')
    #     # ax.plot_surface(X, Y, Z, label='function')
    #     if training_data:
    #         ax.scatter([k[0] for k in training_data.keys()],
    #                    [k[1] for k in training_data.keys()],
    #                    np.array([*training_data.values()]), lw=0.5,
    #                    zorder=3, c='k', marker='o', label='initial points')
    #     plt.legend()
    parametrized_fit = ann(*vars)
    params = parametrized_fit.free_symbols - {*vars}
    best_fit = fit = parametrized_fit.subs({s: s.value for s in params})

    # def make_func(fit):
    #     required_vars = [i for s in fit.free_symbols
    #                      for i, v in enumerate(vars) if v is s]
    #     _fit_func = lambdify(fit)
    #
    #     def fit_func(*args):
    #         passed_args = [args[i] for i in required_vars]
    #         return _fit_func(*passed_args)
    #     return fit_func

    # def get_max_scaled_error(fit):
    #     # NOTE: As zero weights weights may eliminate some input variables we
    #     #       need to manually ensure the function is always callable with
    #     #       all vars!
    #     required_vars = [i for s in fit.free_symbols
    #                      for i, v in enumerate(vars) if v is s]
    #     _fit_func = lambdify(fit)
    #
    #     def fit_func(*args):
    #         passed_args = [args[i] for i in required_vars]
    #         return _fit_func(*passed_args)
    #
    #     max_scaled_error = 0
    #     point = None
    #     for val, res in scaled_data.items():
    #         err = np.abs(res - fit_func(*val))
    #         if err > max_scaled_error:
    #             max_scaled_error = err
    #             point = val
    #
    #     return max_scaled_error, point

    def get_max_error(fit):
        # NOTE: As zero weights weights may eliminate some input variables we
        #       need to manually ensure the function is always callable with
        #       all vars!
        required_vars = [i for i, v in enumerate(vars)
                         if v in fit.free_symbols]
        _fit_func = lambdify(fit, [v for i, v in enumerate(vars)
                                   if i in required_vars])

        def fit_func(*args):
            passed_args = [args[i] for i in required_vars]
            return _fit_func(*passed_args)

        max_error = 0
        point = None
        for val, res in data.items():
            err = np.abs(res - fit_func(*val))
            if err > max_error:
                max_error = err
                point = val

        return max_error, point

    max_error, point = get_max_error(fit)
    best_found = max_error
    max_scaled_err.ub = min(best_found, enforced_error_bound) * stretch
    i = 0
    fits = []
    if point is None:
        print('No point with an error found!')
        if filename:
            # ann.out_lbs[0] = lb
            # ann.out_delta[0] = delta
            ann.dump(filename)
        return fit, 0, data

    # LOOP
    t0 = time()
    while time() - t0 < T_max:
        i += 1
        print(f'Iteration {i}')
        # if n_in == 1:
        #     training_data[point[0]] = scaled_func(*point)
        #     P.timesteps = [range(len(training_data)), 1]
        #     P[x_data[0].name] = training_data.keys()
        # else:
        if point in training_data:
            print("Point of largest violation within previous set of points, "
                  "doubling iteration time!")
            iter_time += iter_time
        training_data[point] = scaled_func(*point)
        P.timesteps = [range(len(training_data)), 1]
        for j, v in enumerate(vars):
            P[v.name] = [k[j] for k in training_data.keys()]
        P['f_data'] = training_data.values()
        print(f'\toptimizing with {len(training_data)} points...')
        # Current maximum of the considered points only
        max_scaled_err.value = max(evaluate(ERR))

        # update time
        options[MAXTIME] = min(iter_time, T_max - (time() - t0))
        ret = solve(**options)
        if ret == -1:
            print("Couldn't find a better solution, doubling iteration time!")
            iter_time += iter_time
            # continue
        fit = parametrized_fit.subs({s: s.value for s in params})
        fits.append(fit)
        max_error, new_point = get_max_error(fit)
        if new_point is None:
            print('No point with an error found!')
            if filename:
                # ann.out_lbs[0] = lb
                # ann.out_delta[0] = delta
                ann.dump(filename)
            return fit, 0, data

        if max_error < best_found:  # Update
            best_found = max_error
            max_scaled_err.ub = min(best_found * stretch, max_scaled_err.ub)
            best_fit = fit
            if filename and best_found < enforced_error_bound:
                # Only update if we actually improved
                # ann.out_lbs[0] = lb
                # ann.out_delta[0] = delta
                ann.dump(filename)
                # ann.out_lbs[0] = -0.5
                # ann.out_delta[0] = 1.
        elif (max_error - best_found) / best_found < 0.01:
            print("Not making progress, doubling iteration time!")
            iter_time += iter_time
        print(f'\tbest error bound: {best_found}')
        print(f'\tmax error: {max_error}')
        if n_in == 1:
            if len(ax.lines) > 1:  # Removing the last fit
                # ax.lines[-1].remove()
                ...
                # ax.collections[-1].remove()
            else:
                ax.scatter(0, float('nan'), zorder=3, marker='o',
                           facecolors='none', edgecolors='k', lw=1,
                           label='tested points')
            from comando.visualization import plot_expr
            scaled_fit = scaled_ann(*vars).subs({s: s.value for s in params})
            plot_expr(scaled_fit, label=f'iteration {i} max_error = '
                      f'{max_error}', ls='dashed')
            ax.scatter(point, training_data[point], zorder=3, marker='o',
                       edgecolors=ax.lines[-1].get_color(), facecolors='none',
                       lw=1)
            plt.title(f'Result with {hidden_layer_neurons} hidden neurons\n'
                      f'best error bound = {best_found}')
            plt.legend()
            plt.pause(0.01)
        # if n_in == 2:
        #     scaled_fit = scaled_ann(*vars).subs({s: s.value for s in params})
        #     z = lambdify(scaled_fit, vars)(x, y)
        #     try:
        #         Z = z.reshape(points_per_variable, -1)
        #     except AttributeError:
        #         Z = np.ones((points_per_variable, points_per_variable)) * z
        #     ax.plot_surface(X, Y, Z,
        #                     label=f'iteration {i} max_error = {max_error}')
        #     plt.title(f'Result with {hidden_layer_neurons} hidden neurons\n'
        #               f'best error bound = {best_found}')
        #     # plt.legend()
        #     plt.pause(0.01)
        point = new_point
    else:
        print('Time limit reached!')

    return best_fit, best_found, {k: unscale(v)
                                  for k, v in training_data.items()}


def make_fit(*vars, data, hidden_layer_neurons, tol=1e-3, T_max=30,
             solver='baron', exact=None, verbose=False, filename=None,
             **options):
    """Create a scalar-valued function fit for the data using an ANN."""
    out_bounds = min(data.values()), max(data.values())
    if out_bounds[0] == out_bounds[1]:
        return comando.S(out_bounds[0]), 0  # constant function, no error
    n_in = len(vars)
    n_data = len(data)
    ann = make_fit_ann([n_in, *hidden_layer_neurons, 1],
                       [v.bounds for v in vars], [out_bounds])
    x_data = [Parameter(f'x_data_{i}') for i in range(n_in)]
    f_data = Parameter('f_data')
    model = ann(*x_data)
    args = sorted(model.free_symbols, key=lambda s: s.name)
    if verbose:
        print(f'# Neural network fit with {n_in} input/s and '
              f'{len(hidden_layer_neurons)} hidden layers, with\n'
              f'# the following nubers of neurons {hidden_layer_neurons}\n'
              f'def f({", ".join(str(s) for s in args)}:\n'
              f'    return {model}')
    # SD
    model_error = (model - f_data) ** 2

    # from comando.utility import smooth_abs
    # model_error = smooth_abs(model - f_data, delta=0)  # RMS
    fit = ann(*vars)
    if exact:
        cons = {var_vals: comando.Eq(fit.subs({var: var_val
                                               for var, var_val
                                               in zip(vars, var_vals)}),
                                     f_val)
                for var_vals, f_val in exact.items()}
    else:
        cons = None

    # DEBUG:
    P = Problem(0, model_error, timesteps=[range(n_data), n_data],
                constraints=cons, states={}, name='ANN_fitting_problem')
    if n_in == 1:
        P['x_data_0'] = data.keys()
    else:
        for i in range(n_in):
            P[f'x_data_{i}'] = [k[i] for k in data.keys()]
    P['f_data'] = data.values()

    # print(f'minimizing {P.objective}')
    # for s in P.objective.free_symbols:
    #     if s.is_Parameter:
    #         continue
    #     print(s, s.bounds)
    # # Maingo requires bounds, and in general appears not to be suitable for this task :(
    # from comando.interfaces.maingo_api import MaingoProblem
    # # P = Problem(SP, 0, {}, {}, [[1], 1])
    # mp = MaingoProblem(P)
    # # sol, stat = mp.solve(outstreamVerbosity=1)
    # sol, stat = mp.solve(outstreamVerbosity=3)
    # print(stat)
    # return sol

    # sol.get_solution_point()
    # print(sol.get_objective_value())
    if solver.lower() == 'baron':
        from comando.interfaces.baron import solve
        solve(P, 'ann_fit.bar', epsa=tol, MaxTime=T_max, reuse=False, **options)
    elif solver.lower() == 'maingo':
        from comando.interfaces.maingo_api import MaingoProblem
        mp = MaingoProblem(P)
        mp.solve(epsilonA=tol, maxTime=T_max, **options)
    elif solver == 'knitro':
        from comando.interfaces.pyomo import to_pyomo
        m = to_pyomo(P)
        m.solve(solver='knitro', remote=True,
                options=dict(feastol_abs=tol, maxtime_cpu=T_max))
    else:
        raise NotImplementedError("No procedure implemented for solver "
                                  f"'{solver}'!")
    expr = P.objective
    RMSE = expr.subs({sym: sym.value for sym in expr.free_symbols}) ** 0.5
    try:
        RMSE = float(RMSE)
    except:
        RMSE = float(RMSE.real)
    if verbose:
        print('RMSE:', RMSE)

    if filename:
        ann.dump(filename)
    return fit.subs({s: s.value for s in fit.free_symbols - {*vars}}), RMSE
