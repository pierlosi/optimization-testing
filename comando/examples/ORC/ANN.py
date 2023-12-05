"""Code to create functions corresponding to artificial neural networks."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu
import numpy as np


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
        return (np.array(data) - self.in_lbs) / self.in_delta * 2 - 1

    def output_transform(self, data):
        """Transform from the interval [-1, 1] to output data."""
        return ((np.array(data) + 1) / 2) * self.out_delta + self.out_lbs

    def __call__(self, *inputs):
        """Compute the outputs of the FFNN with given inputs."""
        inputs = self.input_transform(inputs)
        for layer in self.layers:
            inputs = np.array([neuron(inputs) for neuron in layer])
        return self.output_transform(inputs)


class ScalarFFNN(FFNN):
    """An ANN with only a single output."""

    def __call__(self, *inputs):
        """Compute the output of the ScalarFFNN with given inputs."""
        return FFNN.__call__(self, *inputs)[0]


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
          'T__sat': ('sat', 'p', 'T_sat')}


def make_ann(quantity):
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

    if quantity in MODELS:
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
    return ScalarFFNN(weight_arrays, vals('BW', [1, 3, 5]), af_lists,
                      [*zip(*vals('bounds', [0, 1]))],
                      [*zip(*vals('bounds', [2, 3]))])
