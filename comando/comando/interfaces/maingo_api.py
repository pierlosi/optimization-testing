"""API interface for MAiNGO."""
# This file is part of the COMANDO project which is released under the MIT
# license. See file LICENSE for full license details.
#
# AUTHORS: Marco Langiu

# NOTE: since Python 3.8 %PATH% is no longer searched for DLLs on Windows for
#       security reasons. This can result in failure to find DLLs required by
#       certain modules like maingopy. Since we're very trusting we just
#       recreate the behavior from Python versions prior to 3.8, by manually
#       including all directories from %PATH%.
import sys
if (sys.platform == 'win32'
        and sys.version_info.major >= 3
        and sys.version_info.minor >= 8):
    import os
    from packaging import version
    print('\nLooking for CPLEX versions...\n')
    CPLEX_versions = {
            version.parse(f'{k[16:18]}.{k[18:]}'): v[1:-1] if v.startswith('"') else v
        for k, v in os.environ.items()
        if k.startswith('CPLEX_STUDIO_DIR')
    }
    if len(CPLEX_versions) == 0:
        raise ImportError("CPLEX wasnt found, but maingopy appears to require "
                          "it for its import.")
    if len(CPLEX_versions) == 1:
        version, path = next(iter(CPLEX_versions.items()))
    else:
        print('Found multiple CPLEX versions:')
        for version, path in CPLEX_versions.items():
            print('  Version', version, 'at', path)
        version = max(CPLEX_versions)
        path = CPLEX_versions[version]

    print(f'using CPLEX {version} from {path}')

    os.add_dll_directory(f'{path}cplex/bin/x64_win64')
try:
    import maingopy
except ModuleNotFoundError:
    raise ModuleNotFoundError("I couldn't find the maingopy shared library "
                        "for the Python version you're using "
                        f"({sys.executable})!\nIf you know the path where I "
                        "can find it, you can add it to the PYTHONPATH "
                        "environment variable!")

FEASIBLE_POINT = maingopy.FEASIBLE_POINT
GLOBALLY_OPTIMAL = maingopy.GLOBALLY_OPTIMAL
LANG_ALE = maingopy.LANG_ALE
LANG_GAMS = maingopy.LANG_GAMS
LANG_NONE = maingopy.LANG_NONE

import comando
from comando import REAL, INTEGER, BINARY
from comando.utility import get_index

d_map = {
    REAL: maingopy.VT_CONTINUOUS,
    INTEGER: maingopy.VT_BINARY,
    BINARY: maingopy.VT_INTEGER
}

from functools import reduce

maingo_API_op_map = {
    'Pow': maingopy.pow,
    'acos': maingopy.acos,
    'acquisition_function': maingopy.acquisition_function,
    'arh': maingopy.arh,
    'asin': maingopy.asin,
    'atan': maingopy.atan,
    'bounding_func': maingopy.bounding_func,
    'bstep': maingopy.bstep,
    'cheb': maingopy.cheb,
    'cos': maingopy.cos,
    'cosh': maingopy.cosh,
    'cost_function': maingopy.cost_function,
    'cost_turton': maingopy.cost_turton,
    'coth': maingopy.coth,
    'covariance_function': maingopy.covariance_function,
    'enthalpy_of_vaporization': maingopy.enthalpy_of_vaporization,
    'erf': maingopy.erf,
    'erfc': maingopy.erfc,
    'euclidean_norm_2d': maingopy.euclidean_norm_2d,
    'exp': maingopy.exp,
    'expx_times_y': maingopy.expx_times_y,
    'fabs': maingopy.fabs,
    'fabsx_times_x': maingopy.fabsx_times_x,
    'fstep': maingopy.fstep,
    'gaussian_probability_density_function':
        maingopy.gaussian_probability_density_function,
    'iapws': maingopy.iapws,
    'ideal_gas_enthalpy': maingopy.ideal_gas_enthalpy,
    'inv': maingopy.inv,
    'lb_func': maingopy.lb_func,
    'lmtd': maingopy.lmtd,
    'log': maingopy.log,
    'Max': lambda *args: reduce(maingopy.max, args),
    'Min': lambda *args: reduce(maingopy.min, args),
    'neg': maingopy.neg,
    'norm2': maingopy.norm2,
    'nrtl_G': maingopy.nrtl_G,
    'nrtl_Gdtau': maingopy.nrtl_Gdtau,
    'nrtl_Gtau': maingopy.nrtl_Gtau,
    'nrtl_dGtau': maingopy.nrtl_dGtau,
    'nrtl_dtau': maingopy.nrtl_dtau,
    'nrtl_tau': maingopy.nrtl_tau,
    'p_sat_ethanol_schroeder': maingopy.p_sat_ethanol_schroeder,
    'pos': maingopy.pos,
    'pow': maingopy.pow,
    'regnormal': maingopy.regnormal,
    'rho_liq_sat_ethanol_schroeder': maingopy.rho_liq_sat_ethanol_schroeder,
    'rho_vap_sat_ethanol_schroeder': maingopy.rho_vap_sat_ethanol_schroeder,
    'rlmtd': maingopy.rlmtd,
    'saturation_temperature': maingopy.saturation_temperature,
    'sin': maingopy.sin,
    'sinh': maingopy.sinh,
    'sqr': maingopy.sqr,
    'sqrt': maingopy.sqrt,
    'squash_node': maingopy.squash_node,
    'sum_div': maingopy.sum_div,
    'tan': maingopy.tan,
    'tanh': maingopy.tanh,
    'ub_func': maingopy.ub_func,
    'vapor_pressure': maingopy.vapor_pressure,
    'xexpax': maingopy.xexpax,
    'xlog': maingopy.xlog,
    'xlog_sum': maingopy.xlog_sum
}


def floor_substitute(x, LB, UB):
    from numpy import floor
    lb = int(floor(LB))
    ub = int(floor(UB))
    return maingopy.max(x - 1, lb + sum(maingopy.fstep(x - i)
                                        for i in range(lb + 1, ub + 1)))

maingo_API_op_map['floor_substitute'] = floor_substitute


def MAiNGO_var(v, branching_priorities):
    """Create a MAiNGO variable from a COMANDO variable."""
    try:  # if v is an element of an indexed variable, default to parent prio
        prio = branching_priorities.get(v.parent, 1)
        # prio = branching_priorities.get(v, prio))  # own or parent
        prio = max(branching_priorities.get(v, 1), prio)  # use max
        # prio = branching_priorities.get(v, 0) + prio  # use sum
    except AttributeError:
        prio = branching_priorities.get(v, 1)  # get own priority
    return maingopy.OptimizationVariable(maingopy.Bounds(*v.bounds),
                                         d_map[v.domain], prio, v.name)


class MaingoProblem(maingopy.MAiNGOmodel):
    """Initialize a MAiNGO Model from a COMANDO Problem.

    Arguments
    ---------
    P : comando.Problem
        A COMANDO problem that is to be translated to MAiNGO data structures.
    relaxation_only : set
        constraint ids that are to be treated as relaxation only
    squash : set
        constraint ids that are to be treated as relaxation only
    outputs : Mapping
        Mapping from textual description to expressions of interest for output
    use_cse : bool (default: True)
        Whether to use common subexpression elimination to represent
        reoccurring expressions with intermediate variables
    """

    def __init__(self, P, relaxation_only=None, squash=None, output=None,
                 use_cse=True):
        super().__init__()
        # NOTE: it would make sense if the OptimizationVariable objects and the
        #       vars passed to the evaluate function were the same objects!
        #       This would allow to build a symbol mapping once and just store
        #       the MAiNGO equivalent of objective and constraints.
        # self.sym_map = {}
        # for p in P.parameters:  # translate parameters (scalars and vectors)
        #     try:  # treating the parameter as a vector
        #         mpars = {}
        #         for i, pi in p.expansion.items():
        #             mpars[pi] = self.sym_map[pi] = pi.value
        #         self.sym_map[p] = mpars
        #     except AttributeError:  # Parameter is scalar
        #         self.sym_map[p] = p.value
        # for v in P.design_variables:  # translate scalar variables
        #     sym_map[v] = MAiNGO_var(v)
        # for vv in P.operational_variables:
        #     mvars = {}
        #     for i, v in vv.expansion.items():  # translate vector variables
        #         mvars[v] = self.sym_map[v] = MAiNGO_var(v)
        #     self.sym_map[vv] = mvars

        # Then we could just do all parsing once at this point and be done...:
        # self.obj = parse(P.objective, sym_map)
        # self.cons = {c_id: parse(c, sym_map) for c_id, c in P.constraints}
        # self.states = ...  # TODO: Differential state discretization

        # ... But that's not the case, so we define an order via generators...
        self.P = P

        # def yield_pars():
        #     for p in P.parameters:
        #         if p.indexed:
        #             yield from p
        #         else:
        #             yield p
        # self.yield_pars = yield_pars

        from operator import attrgetter
        name_attr = attrgetter('name')
        self.dvs = sorted(self.P.design_variables, key=name_attr)
        self.ovs = sorted(self.P.operational_variables, key=name_attr)

        def yield_vars():
            yield from self.dvs
            for ov in self.ovs:
                yield from ov.expansion
        self.yield_vars = yield_vars
        # ... and delay sym map and parser creation to _make_sym_map!

        self.relaxation_only = set() if relaxation_only is None \
            else relaxation_only
        self.squash = set() if squash is None else squash
        self.output = dict() if output is None else output

        self.use_cse = use_cse

    def _make_sym_map(self, vars):
        print('_make_sym_map')

        self.sym_map = {}
        for p in self.P.parameters:  # translate parameters (scalars & vectors)
            try:  # treating the parameter as a vector
                mpars = {}
                for i, pi in p.expansion.items():
                    mpars[i] = self.sym_map[pi] = maingopy.FFVar(pi.value)
                self.sym_map[p] = mpars
            except AttributeError:  # Parameter is scalar
                self.sym_map[p] = maingopy.FFVar(p.value)
        # NOTE: Neglects valiable vecor mappings!
        # self.sym_map.update({v: mv for v, mv in zip(self.yield_vars(), vars)})
        var_iter = iter(vars)
        for v in self.dvs:  # translate scalar variables
            self.sym_map[v] = next(var_iter)
        for vv in self.ovs:
            mvars = {}
            for i, v in vv.expansion.items():  # translate vector variables
                mvars[i] = self.sym_map[v] = next(var_iter)
            self.sym_map[vv] = mvars

        def parse(expr, i=None):
            """Parse from COMANDO to MAiNGO."""
            try:
                return comando.utility.parse(expr, self.sym_map,
                                             maingo_API_op_map, i)
            except Exception as e:
                print('\nParsing to MAiNGO raised', e)

                print('\nindex was', i)

                print('\nexpression was:')
                print(expr)

                print('\nSymbols were mapping to:')
                for sym in expr.free_symbols:
                    print(sym, '->', comando.utility.parse(sym, self.sym_map,
                                                           maingo_API_op_map, i))

                raise
        self.parse = parse

    def get_variables(self):
        """Get the MAiNGO equivalent of COMANDO variables."""
        return [MAiNGO_var(v, self.branching_priorities)
                for v in self.yield_vars()]

    def get_initial_point(self):
        """Get the current variable values as the initial point."""
        return [v.value for v in self.yield_vars()]

    def evaluate(self, vars):
        """Get an evaluation container representing objective and constraints.

        The evaluation container contains entries for the following
        expressions:

        * obj: objective
        * ineq: inequality constraints
        * eq: equality constraints
        * ineqSquash: squash inequality constraints (vegans only!?)
        * ineqRO: relaxation-only inequality constraints
        * ineqRO: relaxation-only inequality constraints
        * eqRO: relaxation-only equality constraints
        * out: output

        All of these expressions contain FFVar Objects for the variables (NOT
        OptimizationVariable objects!).
        The `out` entry is a list of OutputVariable Objects for display.
        OutputVariables combine an expression in the form of an FFVar object
        with a textual description.
        """

        print('EVALUATE')
        self._make_sym_map(vars)

        P = self.P
        if self.use_cse is True:
            # Using Common Subexpression Elimination to determine intermediate
            # variables for reoccurring subexpressions within the problem
            reps, exprs = comando.cse((P.design_objective,
                                       P.operational_objective,
                                       *P.constraints.values()))

            defs = {}
            # self = mp; sym, rep = next(iter(reps)); i = 'nominal'
            for sym, rep in reps:
                e = rep.subs(defs)
                index = get_index(e)
                if index is None:
                    x = comando.Variable(sym.name)
                    self.sym_map[x] = self.parse(e)
                else:
                    x = comando.VariableVector(sym.name)
                    x.instantiate(index)
                    mreps = {}
                    for i, xi in x.expansion.items():
                        mreps[i] = self.sym_map[xi] = self.parse(e, i)
                    self.sym_map[x] = mreps
                defs[sym] = x

            print('CSE_END')

            obj = exprs[0].subs(defs) + P.weighted_sum(exprs[1].subs(defs))
            cons = {con_id: expr.subs(defs) for con_id, expr in
                    zip(P.constraints, exprs[2:])}
        else:
            obj, cons = P.objective, P.constraints

        print('EvaluationContainer')

        result = maingopy.EvaluationContainer()
        result.obj = self.parse(obj)
        for c_id, con in cons.items():
            try:  # treating con as ≤ or ≥ side
                norm_con = con.lts - con.gts
                if c_id in self.relaxation_only:
                    entry = result.ineqRO
                elif c_id in self.squash:
                    entry = result.ineqSquash
                else:
                    entry = result.ineq
            except AttributeError:  # Eq constraint
                norm_con = con.lhs - con.rhs
                if c_id in self.relaxation_only:
                    entry = result.ineqRO
                else:
                    entry = result.eq
            index = get_index(norm_con)
            if index is None:
                entry.push_back(self.parse(norm_con), c_id)
            else:
                for i in index:
                    entry.push_back(self.parse(norm_con, i), f'{c_id}_{i}')

            for descr, expr in self.output.items():
                index = get_index(expr)
                if index is None:
                    result.out.append(
                        maingopy.OutputVariable(descr, self.parse(expr)))
                else:
                    result.out.extend(
                        [maingopy.OutputVariable(descr,  self.parse(expr, i))
                         for i in index])

        return result

    def _check_prios(self, branching_priorities):
        vars = {*self.yield_vars(), *self.P.operational_variables}
        for var, val in branching_priorities.items():
            if var not in vars:
                raise ValueError(f'Entry "{var}" is neither a design, nor an '
                                 f'operational variable of "{self.P.name}"')
        return branching_priorities

    def solve(self, **options):
        """Solve the problem with the given options.

        Options
        -------
        epsilonA : double
            Absolute optimality tolerance, i.e., termination when (UBD-LBD) <
            BAB_epsilon_a.

        epsilonR : double
            Relative optimality tolerance, i.e., termination when (UBD-LBD) <
            BAB_epsilon_r * UBD.

        deltaIneq : double
            Absolute feasibility tolerance for inequality constraints (i.e.,
            constraint is considered satisfied if gi_(x)<=UBP_delta_ineq.

        deltaEq : double
            Absolute feasibility tolerance for equality constraints (i.e.,
            constraint is considered satisfied if abs(hi_(x))<=UBP_delta_eq.

        relNodeTol : double
            Relative tolerance for minimum node size.

        BAB_maxNodes : unsigned
            Maximum number of nodes (i.e., solver terminates when more than
            BAB_maxnodes are held in memory; used to avoid excessive branching)

        BAB_maxIterations : unsigned
            Maximum number of iterations (i.e., maximum number of nodes visited
            in the Branch-and-Bound tree)

        maxTime : unsigned
            CPU time limit in seconds.

        confirmTermination : bool
            Whether to ask the user before terminating when reaching time,
            node, or iteration limits.

        terminateOnFeasiblePoint : bool
            Whether to terminate as soon as the first feasible point was found
            (no guarantee of global or local optimality!)

        targetLowerBound : double
            Target value for the lower bound on the optimal objective. MAiNGO
            terminates once LBD>=targetLowerBound (no guarantee of global or
            local optimality!)

        targetUpperBound : double
            Target value for the upper bound on the optimal objective. MAiNGO
            terminates once UBD<=targetUpperBound (no guarantee of global or
            local optimality!)

        infinity : double
            User definition of infinity (used to initialize UBD and LBD)
            [currently cannot be set by the user via set_option].

        PRE_maxLocalSearches : unsigned
            Number of local searches in the multistart heuristic during
            preprocessing at the root node.

        PRE_obbtMaxRounds : unsigned
            Maximum number of rounds of optimization-based range reduction
            (OBBT; cf., e.g., Gleixner et al., J. Glob. Optim. 67 (2017) 731;
            maximizing and minimizing each variable subject to relaxed
            constraints) at the root node. If >=1 and a feasible point is
            found during multistart, one round of OBBT using an objective cut
            (f_cv<=UBD) is conducted as well.

        PRE_pureMultistart : bool
            Whether to perform a multistart only. A B&B tree will not be
            constructed and no lower bounding problems will be solved.

        BAB_nodeSelection : babBase::enums::NS
            How to select the next node to process. See documentation of
            babBase::enums::NS for possible values.

        BAB_branchVariable : babBase::enums::BV
            Which dimension to branch in for the current node. See
            documentation of babBase::enums::BV for possible values.

        BAB_alwaysSolveObbt : bool
            Whether to solve OBBT (feasibility- and, once a feasible point has
            been found, also optimality-based) at every BaB node.

        BAB_dbbt : bool
            Whether to do a single round of duality based bound tightening
            (DBBT, cf. Ryoo&Sahinidis, Comput. Chem. Eng. 19 (1995) 551). If
            false, no DBBT is used. If true, multipliers from CPLEX are used to
            tighten bounds (essentially for free). we tried additional rounds
            but without reasonable improvement.

        BAB_probing : bool
            Whether to do probing (cf. Ryoo&Sahinidis, Comput. Chem. Eng. 19
            (1995) 551) at every node (can only be done if BAB_DBBT_maxrounds
            >= 1)

        BAB_constraintPropagation : bool
            Whether to do constraint propagation. If false, no constraint
            propagation is executed.

        LBP_solver : lbp::LBP_SOLVER
            Solver for solution of lower bounding problems.

        LBP_linPoints : lbp::LINP
            At which points to linearize for affine relaxation. See
            documentation of lbp::LINP for possible values.

        LBP_subgradientIntervals : bool
            Whether to use the heuristic to improve McCormick relaxations by
            tightening the range of each factor with the use of subgradients
            (cf. Najman & Mitsos, JOGO 2019)

        LBP_obbtMinImprovement : double
            How much improvement needs to be achievable (relative to initial
            diameter) to conduct OBBT for a variable.

        LBP_activateMoreScaling : unsigned
            Number of consecutive iterations without LBD improvement needed to
            activate more aggressive scaling in LP solver (e.g., CPLEX)

        LBP_addAuxiliaryVars : bool
            Whether to add auxiliary variables for common factors in the lower
            bounding DAG/problem.

        LBP_minFactorsForAux : unsigned
            Minimum number of common factors to add an auxiliary variable.

        LBP_maxNumberOfAddedFactors : unsigned
            Maximum number of added factor as auxiliaries.

        MC_mvcompUse : bool
            Whether to use multivariate composition theorem for computing
            McCormick relaxations (see MC++ documentation for details)

        MC_mvcompTol : double
            (see MC++ documentation for details)

        MC_envelTol : double
            (see MC++ documentation for details)

        UBP_solverPreprocessing : ubp::UBP_SOLVER
            Solver to be used during pre-processing (i.e., multistart). See
            documentation of ubp::UBP_SOLVER for possible values.

        UBP_maxStepsPreprocessing : unsigned
            Maximum number of steps the local solver is allowed to take in each
            local run during multistart in pre-processing.

        UBP_maxTimePreprocessing : double
            Maximum CPU time the local solver is allowed to take in each local
            run during multistart in pre-processing. Usually, this should only
            be a fall-back option to prevent truly getting stuck in local
            solution.

        UBP_solverBab : ubp::UBP_SOLVER
            Solver to be used during Branch-and-Bound. See documentation of
            ubp::UBP_SOLVER for possible values.

        UBP_maxStepsBab : unsigned
            Maximum number of steps the local solver is allowed to take at each
            BaB node.

        UBP_maxTimeBab : double
            Maximum CPU time the local solver is allowed to take at each BaB
            node. Usually, this should only be a fall-back option to prevent
            truly getting stuck in local solution.

        UBP_ignoreNodeBounds : bool
            Flag indicating whether the UBP solvers should ignore the box
            constraints of the current node during the B&B (and consider only
            the ones of the root node instead).

        EC_nPoints : unsigned
            Number of points on the Pareto front to be computed in
            epsilon-constraint method (only available via the C++ API)

        BAB_verbosity : VERB
            How much output to print from Branch & Bound solver. Possible
            values are VERB_NONE (=0), VERB_NORMAL (=1), VERB_ALL (=2)

        LBP_verbosity : VERB
            How much output to print from Lower Bounding Solver. Possible
            values are VERB_NONE (=0), VERB_NORMAL (=1), VERB_ALL (=2)

        UBP_verbosity : VERB
            How much output to print from Upper Bounding Solver. Possible
            values are VERB_NONE (=0), VERB_NORMAL (=1), VERB_ALL (=2)

        BAB_printFreq : unsigned
            After how many iterations to print progress on screen
            (additionally, a line is printed when a new incumbent is found)

        BAB_logFreq : unsigned
            Like BAB_printFreq, but for log.

        writeLog : bool
            Whether to write a log file (named bab.log)

        writeToLogSec : unsigned
            Write to log file after a given ammount of CPU seconds.

        writeResFile : bool
            Whether to write an additional file containing non-standard
            information about the solved model.

        writeCsv : bool
            Whether to write a csv-log file (named bab.csv). Currently, this
            only include time, LBD, UBD, and final output.

        PRE_printEveryLocalSearch : bool
            Whether to print every run during multistart at the root node.

        writeToOtherLanguage : PARSING_LANGUAGE
            Write to a file in a different modeling language.

        File name options added for convenience
        ---------------------------------------

        iterations_csv_file_name, json_file_name, log_file_name,
        result_file_name, solution_and_statistics_csv_file_name : str
            names for the respective files generated by MAiNGO, paths are
            interpreted as relative to the current working directory.

        Returns
        -------
        solver : MAiNGO solver object
            A solver object that can be queried for solve related information
            and adjust different settings:

            * solver.evaluate_additional_outputs_at_point(point)
            * solver.evaluate_additional_outputs_at_solution_point()
            * solver.evaluate_model_at_point(point)
            * solver.evaluate_model_at_solution_point()
            * solver.get_LBP_count()
            * solver.get_UBP_count()
            * solver.get_cpu_solution_time()
            * solver.get_final_LBD()
            * solver.get_final_abs_gap()
            * solver.get_final_rel_gap()
            * solver.get_iterations()
            * solver.get_max_nodes_in_memory()
            * solver.get_objective_value()
            * solver.get_solution_point()
            * solver.get_status()
            * solver.get_wallclock_solution_time()
            * solver.read_settings('settings.txt')
            * solver.set_iterations_csv_file_name('iterations.csv')
            * solver.set_json_file_name('results.json')
            * solver.set_log_file_name('results.log')
            * solver.set_model(myMAiNGOmodel)
            * solver.set_option(option, value)
            * solver.set_result_file_name('res.txt')
            * solver.set_solution_and_statistics_csv_file_name('sol.csv')
            * solver.solve()
            * solver.write_model_to_file_in_other_language('ALE', 'prob.ale')

        status : MAiNGO RETCODE
            Return code for the solution, possible values are:

            * GLOBALLY_OPTIMAL
            * INFEASIBLE
            * FEASIBLE_POINT
            * NO_FEASIBLE_POINT_FOUND
            * BOUND_TARGETS
            * NOT_SOLVED_YET
            * JUST_A_WORKER_DONT_ASK_ME
        """
        self.branching_priorities = \
            self._check_prios(options.pop('branching_priorities', {}))
        solver = maingopy.MAiNGO(self)
        lang = options.pop('writeToOtherLanguage', LANG_NONE)
        if lang is None:
            lang = LANG_NONE
        if lang not in {LANG_ALE, LANG_GAMS, LANG_NONE}:
            try:  # whether a string was given
                lang = globals().get(f'LANG_{lang.upper()}')
            except KeyError:
                raise ValueError(f'Language {lang} is not implemented! '
                                 'Possible values for writeToOtherLanguage are'
                                 ' ALE, GAMS or NONE!')
        if lang != LANG_NONE:
            ending = {LANG_ALE: '.ale', LANG_GAMS: '.gms'}[lang]
            solver.write_model_to_file_in_other_language(lang,
                                                         self.P.name + ending)

        # Handle special options for adjusting default file names
        for file_name_option in [
            'iterations_csv_file_name',
            'json_file_name',
            'log_file_name',
            'result_file_name',
            'solution_and_statistics_csv_file_name',
        ]:
            file_name = options.pop(file_name_option, '')
            if file_name:
                getattr(solver, 'set_' + file_name_option)(file_name)

        for option, value in options.items():
            if not solver.set_option(option, value):
                raise ValueError(f'Option "{option}" is not recognized!')
        status = solver.solve()
        if status in {FEASIBLE_POINT, GLOBALLY_OPTIMAL}:
            for var, val in zip(self.yield_vars(),
                                solver.get_solution_point()):
                var.value = val
        return solver, status
