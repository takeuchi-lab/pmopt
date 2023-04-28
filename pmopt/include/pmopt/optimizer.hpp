#pragma once

#include <chrono>
#include <random>
#include "pmopt/optimization/basic.hpp"
#include "pmopt/optimization/screening.hpp"
#include "pmopt/pruners/multi_safe_pruner.hpp"
#include "pmopt/utils/containers.hpp"
#include "pmopt/logger.hpp"
#include "pmopt/miner.hpp"

namespace PMOpt
{

    template< class Enumerator, class Loss, class XGenerator, class WGenerator >
    struct Optimizer
    {
        using Symbol = typename Enumerator::Symbol;

        const std::string LOG_PREFIX = ">> ";
        const static unsigned int LOG_PRECISION = 6;
        const static unsigned int LOGW = 6;
        constexpr static double ZERO_TOL = 1e-12;
        constexpr static double ZERO_TOL_GAP = 1e-8;
        constexpr static double GAP_INF = 1e+10;

        Miner<Enumerator> & _miner;
        Activeset<Symbol> _activeset;
        std::mt19937 _random_generator;
        std::vector<unsigned int> _active_indices;
        std::vector<Variables> _variables;
        unsigned int _maxiter = 1000;
        unsigned int _init_dynamic_screen = 0;
        unsigned int _multi_dynamic_screen = 0;
        unsigned int _screening_freq = 10;
        unsigned int _screening_freq_init = 2;
        double _tolerance = 1e-4;
        double _time;
        bool _allow_multi_screen = true;
        bool _allow_multi_update = false;
        bool _sequential_update = false;


        /**
         * @brief Construct a new Optimizer object
         * @param miner 
         * @param seed
         */
        Optimizer(
            Miner<Enumerator> & miner,
            unsigned int seed
        ) : _miner(miner), _random_generator(seed) {}


        /**
         * @brief optimize
         * @param[out] solution
         * @param references 
         * @param params 
         * @return bool
         */
        bool optimize(
            Reference & solution,
            std::vector<const Reference *> references,
            const Parameters & params
        );


        /**
         * @brief screen inactive features by mining algorithm (i.e. safe pattern pruning)
         * @param params
         * @param size number of current variables
         * @param iter the current number of iteration used for logging
         */
        void screen_mining(const Parameters &, unsigned int, unsigned int) noexcept;


        /**
         * @brief screen inactive features by linear search (i.e. safe sceening)
         * @param params
         * @param size number of current variables
         * @param iter the current number of iteration used for logging
         */
        void screen(const Parameters &, unsigned int, unsigned int) noexcept;


        /**
         * @brief search max score by mining algorithm
         * @param vars
         * @param params
         * @return double 
         */
        double search_max_mining(const Variables &, const Parameters &) noexcept;


        /**
         * @brief search max score by linear search
         * @param[out] vars not const because score is also computed
         * @param params
         * @return double 
         */
        double search_max(Variables &, const Parameters &) const noexcept;


        /**
         * @brief compute norm of dual problem by mining algorithm
         * @param vars
         * @param params
         * @return double 
         */
        double compute_dual_norm_mining(const Variables &, const Parameters &) noexcept;


        /**
         * @brief compute norm of dual problem by linear search
         * @param[out] vars
         * @param params
         * @return double 
         */
        double compute_dual_norm(Variables &, const Parameters &) const noexcept;


        /**
         * @brief compute gap from primal/dual loss, norm , e.t.c.
         * @param[out] vars
         * @param params
         * @return double 
         */
        double compute_gap(Variables &, const Parameters &) const noexcept;


        /**
         * @brief update primal and dual variables
         * @param[out] vars,
         * @param params
         */
        void update(Variables &, const Parameters &) noexcept;


        /**
         * @brief convert iteration to string in order to log progress
         * @param iter 
         * @return std::string 
         */
        std::string logit(unsigned int iter) const noexcept
        {
            std::stringstream ss;
            ss << " [Iter: " << std::setw(std::log10(_maxiter) + 1) << iter << "]";
            return ss.str();
        }


        /**
         * @brief convert iteration to string in order to log progress
         * @param iter 
         * @return std::string 
         */
        std::string logit(const std::string & iter) const noexcept
        {
            std::stringstream ss;
            ss << ' ' << std::setw(std::log10(_maxiter) + std::string("[Iter: ]").size() + 1) << iter;
            return ss.str();
        }
    };
}


template< class Enumerator, class Loss, class XGenerator, class WGenerator >
bool PMOpt::Optimizer<Enumerator, Loss, XGenerator, WGenerator>::optimize(
    Reference & solution,
    std::vector<const Reference *> references,
    const Parameters & params
)
{
    // flow of optimization
    //
    // input
    //  - current coef
    //  - current intercept
    //  - current dual_raw (maybe infeasible)
    //  - weigth of regularization for each coef
    //
    // procedure:
    //  1. initialize variables
    //  2. compute max score used to scale dual variable in order to make feasible
    //  3. compute duality gap
    //  4. check convergence
    //  5. screen feature
    //  6. update variables which are not screened
    //  7. go to 2.
    //
    // note
    //  mining algorithm is used only in the first iteration
    //  in later iteration, simple linear search is done

    using namespace std::chrono;
    using seconds = std::chrono::duration<double>;
    auto start = system_clock::now();
    auto n_vars = references.size();
    _variables.resize(std::max(_variables.size(), (std::size_t) n_vars));

    Logger logger(std::cout);

    // for each reference 
    // compute duality gap, this will be used for screening
    double gap_min = GAP_INF;
    unsigned int index_gap_min = 0;
    for (unsigned int i_vars = 0; i_vars < n_vars; ++i_vars)
    {
        // construct variables from reference
        auto & vars = _variables[i_vars];
        vars.load(* references[i_vars]);

        // update intercept to make optimal for intercept
        // new intercept
        double intercept;
        if (vars._coef.empty())
            intercept = Loss::intercept_init(vars, params);
        else
            intercept = Loss::intercept(vars, params);

        // update dual variables
        double optimality = 0.0;
        for (auto i_instance : * params._instance_indices)
        {
            vars._y_pred[i_instance] += intercept - vars._intercept; // no masking
            optimality += Loss::dual(vars, params, i_instance) * params._mask[i_instance];
        }

        assert(std::abs(optimality) < ZERO_TOL);

        // update primal variable
        vars._intercept = intercept;


        logger(Logger::INFO+1, logit(0), " mining scores ...");


        // if lambda_l2 == 0 then search maximum score in order to scale dual variables
        // else compute the norm of the dual objective
        // TODO
        // is it necessary to mined lammax multiple times ????
        if (params._l2_scale || params._lambda_l2 < ZERO_TOL)
            vars._max_score = search_max_mining(vars, params);
        else
            vars._norm_dual = compute_dual_norm_mining(vars, params);


        logger(Logger::INFO+1, logit(0), " mining finished. max: ", vars._max_score, " sum: ", vars._norm_dual);
        logger(Logger::INFO, logit(0), " mining time: ", _miner._time, " [sec]");


        // evaluate loss function (with scaling dual variables)
        // compute: vars._loss_primal vars._loss_dual, vars._scale
        Loss::evaluate(vars._loss_primal, vars._loss_dual, vars._scale, vars, params);

        // compute primal regularization term with the current variable
        // compute: vars._norm_l1, vars._norm_l2
        Optimization::norm(vars._norm_l1, vars._norm_l2, vars);

        // compute duality gap
        // compute: vars._obj_primal/dual, vars._gap, vars._radius, vars._relative_gap, vars._norm_primal
        auto gap = compute_gap(vars, params);

        // check minimum gap
        if (gap < gap_min)
        {
            gap_min = gap;
            index_gap_min = i_vars;
        }


        logger(Logger::INFO, std::scientific, std::setprecision(LOG_PRECISION),
            logit(0), ' ', i_vars+1, '/', n_vars, " #A: ", std::setw(LOGW), "N/A", " relgap: ", gap);
        logger(Logger::INFO+1, vars.str(LOG_PRECISION));


        if (gap < - ZERO_TOL_GAP)
        {
            std::cerr << "Error: Negative gap detected: " << gap << "\n";
            exit(1);
        }
    }


    // check convergence
    if (gap_min < _tolerance) // converged only with intercept update
    {
        logger(Logger::INFO, logit(0), ' ', index_gap_min+1, '/', n_vars,
            " #A: ", std::setw(LOGW), std::setw(LOGW), _active_indices.size(),
            " gap: ", _variables[index_gap_min]._relative_gap);

        auto end = system_clock::now();
        _time = std::max(duration_cast<seconds>(end - start).count(), 0.0);

        solution = std::move(_variables[index_gap_min].dump(* references[index_gap_min]));
        return true;
    }


    // narrow down to one reference if multi screening is not allowed
    if (n_vars > 1 && !_allow_multi_screen)
    {
        logger(Logger::INFO, logit(0), " selected: ", index_gap_min + 1, '/', n_vars);

        std::swap(_variables[0], _variables[index_gap_min]);
        std::swap(references[0], references[index_gap_min]);
        n_vars = 1;
    }


    logger(Logger::INFO+1, logit(0), " mining activeset ...");


    // safe pattern pruning
    screen_mining(params, n_vars, 0);


    logger(Logger::INFO, logit(0), " mining finished. #A: ", _activeset._keys.size());
    logger(Logger::INFO, logit(0), " mining time: ", _miner._time, " [sec]");


    // if only use one variables for optimization
    if (n_vars > 1 && !_allow_multi_update)
    {
        logger(Logger::INFO, logit(0), " selected: ", index_gap_min + 1, '/', n_vars);

        std::swap(_variables[0], _variables[index_gap_min]);
        std::swap(references[0], references[index_gap_min]);
        n_vars = 1;
    }

    // set indices and coefficient values to variables
    for (unsigned int i_vars = 0; i_vars < n_vars; ++i_vars)
        _variables[i_vars].reload(_activeset, * references[i_vars]);

    // run optimization ...
    bool is_optimal = false;
    unsigned int iter = 1;
    unsigned int n_dynamic = 0;
    for (; iter < _maxiter; ++iter)
    {
        // model parameter update
        for (unsigned int i_vars = 0; i_vars < n_vars; ++i_vars)
            update(_variables[i_vars], params);

        // screening
        // screening and convergence check is done every _screening_freq iterations
        // screening will be done with the frequency _screening_freq if iter > _maxiter_multi
        // or with the frequency _screening_freq_multi if iter <= _maxiter_multi
        bool is_initial = n_dynamic < _init_dynamic_screen;
        if (is_initial && iter % _screening_freq_init != 0)
            continue;
        if (!is_initial && iter % _screening_freq != 0)
            continue;

        // compute gap for each variable
        // it is used in safe screening
        gap_min = GAP_INF;
        for (unsigned int i_vars = 0; i_vars < n_vars; ++i_vars)
        {
            auto & vars = _variables[i_vars];

            // compute score linearly 
            // if lambda_l2 is zero then search max in order to scale the dual variables
            // else compute the norm of the dual objective
            if (params._l2_scale || params._lambda_l2 < ZERO_TOL)
                vars._max_score = search_max(vars, params);
            else
                vars._norm_dual = compute_dual_norm(vars, params);

            // evaluate loss function (with scaling dual variables)
            // compute: vars._loss_primal vars._loss_dual, vars._scale
            Loss::evaluate(vars._loss_primal, vars._loss_dual, vars._scale, vars, params);

            // compute primal regularization term with the current variable
            // compute: vars._norm_l1, vars._norm_l2
            Optimization::norm(vars._norm_l1, vars._norm_l2, vars);

            // compute duality gap
            // compute: vars._obj_primal/dual, vars._gap, vars._radius, vars._relative_gap, vars._norm_primal
            auto gap = compute_gap(vars, params);

            // check minimum gap
            if (gap < gap_min)
            {
                gap_min = gap;
                index_gap_min = i_vars;
            }

            logger(Logger::INFO, logit(iter), ' ', i_vars+1, '/', n_vars,
                " #A: ", std::setw(LOGW), std::setw(LOGW), _active_indices.size(), " relgap: ", gap);
            logger(Logger::INFO+1, vars.str(LOG_PRECISION));

            if (gap < - ZERO_TOL_GAP)
            {
                std::cerr << "Error: Negative gap detected: " << gap << "\n";
                exit(1);
            }
        }
        ++n_dynamic;

        // check convergence
        is_optimal = gap_min < _tolerance;
        if (is_optimal)
            break;

        // safe screening
        screen(params, n_vars, iter);
        logger(Logger::INFO, logit(iter), " screening finished. #A: ", _active_indices.size());

        // narrow down to one reference
        if (n_vars > 1 && n_dynamic >= _multi_dynamic_screen)
        {
            logger(Logger::INFO, logit(0), " selected: ", index_gap_min + 1, '/', n_vars);
            std::swap(_variables[0], _variables[index_gap_min]);
            n_vars = 1;
        }
    }


    logger(Logger::INFO, logit(iter), ' ', index_gap_min+1, '/', n_vars,
        " #A: ", std::setw(LOGW), std::setw(LOGW), _active_indices.size(),
        " gap: ", _variables[index_gap_min]._relative_gap);


    // calculate computing time
    auto end = system_clock::now();
    _time = std::max(duration_cast<seconds>(end - start).count(), 0.0);

    // copy variables to solution
    solution = std::move(_variables[index_gap_min].dump(_activeset));
    return is_optimal;
}


// screen by mining algorithm
template< class Enumerator, class Loss, class XGenerator, class WGenerator >
void PMOpt::Optimizer<Enumerator, Loss, XGenerator, WGenerator>::screen_mining(
    const Parameters & params,
    unsigned int size,
    unsigned int iter
)
noexcept
{
    // safe pattern pruning
    // define pruner
    MultiSafePruner<Loss, XGenerator, WGenerator> pruner(_variables, params, size, logit(iter));

    // pruning and screening
    _activeset.clear();

    // for numerical feature (only screening)
    for (unsigned int j = 0; j < params._X_num.size(); ++j)
        pruner.evaluate_numeric(_activeset, j); 

    // for pattern feature
    _miner.mine(_activeset, pruner); 

    // post processing for keys
    for (auto pattern : _activeset._patterns)
        _activeset._keys.emplace_back(_miner._enumerator.pat_to_str(* pattern));   

    // initialize mined features
    // set reference _coef into variables _coef
    _active_indices.resize(_activeset._keys.size());

    // fill active indices
    std::iota(_active_indices.begin(), _active_indices.end(), 0);
}


// screen by linear search
template< class Enumerator, class Loss, class XGenerator, class WGenerator >
void PMOpt::Optimizer<Enumerator, Loss, XGenerator, WGenerator>::screen(
    const Parameters & params,
    unsigned int n_vars,
    unsigned int iter
)
noexcept
{
    // compute distance matrix and inclusion relationship
    std::vector<std::vector<double>> distance;
    std::vector<unsigned int> including;
    Optimization::check_inclusion<Loss>(distance, including, _variables, params, n_vars);

    Logger logger(std::cout);
    if (Logger::_verbose <= Logger::INFO)
    {
        for (unsigned int i_vars_1 = 0; i_vars_1 < n_vars; ++i_vars_1)
        {
            for (unsigned int i_vars_2 = i_vars_1 + 1; i_vars_2 < n_vars; ++i_vars_2)
            {
                const auto & vars_1 = _variables[i_vars_1];
                const auto & vars_2 = _variables[i_vars_2];
                logger(Logger::INFO, logit(iter), " INTERSECTION (", i_vars_1, ", ", i_vars_2, ") ",
                    " include: ", including[i_vars_1], ' ', including[i_vars_2],
                    " radius: ", vars_1._radius, ' ', vars_2._radius,
                    " distance: ", distance[i_vars_1][i_vars_2]);
            }
        }
    }

    // simple safe screening (not safe pattern pruning)
    auto n_active = _active_indices.size();
    for (unsigned int i_active = 0; i_active < n_active; ++i_active)
    {
        auto n = params._n_training;
        auto i_feature = _active_indices[i_active];
        auto x_sum_sq = _activeset._x_sum_sq[i_feature];
        auto x_sum = _activeset._x_sum[i_feature];
        bool is_screened = false;

        // single reference screening
        for (unsigned int i_vars = 0; i_vars < n_vars; ++i_vars)
        {
            // finish check screening if feature is screened by at least one variables
            if (is_screened)
                break;

            // skip inactive sphere
            if (including[i_vars])
                continue;
            
            const auto & vars = _variables[i_vars];

            // compute screenign criterion
            double criterion;
            Optimization::safe_screening(criterion,
                vars._scale, vars._radius, vars._score[i_feature], x_sum, x_sum_sq, n);

            // update flag that the feature is screened or not 
            is_screened = criterion < params._lambda_l1 * vars._reg_weight[i_feature];
        }

        // two references screening
        for (unsigned int i_vars_1 = 0; i_vars_1 < n_vars; ++i_vars_1)
        {
            // finish check screening if feature is screened by at least one variables
            if (is_screened)
                break;

            // skip inactive sphere
            if (including[i_vars_1])
                continue;

            for (unsigned int i_vars_2 = i_vars_1 + 1; i_vars_2 < n_vars; ++i_vars_2)
            {
                // finish check screening if feature is screened by at least one variables
                if (is_screened)
                    break;

                // skip inactive sphere
                if (including[i_vars_2])
                    continue;

                const auto & vars_1 = _variables[i_vars_1];
                const auto & vars_2 = _variables[i_vars_2];
                double d = distance[i_vars_1][i_vars_2];

                double criterion;
                Optimization::safe_screening_multi(criterion,
                    d, vars_1._scale, vars_2._scale, vars_1._radius, vars_2._radius,
                    vars_1._score[i_feature], vars_2._score[i_feature], x_sum, x_sum_sq, n);

                // update flag
                is_screened = criterion < params._lambda_l1 * _activeset._reg_weight[i_feature];
            }
        }

        if (!is_screened)
            continue;

        // remove inactive coef
        // coef is not actually deleted, only i_vars of coef is deleted
        const auto & x_i = * _activeset._x_indices[i_feature];
        const auto & x_v = * _activeset._x_values[i_feature];
        for (unsigned int i_vars = 0; i_vars < n_vars; ++ i_vars)
            for (unsigned int i_occ = 0; i_occ < x_i.size(); ++i_occ)
                _variables[i_vars]._y_pred[x_i[i_occ]] -= x_v[i_occ] * _variables[i_vars]._coef[i_feature]; // no masking

        // remove from indices
        --n_active;
        std::swap(_active_indices[i_active], _active_indices[n_active]);
        --i_active;
    }
    _active_indices.resize(n_active);

    // re-initailize coef indices of variables
    for (unsigned int i_vars = 0; i_vars < n_vars; ++ i_vars)
    {
        auto & vars = _variables[i_vars];
        vars._coef_indices.resize(n_active);
        std::copy(_active_indices.begin(), _active_indices.end(), vars._coef_indices.begin());
    }
}


// search max by mining algorithm
template< class Enumerator, class Loss, class XGenerator, class WGenerator >
double PMOpt::Optimizer<Enumerator, Loss, XGenerator, WGenerator>::search_max_mining(
    const Variables & vars,
    const Parameters & params
)
noexcept
{
    // define lammax calculator and result object
    LammaxCalculator<Loss, XGenerator, WGenerator> pruner(vars, params);
    LammaxObject<Symbol> obj;

    // mine lammax
    for (unsigned int i_feature = 0; i_feature < params._X_num.size(); ++i_feature)
        pruner.evaluate_numeric(obj, i_feature); // for numerical feature
    _miner.mine(obj, pruner); // for pattern feature
    return obj._max_score;
}


// search max by linear search
template< class Enumerator, class Loss, class XGenerator, class WGenerator >
double PMOpt::Optimizer<Enumerator, Loss, XGenerator, WGenerator>::search_max(
    Variables & vars,
    const Parameters & params
)
const noexcept
{
    double max_score = 0.0;
    for (auto i_feature : vars._coef_indices)
    {
        const auto & x_i = * _activeset._x_indices[i_feature];
        const auto & x_v = * _activeset._x_values[i_feature];
        Optimization::dot<Loss>(vars._score[i_feature], x_i, x_v, vars, params);
        max_score = std::max(max_score, std::abs(vars._score[i_feature]) / vars._reg_weight[i_feature]);
    }
    return max_score;
}


// compute dual norm by mining algorithm
template< class Enumerator, class Loss, class XGenerator, class WGenerator >
double PMOpt::Optimizer<Enumerator, Loss, XGenerator, WGenerator>::compute_dual_norm_mining(
    const Variables & vars,
    const Parameters & params
)
noexcept
{
    // define calculator and result object
    DualNormCalculator<Loss, XGenerator, WGenerator> pruner(vars, params);
    DualNormObject obj;

    // sum over all features
    for (unsigned int i_feature = 0; i_feature < params._X_num.size(); ++i_feature)
        pruner.evaluate_numeric(obj, i_feature); // for numerical feature
    _miner.mine(obj, pruner); // for pattern feature
    return obj._sum_score / 2 / params._lambda_l2;
}


// compute dual norm by lienar search
template< class Enumerator, class Loss, class XGenerator, class WGenerator >
double PMOpt::Optimizer<Enumerator, Loss, XGenerator, WGenerator>::compute_dual_norm(
    Variables & vars,
    const Parameters & params
)
const noexcept
{
    double norm_dual = 0.0;
    for (auto i_feature : vars._coef_indices)
    {
        const auto & x_i = * _activeset._x_indices[i_feature];
        const auto & x_v = * _activeset._x_values[i_feature];
        Optimization::dot<Loss>(vars._score[i_feature], x_i, x_v, vars, params);
        auto score = std::abs(vars._score[i_feature]) - params._lambda_l1 * vars._reg_weight[i_feature];
        if (score > 0)
            norm_dual -= score * score / vars._reg_weight[i_feature];
    }
    return norm_dual / 2 / params._lambda_l2;
}


// compute duality gap
template< class Enumerator, class Loss, class XGenerator, class WGenerator >
double PMOpt::Optimizer<Enumerator, Loss, XGenerator, WGenerator>::compute_gap(
    Variables & vars,
    const Parameters & params
)
const noexcept
{
    // norm of primal problem
    vars._norm_primal = vars._norm_l1 * params._lambda_l1 + vars._norm_l2 * params._lambda_l2 * 0.5;

    // primal objective function
    vars._obj_primal = vars._loss_primal + vars._norm_primal;

    // dual objective function
    vars._obj_dual = vars._loss_dual + vars._norm_dual;

    // duality gap and screening radius
    vars._gap = vars._obj_primal - vars._obj_dual;
    vars._relative_gap = vars._gap / vars._obj_primal;
    vars._radius = std::sqrt(2 * vars._gap * Loss::SMOOTHNESS * params._w_max);

    return vars._relative_gap;
}


template< class Enumerator, class Loss, class XGenerator, class WGenerator >
void PMOpt::Optimizer<Enumerator, Loss, XGenerator, WGenerator>::update(
    Variables & vars,
    const Parameters & params
)
noexcept
{
    // shuffle order of updating coef
    if (!_sequential_update)
        std::shuffle(vars._coef_indices.begin(), vars._coef_indices.end(), _random_generator);

    // for coefficient of model
    for (auto i_feature : vars._coef_indices)
    {
        // TODO check whether update all instance or training instance ?
        const auto & x_i = * _activeset._x_indices[i_feature];
        const auto & x_v = * _activeset._x_values[i_feature];

        // new coefficient
        auto coef = Loss::coefficient(x_i, x_v, vars, params, i_feature);

        // update dual variables
        for (unsigned int i_occ = 0; i_occ < x_i.size(); ++i_occ)
            vars._y_pred[x_i[i_occ]] += x_v[i_occ] * (coef - vars._coef[i_feature]); // no masking

        // update primal variable
        vars._coef[i_feature] = coef;
    }


    // new intercept
    auto intercept = Loss::intercept(vars, params);

    // update dual variables
    double optimality = 0.0;
    for (auto i_instance : * params._instance_indices)
    {
        vars._y_pred[i_instance] += intercept - vars._intercept; // no masking
        optimality += Loss::dual(vars, params, i_instance) * params._mask[i_instance];
    }
    assert(std::abs(optimality) < ZERO_TOL);

    // update primal variable
    vars._intercept = intercept;
}
