#include <filesystem>
#include <iostream>
#include "pmopt/databases/graph_database.hpp"
#include "pmopt/databases/itemset_database.hpp"
#include "pmopt/databases/sequence_database.hpp"
#include "pmopt/enumerators/graph_enumerator.hpp"
#include "pmopt/enumerators/itemset_enumerator.hpp"
#include "pmopt/enumerators/sequence_enumerator.hpp"
#include "pmopt/losses/logistic_loss.hpp"
#include "pmopt/losses/squared_loss.hpp"
#include "pmopt/losses/squaredhinge_loss.hpp"
#include "pmopt/optimization/generator.hpp"
#include "pmopt/pruners/frequent_pruner.hpp"
#include "pmopt/pruners/safe_pruner.hpp"
#include "pmopt/utils/containers.hpp"
#include "pmopt/utils/cross_validation.hpp"
#include "pmopt/utils/help.hpp"
#include "pmopt/utils/misc.hpp"
#include "pmopt/utils/version.hpp"
#include "pmopt/utils/write.hpp"
#include "pmopt/loader.hpp"
#include "pmopt/logger.hpp"
#include "pmopt/miner.hpp"
#include "pmopt/optimizer.hpp"
#include "pmopt/parameters.hpp"
#include "pmopt/variables.hpp"

using namespace PMOpt;

// constant and/or static variables
const std::string FILE_FORMAT = ".csv";
const std::string SUBCMD_LOAD = "load";
const std::string SUBCMD_MINE = "mine";
const std::string SUBCMD_PREDICT = "predict";
const std::string MINING_FREQUENT = "frequent";
const std::string MINING_PREDICTIVE = "predictive";
const std::string MINING_SIGNIFICANT = "significant";
const std::string STRUCT_GRAPH = "graph";
const std::string STRUCT_ITEMSET = "itemset";
const std::string STRUCT_SEQUENCE = "sequence";
const std::string STRUCT_STRING = "string";
const std::string LOSS_LOGISTIC = "logistic";
const std::string LOSS_SQUARED = "squared";
const std::string LOSS_SQUAREDHINGE = "squaredhinge";
const unsigned int WIDTH_LAMBDA = 8;
const double LAMBDA_INFTY = 1e+16;
const double REG_L1_BASE = 10.0;

// static and global variables
unsigned int Logger::_verbose = 1;
unsigned int LogisticLoss::_maxiter_inner = 10;
unsigned int LogisticLoss::_maxiter_middle = 10000;
unsigned int SquaredHingeLoss::_maxiter_inner = 10;
unsigned int SequenceEnumerator::_max_skip = 0;
double LogisticLoss::_criterion_scale = 0.5;
double LogisticLoss::_linsearch_update = 0.5;
double SquaredHingeLoss::_criterion_scale = 0.5;
double SquaredHingeLoss::_linsearch_update = 0.5;
std::string subcommand;
std::string input;
std::string outputs = "outputs";
std::string mining_task;
std::string structure_type; 
std::string loss_type = LOSS_SQUARED;
unsigned int seed = 0;
unsigned int maxpat = 3;
unsigned int maxiter = 5;
unsigned int maxtrial = 10;
unsigned int max_folds = 0;
unsigned int init_dynamic_screen = 5;
unsigned int multi_dynamic_screen = 0;
unsigned int screening_freq = 10;
unsigned int screening_freq_init = 2;
unsigned int n_folds = 5;
double reg_l1_step = 0.1;
double reg_l1_min = 1e-2;
double reg_l2_step = 10.0;
double reg_l2_max = 1e+2;
double reg_l2_min = 1e-2;
double significance_level = 0.05;
double tolerance = 1e-4;
double minsup = 0.9;
bool dry_run = false;
bool nol2 = false;
bool allow_gap_select = true;
bool allow_multi_screen = true;
bool allow_multi_update = false;
bool select_model = false;
bool leave_one_out = false;
bool leave_group_out = false;
bool do_l2_scale = false;
bool sequential_update = true;
bool shuffle_loocv_order = false;


// parse command-line arguments
void parse_args(int, char **);


// execute specified task
template < class Enumerator, class Loss, class XGenerator, class WGenerator >
void run();


// compuate regularization path
template < class Enumerator, class Loss, class XGenerator, class WGenerator >
void compute_path(Miner<Enumerator> &, Parameters &);


int main(int argc, char ** argv)
{
    // parse command-line arguments and initialize global variables
    parse_args(argc, argv);

    // create directory
    std::filesystem::create_directory(outputs);

    // run program for specified structure
    if (structure_type == STRUCT_ITEMSET)
    {
        if (loss_type == LOSS_SQUARED)
            run<ItemsetEnumerator, SquaredLoss, BinaryXGenerator, LengthWGenerator>();
        else if (loss_type == LOSS_SQUAREDHINGE)
            run<ItemsetEnumerator, SquaredHingeLoss, BinaryXGenerator, LengthWGenerator>();
        else if (loss_type == LOSS_LOGISTIC)
            run<ItemsetEnumerator, LogisticLoss, BinaryXGenerator, LengthWGenerator>();
    }
    else if (structure_type == STRUCT_SEQUENCE || structure_type == STRUCT_STRING)
    {
        if (loss_type == LOSS_SQUARED)
            run<SequenceEnumerator, SquaredLoss, BinaryXGenerator, LengthWGenerator>();
        else if (loss_type == LOSS_SQUAREDHINGE)
            run<SequenceEnumerator, SquaredHingeLoss, BinaryXGenerator, LengthWGenerator>();
        else if (loss_type == LOSS_LOGISTIC)
            run<SequenceEnumerator, LogisticLoss, BinaryXGenerator, LengthWGenerator>();
    }
    else if (structure_type == STRUCT_GRAPH)
    {
        if (loss_type == LOSS_SQUARED)
            run<GraphEnumerator, SquaredLoss, BinaryXGenerator, LengthWGenerator>();
        else if (loss_type == LOSS_SQUAREDHINGE)
            run<GraphEnumerator, SquaredHingeLoss, BinaryXGenerator, LengthWGenerator>();
        else if (loss_type == LOSS_LOGISTIC)
            run<GraphEnumerator, LogisticLoss, BinaryXGenerator, LengthWGenerator>();
    }
}


// main function
template < class Enumerator, class Loss, class XGenerator, class WGenerator >
void run()
{
    using Database = typename Enumerator::Database;
    using Symbol = typename Enumerator::Symbol;

    // log into standard output
    Logger logger(std::cout);
    logger(10, "Input: ", input);

    // load database
    Database database(structure_type);
    Parameters params;
    params._l2_scale = do_l2_scale;
    Loader::load(database, params, input);

    logger(Logger::INFO, "Database loaded\n -- size: ", database.size(), " instances");
    logger(Logger::DEBUG, " -- y: ", params._y.empty() ? "no" : "yes");
    logger(Logger::DEBUG, " -- w: ", params._w.empty() ? "no" : "yes");
    logger(Logger::DEBUG, " -- ", params._X_num.size(), " numerical features");

    // initialize parameters
    params.initialize(database.size(), Loss::NUMERIC_OBJECTIVE);

    if (leave_group_out && params._group.size() == params._y.size())
    {
        std::cerr << "Error: Missing group information in leave-group-out cv\n";
        exit(1);
    }

    if (subcommand == SUBCMD_LOAD)
    {
        for (const auto & data : database)
            logger(Logger::DEBUG, data);
        return;
    }

    if (subcommand == SUBCMD_MINE)
    {
        // initialize miner
        Enumerator enumerator(database);
        Miner miner(enumerator, maxpat);

        // frequent mining
        if (mining_task == MINING_FREQUENT)
        {
            // mine frequent patterns
            FrequentPruner pruner(std::ceil(minsup * database.size()));
            FrequentObject<Symbol> mined;
            miner.mine(mined, pruner);

            // save mined pattern
            auto path = std::filesystem::path(outputs).append("patterns.csv");
            auto time = miner._time;
            auto n_patterns = mined._patterns.size();
            Utils::write_frequent(mined, enumerator, path);
            logger(Logger::INFO, "Finished mining: ", n_patterns, " patterns, ", time, " [sec]");
            logger(Logger::INFO, "Saved: ", path);

            return;
        }

        // predictive mining
        if (mining_task == MINING_PREDICTIVE)
        {

            if (params._y.empty())
            {
                std::cerr << "Erorr: missing objective variables in predictive mining.\n";
                exit(1);
            }

            compute_path<Enumerator, Loss, XGenerator, WGenerator>(miner, params);
        }
    }
}


template < class Enumerator, class Loss, class XGenerator, class WGenerator >
void compute_path(Miner<Enumerator> & miner, Parameters & params)
{
    // l2_scale option is only for testing
    // do not use this option !!!
    assert(!params._l2_scale);

    Logger logger(std::cout);

    // initialize optimizer
    Optimizer<Enumerator, Loss, XGenerator, WGenerator> optimizer(miner, seed);
    optimizer._maxiter = maxiter;
    optimizer._init_dynamic_screen = init_dynamic_screen;
    optimizer._multi_dynamic_screen = multi_dynamic_screen;
    optimizer._screening_freq = screening_freq;
    optimizer._screening_freq_init = screening_freq_init;
    optimizer._tolerance = tolerance;
    optimizer._allow_multi_screen = allow_multi_screen;
    optimizer._allow_multi_update = allow_multi_update;
    optimizer._sequential_update = sequential_update;

    // initialize fold indices
    std::vector<unsigned int> fold_ids;
    if (!select_model)
    {
        n_folds = 1;
        fold_ids.clear();
        fold_ids.resize(params._y.size(), 1);
    }
    else
    {
        if (leave_one_out)
            n_folds = Utils::make_leave_one_out(fold_ids, params._y.size(), seed, shuffle_loocv_order);
        else if (leave_group_out)
            n_folds = Utils::make_leave_group_out(fold_ids, params._group);
        else
            n_folds = Utils::make_k_fold(fold_ids, params._y.size(), n_folds, seed);
    }


    // if max folds set, only some parts of fold will be optimized
    if (max_folds > 0)
        n_folds = std::min(n_folds, max_folds);


    // check number of lambda l2
    // n_l2 = max_k { reg_l2_min * reg_l2_step ** k <= reg_l2_max } + 1
    unsigned int n_l2 = 1;
    if (!nol2)
        n_l2 += std::floor(std::log(reg_l2_max / reg_l2_min) / std::log(reg_l2_step) + 1);

    // initialize solution history
    std::vector<std::vector<std::vector<Reference>>> history;
    history.resize(
        n_folds,
        std::vector<std::vector<Reference>>(
            n_l2,
            std::vector<Reference>(2) // no need to store all references of lambda_l1
        )
    );

    // initialize ROOT solution
    // this is used in optimization at lambda_l1 is lambda_max
    Reference root_reference;
    root_reference._y_pred.resize(params._y.size());

    // pathwise computation
    std::vector<unsigned int> test_indices;
    std::vector<double> lammax_list(n_folds, LAMBDA_INFTY);
    unsigned int trial = 0;
    unsigned int i_l1 = 0;
    double reg_l1 = 1.0;
    double total_time = 0.0;
    while (reg_l1 >= reg_l1_min) // for each l1 lambda
    {
        unsigned int i_l2 = 0;
        double reg_l2 = 0.0;
        while (i_l2 < n_l2) // for each l2 lambda
        {
            for (unsigned int i_fold = 0; i_fold < n_folds; ++i_fold) // for each fold of cv
            {
                // if lambda_l1 is lammax then lambda_l2 must be zero
                assert(i_l1 != 0 || i_l2 == 0);

                // split train and test indices
                Utils::make_test_mask(params._mask, params._n_training, fold_ids, i_fold);

                // set regularization parameters
                params._lambda_l1 = reg_l1 * lammax_list[i_fold];
                params._lambda_l2 = reg_l2 * params._lambda_l1;


                auto log_prefix = " [Iter:" + std::string(std::log10(maxiter) + 1, ' ') + "-] ";
                logger(Logger::INFO, 
                    "\n<< Trial: ", std::setw(std::log10(maxtrial)), trial,
                    " Fold: ", std::setw(std::log10(n_folds)), i_fold,
                    " L1: ", std::setw(WIDTH_LAMBDA), reg_l1,
                    " L2: ", std::setw(WIDTH_LAMBDA), reg_l2, " >>");
                logger(Logger::INFO, log_prefix, "lambda l1: ", params._lambda_l1, " lambda l2: ", params._lambda_l2);


                // if the previous lambda_l1 is not lambda_max
                //  then the solution of previous lambda_l1 exists,
                // even if previous lambda_l1 is lambda_max,
                //  the solution of previous lambda_l1 exists if lambda_l2 is zero.
                // i_l1 % 2 is used instead of i_l1
                //  because it is only necessary to store current and previous lambda_l1 reference.
                std::vector<const Reference *> references;
                if (i_l1 > 1 || (i_l1 == 1 && i_l2 == 0))
                {
                    // there exists a reference of previous lambda_l1
                    logger(Logger::INFO, log_prefix, "load prev. l1   (", i_fold, ',', i_l2, ',', i_l1 - 1, ')');
                    references.push_back(& history[i_fold][i_l2][(i_l1 - 1) % 2]);
                }

                // if lambda_l2 is not zero, then the solution of previous lambda_l2 exists.
                // if gap selection not allowed,
                //  the reference will be loaded only if no references have been loaded yet.
                if (i_l2 > 0 && (allow_gap_select || references.empty()))
                {
                    // there exists a reference of previous lambda_l2
                    logger(Logger::INFO, log_prefix, "load prev. l2   (", i_fold, ',', i_l2 - 1, ',', i_l1, ')');
                    references.push_back(& history[i_fold][i_l2 - 1][i_l1 % 2]);
                }

                // if fold is not zero, i.e. not optimization for entire dataset, 
                //  and lambda_l1 is not lambda_max then the solution of the entire dataest exists. 
                // if gap selection not allowed,
                //  the reference will be loaded only if no references have been loaded yet.
                if (i_fold > 0 && i_l1 > 0 && (allow_gap_select || references.empty()))
                {
                    // there exists a reference of privious dataaset
                    logger(Logger::INFO, log_prefix, "load prev. fold (", 0, ',', i_l2, ',', i_l1, ')');
                    references.push_back(& history[0][i_l2][i_l1 % 2]);
                }

                // if no solutions of the previous optimization
                //  then the ROOT is added
                // this occurs only if i_l1 == 0
                if (references.empty())
                {
                    logger(Logger::INFO, log_prefix, "load prev. ROOT");
                    references.push_back(& root_reference);
                }

                // if --dry-run specified in commandline then only show the path
                //  and not calculated actually
                if (dry_run)
                    continue;
   
                // standardize parameters
                params.standardize(Loss::NUMERIC_OBJECTIVE);

                // solve optimization and save as history
                auto & solution = history[i_fold][i_l2][i_l1 % 2];
                auto is_optimal = optimizer.optimize(solution, references, params);
                total_time += optimizer._time;

                // If at lambda_l1 = lambda_max then lammbda_max will be set.
                if (i_l1 == 0)
                {
                    // TODO : chcek validity of lambda_max
                    // Here, we set a value that is slightly larger than calculated lambda_max
                    // in order to ensure that all features are shrinked to be zero at lambda_l1 = lambda_max.
                    // lammax_list[i_fold] = solution._max_score * std::pow(REG_L1_BASE, reg_l1_step / 2);
                    lammax_list[i_fold] = solution._max_score;
                    params._lambda_l1 = solution._max_score;
                }

                // log covergence result
                if (is_optimal)
                    logger(Logger::INFO, " Successfully converged.");
                else
                    std::cerr << "Warning: Failed to converge in trial: "
                        << trial << " of fold: " << i_fold <<  ".\n";

                logger(Logger::INFO, " Time: ", optimizer._time, " [sec]");

                // save optimization result
                auto filename = "models_" +  std::to_string(i_fold) + ".csv";
                auto path = std::filesystem::path(outputs).append(filename);
                Utils::write_model(history[i_fold][i_l2][i_l1 % 2], params, path, trial, nol2);

                logger(Logger::INFO, " Saved: ", path);
            }

            ++trial;
            if (trial == maxtrial)
                break;

            // if l1 is lammax then there are no model parameters
            if (i_l1 == 0)
                break;

            // next l2 lambda
            if (i_l2 == 0)
                reg_l2 = reg_l2_min;
            else
                reg_l2 *= reg_l2_step;
            ++i_l2;
        }

        if (trial == maxtrial)
            break;

        // next l1 lambda
        ++i_l1;
        reg_l1 = std::pow(REG_L1_BASE, - (i_l1 * reg_l1_step));
    }

    logger(Logger::INFO, "Total Time: ", total_time, " [sec]");
}




// parse command-line arguments
void parse_args(int argc, char ** _argv)
{

    std::vector<std::string> argv((std::size_t) argc);
    unsigned int i = 0;
    for (i = 0; i < argv.size(); ++i)
        argv[i] = _argv[i];

    i = 1;

    // Set subcommand
    if (i >= (unsigned int) argc)
        Utils::exit_with_help();
    subcommand = argv[i++];


    if (subcommand == "--version" || subcommand == "-v")
        Utils::exit_with_version();
    if (subcommand == "--help" || subcommand == "-h")
        Utils::exit_with_help();


    // Check subcommand is valid
    if (!Utils::contains(subcommand, {SUBCMD_LOAD, SUBCMD_MINE, SUBCMD_PREDICT}))
    {
        std::cerr << "Error: Undefined subcommand `" << subcommand << "`.\n";
        exit(1);
    }


    // Set mining task name if subcommand is `mine`
    if (subcommand == SUBCMD_MINE)
    {   
        if (i >= (unsigned int) argc)
        {
            std::cerr << "Error: Missing mining task in commandline.\n";
            exit(1);
        }

        mining_task = argv[i++];
        if (!Utils::contains(mining_task, {
            MINING_FREQUENT,
            MINING_PREDICTIVE,
            MINING_SIGNIFICANT
        }))
        {
            std::cerr << "Error: Undefined mining task `" << mining_task << "`.\n";
            exit(1);
        }
    }

    // Set other options and input path
    for (; i < (unsigned int) argc; ++i)
    {
        std::string key(argv[i]);

        // Set input path
        if (key[0] != '-')
        {
            if (input != std::string())
            {
                std::cerr << "Error: Multiple input files given.\n";
                std::cerr << "\t`" << input << "` is already given\n\t`" << key << "` is newly given.\n";
                exit(1);
            }

            input = key;

            auto extension = std::filesystem::path(input).extension();
            if (input.size() < FILE_FORMAT.size() || extension != FILE_FORMAT)
            {
                std::cerr << "Error: Input file must be .csv.\n";
                exit(1);
            }

            structure_type = Utils::get_structure_type({
                STRUCT_GRAPH,
                STRUCT_ITEMSET,
                STRUCT_SEQUENCE,
                STRUCT_STRING
            }, input);

            if(!std::ifstream(input))
            {
                std::cerr << "Error: File cannot open `" << input << "`.\n";
                exit(1);
            }

            continue;
        }

        // Set other options
        std::string key_name(key.substr(2, key.size()));
        std::string key_type;
        bool is_undefined_key = false;
        try
        {
            if (key == "--maxpat")
            {
                key_type = "positive integer (must be >0)";
                maxpat = Utils::stopi(argv.at(++i));
            }
            else if (key == "--maxiter")
            {
                maxiter = Utils::stopi(argv.at(++i));
            }
            else if (key == "--minsup")
                minsup = std::stod(argv.at(++i));
            else if (key == "--tol")
                tolerance = std::stod(argv.at(++i));
            else if (key == "--significance-level")
                significance_level = std::stod(argv.at(++i));
            else if (key == "--output-dir" || key == "-o")
                outputs = argv.at(++i);
            else if (key == "--verbose")
                Logger::_verbose = Utils::stoui(argv.at(++i));
            else if (key == "--loss")
            {
                loss_type = argv.at(++i);
                if (!Utils::contains(loss_type, {
                    LOSS_SQUARED,
                    LOSS_SQUAREDHINGE,
                    LOSS_LOGISTIC
                }))
                {
                    std::cerr << "Error: Undefined loss `" << loss_type << "`\n";
                    exit(1);
                }
            }
            else if (key == "--maxtrial")
                maxtrial = Utils::stopi(argv.at(++i));
            else if (key == "--init-dynamic-screen")
                init_dynamic_screen = Utils::stoui(argv.at(++i));
            else if (key == "--multi-dynamic-screen")
                multi_dynamic_screen = Utils::stoui(argv.at(++i));
            else if (key == "--l1-step")
                reg_l1_step = std::stod(argv.at(++i));
            else if (key == "--l2-step")
                reg_l2_step = std::stod(argv.at(++i));
            else if (key == "--l1-min")
                reg_l1_min = std::stod(argv.at(++i));
            else if (key == "--l2-min")
                reg_l2_min = std::stod(argv.at(++i));
            else if (key == "--l2-max")
                reg_l2_max = std::stod(argv.at(++i));
            else if (key == "--screening-freq")
                screening_freq = Utils::stoui(argv.at(++i));
            else if (key == "--screening-freq-init")
                screening_freq_init = Utils::stoui(argv.at(++i));
            else if (key == "--seed")
                seed = Utils::stoui(argv.at(++i));
            else if (key == "--dry-run")
                dry_run = true;
            else if (key == "--no-l2")
                nol2 = true;
            else if (key == "--disable-multi-screen")
                allow_multi_screen = false;
            else if (key == "--enable-multi-screen")
                allow_multi_screen = true;
            else if (key == "--disable-multi-update")
                allow_multi_update = false;
            else if (key == "--enable-multi-update")
                allow_multi_update = true;
            else if (key == "--disable-gap-select")
                allow_gap_select = false;
            else if (key == "--enable-gap-select")
                allow_gap_select = true;
            else if (key == "--select-model")
                select_model = true;
            else if (key == "--n-folds")
                n_folds = Utils::stopi(argv.at(++i));
            else if (key == "--leave-one-out")
                leave_one_out = true;
            else if (key == "--leave-group-out")
                leave_group_out = true;
            else if (key == "--sequential-update")
                sequential_update = true;
            else if (key == "--max-folds")
                max_folds = Utils::stopi(argv.at(++i));
            else if (key == "--shuffle-loocv-order")
                shuffle_loocv_order = true;
            else
                is_undefined_key = true;
        }
        catch (const std::out_of_range &)
        {
            // If there no value after argument
            std::cerr << "Error: Missing value of commandline argmument `" << key_name << "`.\n";
            exit(1);
        }
        catch (const std::invalid_argument &)
        {
            // If there is type mismatch
            std::cerr << "Error: Type of commandline argument `" << key_name
                << "` must be " << key_type << ".\n";
            exit(1);
        }

        if (is_undefined_key)
        {
            std::cerr << "Error: Undefined commandline option `" << key_name << "`.\n";
            exit(1);
        }

        if (leave_one_out && leave_group_out)
        {
            std::cerr << "Error: both --leave-one-out and --leave-group-out are specified.\n";
            exit(1);
        }
    }

    if (input == std::string())
    {
        std::cerr << "Error: Missing input database file.\n";
        exit(1);
    }

    if (!allow_gap_select && allow_multi_update)
    {
        std::cerr << "Warning: Both multi-screening and multi-update will be disabled"
            << " if duality gap based selection is disabled.\n";
        allow_multi_screen = false;
        allow_multi_update = false;
    }

    if (!allow_multi_screen && allow_multi_update)
    {
        std::cerr << "Warning: Multi-screening will be disabled"
            << " if multi-update is disabled.\n";
        allow_multi_update = false;
    }

    if (multi_dynamic_screen == 0 && allow_multi_update)
    {
        std::cerr << "Warning: Multi-update will be disabled"
            << " if # of times of multi dynamic screening is zero.\n";
        allow_multi_update = false;
    }
}