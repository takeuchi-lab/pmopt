#pragma once

namespace PMOpt
{
    namespace Utils
    {

        // print help and exit
        inline void exit_with_help()
        {
            std::cout << 
                "PMOpt: Pattern Mining and Optimization tool\n"
                "Usage: ./pmoptz <subcommand> [options] ...\n"
                "\n"
                "<subcommand>:\n"
                "    mine           Mining patterns from input database.\n"
                "    predict        Predict objective variables of input database\n"
                "                   using pattern-based prediction model created by mine subcommand.\n"
                "    -v --version   Print the version of this program\n"
                "    -h --help      Print this message\n"
                "\n"
                "Usage of mine: ./pmoptz mine <mining_task> <database> [options] ...\n"
                "Usage of predict: ./pmoptz predict <model> <database> [options] ...\n"
                "\n"
                "<database>:        Input database from which patterns are mined.\n"
                "                   It must be csv format and have the input structure column.\n"
                "                   It also requires the column of objective variables or classes in some mining task.\n"
                "                   You can include columns of instance weight or other numerical features.\n"
                "\n"
                "<model>:           Model csv file which include pattern structure and its coefficient.\n"
                "                   Please use the file that the program retuerns in predictive mining.\n"
                "                   Structure of pattern must equal to the structure of input database.\n"
                "\n"
                "<mining_task>:\n"
                "    frequent       Mining patterns which frequently occurred in database.\n"
                "    predictive     Mining patterns which can be used to predict objective variables of database.\n"
                "    significant    Not implemented.\n"
                "\n"
                "[options]:\n"
                "    --maxpat       Maximum length of pattern which is mined by the program. (default 3)\n"
                "    --minsup       Frequency level of patterns in frequent mining\n"
                "                   If zero then all substructure are determined to be frequent.\n"
                "                   and if one then all substructure are not determined to be frequent (default 0.9)\n"
                "    --verbose      Verbosity level of standard output. (default 5)\n"
                "    --maxiter      Maximum number of iteration of the optimization problem in predictive mining (default 1000000)\n"
                "\n"
                "    --screening-freq       Frequency of safe screening in predictive mining (default 10)\n"
                "    --screening-freq-init  Frequency of safe screening in predictive mining, used only in first initial dynamic screening (default 2)\n"
                "    --tolerance            Convergence tolerance of the optimizatio problem in predictive mining (default 1e-4)\n"
                "    --significance-level   Significance level of statistical test used in significant mining (default 0.05)\n"
                "    -o --output-dir        Directory in which the output files will be stored. (default outputs)\n"
                "    --loss                 Type of loss function. choose from following key words (default squared)\n"
                "\n"
                "       squared         for regression problem\n"
                "       squaredhinge    for binary classification\n"
                "       logistic        for binary probability regression or classification\n"
                "\n"
                "    --maxtrial     Miximum number of trials in pathwise computation (default 10)\n"
                "    --l1-step      Log scale step size of l1-regularization parameter (default 0.1)\n"
                "    --l2-step      Step size of l2-regularization parameter (default 10.0)\n"
                "    --l1-min       Mimimum value of l1-regularization parameter (default 1e-2)\n"
                "    --l2-min       Mimimum value of l2-regularization parameter (default 1e-2)\n"
                "    --l2-max       Maximum value of l2-regularization parameter (default 1e+2)\n"
                "    --seed         Seed for random number generation (default 0)\n"
                "    --dry-run      Check pathwise compute schedule without actual optimization\n"
                "    --no-l2        Disable l2 regularization\n"
                "\n"
                "    --init-dynamic-screen      # of dynamic screening that screening will be more frequenty done (default 5)\n"
                "    --multi-dynamic-screen     # of dynamic screening with multiple references (default 0)\n"
                "    --disable-multi-screen     Disable multi reference safe screening (default enabled)\n"
                "    --enable-multi-screen      Enable multi reference safe screening (default enabled)\n"
                "    --disable-multi-update     Disable update of multi reference feasible solution (default disabled)\n"
                "    --enable-multi-update      Enable update of multi reference feasible solution (default disabled)\n"
                "    --disable-gap-select       Disable selection of the references of minimum duality gap (default enabled)\n"
                "    --enable-gap-select        Enable selection of the references of minimum duality gap (default enabled)\n"
                "    --select-model             Select best regularization parameter by K-fold cross validation.\n"
                "    --n-folds                  # of folds of cross validation in the model selection (default 5).\n"
                "    --max-folds                # of maximum folds of cross validation. If it is less than # of folds,\n"
                "                               then a fold over # of folds is omitted. (default 0, means infity).\n"
                "    --leave-one-out            Use leave-one-out cross validation if the model selection is enabled.\n"
                "    --leave-group-out          Use leave-group-out cross validation if the model selection is enabled,\n"
                "                               and `group` column is given in input csv file.\n"
                "    --sequential-update        If specified, update order randomization in coordinate descent will be disabled.\n"
                "    --shuffle-loocv-order      If specified, the order of trials in loocv will be shuffled.\n"
            ;

            exit(1);
        }

    }
}