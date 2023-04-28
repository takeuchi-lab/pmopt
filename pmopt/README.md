# PMOpt: Pattern Mining and Optimization tools

## Requirements

- `std++20`
    - `gcc`: `>=10`
    - `clang`: `>=12.0`
- `cmake`: `>=3.25.1`
- [Boost C++ Libraries](https://www.boost.org/)
    - Please download [Boost C++ Libraries](https://www.boost.org/), and copy all contents in the `boost` directory in the downloaded files to the directory `sysinclude/boost`.

## Compilation

You can compile by running the following command.
```
make install
```

To check if it is successfully compiled, Run the following commands.
```
./pmopt -v
```

## Usage

```
./pmopt SUBCOMMAND INPUT [options] ...
```

For detail options, please refer to the help.
```
./pmopt -h
```

### Subcommands

There are two subcommands in PMOpt

- `mine`: Mining patterns from input database. In this subcommand, you can specify the property of the patterns to be mine.
    ```
    ./pmopt mine PROPERTY INPUT [options]
    ```
    You can choose `PROPERTY` from the following keywords:
    - `frequent`
    - `preditive`

    For example, the following command enumerates frequent patterns and writes them to `patterns.csv`.
    ```
    ./pmopt mine frequent ./input.csv -o patterns.csv
    ```

- `predict`: Predict objective variables of new instances using preditive mining model obtained by `mine predictive`. For example, you need to o
    1. Obtain the predictive patterns from training dataset.
        ```
        ./pmopt mine predictive ./train.csv -o models.csv
        ```
    2. Predict objective variables of test instances.
        ```
        ./pmopt predict ./models.csv ./test.csv -o ./predictions.csv
        ```

### Inputs

We expect input file to be CSV. It must contain a structure column, where each entry represent structured input. The structure type (graphs, sequences, or sets) will be automatically determined by the name of the column. The name is chosen from the following keywords.
- `graph`
- `itemset`
- `sequence`
- `string`
