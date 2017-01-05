# CoIoTe Project Solution

This is the CoIoTe assignment solution proposed by group 34. The group members are:
 * Moriello Paolo
 * Coccia Giuseppe
 * Guerriero Mario
 * Marallo Graziano

A `Makefile` is provided along the code in order to make it easier to compile.

Before compiling the users must specify how many threads the program should be allowed to use. It can be done by simply changing the value of the `THREADS_NO` variable from the mentioned `Makefile`.

Before executing the user should specify which instance the program should solve. It can be done by simply changing the `INSTANCE` variable from the mentioned `Makefile`.

Once the `THREADS_NO` and the `INSTANCE` variables are properly set, the compile and execution steps can be done with the following two steps:
```
make
make run
```

The user can also check the feasibility of the proposed solution by mean of the following command:
```
make feasibility-check
```

All the instances include in the `input` directory could be solved at once by executing the `solveAll.sh` script.

The file `data/results.xlsx` contains a detailed description of the results obtained by resolving the instances given by the problem.
