# P4EST Search Partition

This is a C++ interface to test integration of p4est functionality 
into the Deal.II library.

The lates version of p4est can be found [here](https://github.com/cburstedde/p4est).

In particular function in these files are relevant: 
* [p4est_search_partition.h](https://github.com/cburstedde/p4est/blob/master/src/p4est_search.h)
* [p8est_search_partition.h](https://github.com/cburstedde/p4est/blob/master/src/p8est_search.h)

Mathematical documentation can be found in the following papers:
* ["Parallel tree algorithms for AMR and non-standard data access"](https://arxiv.org/pdf/1803.08432.pdf), 2018
* [Original p4est paper](https://ins.uni-bonn.de/media/public/publication-media/BursteddeWilcoxGhattas11.pdf),  2011
* ["Recursive Algorithms for Distributed Forests of Octrees"](https://arxiv.org/pdf/1406.0089.pdf), 2015

### How to compile?

Requirements:

* **cmake** v2.8.12 or higher	
*  **[Deal.II](www.dealii.org)** v9.3.0 or linked against **MPI** and [**p4est**](https://www.p4est.org/).
* **doxygen**	(optional)
* **clang-format-6.0** (optional)

	

We can switch between build modes:
```
make debug
```
switches the compile mode to debug mode (but does not compile anything) and

```
make release
```
switches to optimized mode. 


### How to extend the code (conventions)?

For small projects I usually adhere to this:

* header (`.h` files) are placed in the `include/` folder
* implementations (`.cpp` files) go into the `/source` folder
* if headers include templates they should be separated
  from the implementation as well in the following way:
  - headers with templates go into the `include/` folder
    as well as the implementation (`.tpp` files)
  - each header contains a list of instantiations, 
    example: `extern template class PartitionSearch<2>;`. This 
    tells the compiler that this class is instantiated somwhere else
  - a *.inst.cpp* file in the `source/` folder containing the
    actual instantiation, example: `template class PartitionSearch<2>;`.
    This file needs to include both the `.h` and the `.tpp` file
* executables should have the ending `.cxx` and be located in 
  the `source/` folder
* do not forget to add sources to the `source/CMakeLists.txt`
* a test would be nice


### Unit Tests

The Project comes with a unit test suite based on Deal.II's test suite.
Each test is a `.cpp` file in the test folder under the casis_code directory.
For how to set up a test see the documentation [here](https://www.dealii.org/current/users/testsuite.html).
A specific test can be run through

```
ctest -V -R "test/name_of_test"
```

or similarly through

```
ctest -V -R
```


### Code Indentation

The code should be indented. We use the style that is also used by Deal.II
to make code review easier. The code is automatically indented through

```
make indent
```


### Documentation

A documantation that can sometimes make code easier to understand
can be built using

```
make indent
```