# Optimizing B+ tree Node Layout for Improved Search Performance Using SIMD Instructions

Guided research project conducted at the [Chair for Database Systems](https://db.in.tum.de) at the Technical University of Munich, under supervision of [Altan Birler](https://db.in.tum.de/people/sites/birler) and Prof. [Thomas Neumann](https://db.in.tum.de/people/sites/neumann).

## Abstract

B<sup>+</sup>-tree is a ubiquotous data structure in database systems where it is used to index large datasets which do not fit into main memory. Its initial design was based on two premises: first, external memory has several orders of magnitude higher latency than main memory; second, unit of transfer between external memory and main memory is often much bigger than a single data record. In the following work, we aim to further improve the performance of B<sup>+</sup>-trees also by drawing inspiration from hardware characteristics– specifically the latency discrepancy between CPU caches and main memory, as well as availability of Single-Instruction-Multiple-Data (SIMD) instructions on modern processors. We experiment with different ways of organizing data within a single B<sup>+</sup>-tree node, conceived to more effectively leverage the aforementioned capabilities, with primary goal of optimizing lookups without worsening other operations like inserts, erasures and range-scans.
Our benchmarks show that lookups become 75% faster on average when we use Eytzinger layout combined with SIMD instructions, whereas inserts and erasures are at most 15% slower and often faster than with standard layout.

For a complete description of the project, take a look at the full [report](/report.pdf).

## Running the Code on Linux
Follow these steps to build and run the code on Linux:
```shell
# create build directory
$ mkdir build && cd build
# generate build files
$ cmake -DCMAKE_BUILD_TYPE=Release ..
# build project
$ make
# run tests
$ ./tester
# run benchmarks
$ ./bench
```
Make sure you have CMake and a C++ compiler installed on your system before running these commands.