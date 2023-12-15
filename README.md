# Instructions

>This program is intended to generate the referrence answers to test the ongoing embedded rounding feature enabling project in C# JIT backend.

As a quick start, please compile the project with `Visual Studio 2022` under `Debug` build on a machine with `AVX512F` support. 

Run the executable and the answer table will be generated in the `output.txt`

Currently, this program would only generate the arithmetic intrisics with embedded rounding under `ToNegativeInfinity`, `ToPositiveInfinity` and `ToZero` with fixed inputs.

