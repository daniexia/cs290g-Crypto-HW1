cs290g-Crypto-HW1
=================

modular exponentiation

Apr 22, 2013

Four files are contained in the directory.
basic_impl.cpp: standard modular multiplication, binary method;
exp_opt.cpp:    standard modular multiplication, m-ary method;
mult_opt.cpp:   Blakley's shift-add method, binary method;
modexp.cpp: all the above three implementations.


1 Abstract
In this report, I will show an implementation of basic modular exponentiation, an implementation of exponentiation optimization using m-ary method, and an implementation of modular multiplication optimization using shift-add Blakley method.

The platform is the CSIL computer. It is a 4-core CPU running 64-bit Fedora Core 18 Linux. The language is C++. Below are the source files.
1) basic_impl.cpp - basic implementation (standard mult, modular operation, binary method);
2) exp_opt.cpp - m-ary method, standard multiplication and modular operation;
3) mult_opt.cpp - Blakley’s method (shift-add), binary method for exponentiation;
4) modexp.cpp - contains all the above three implementations. This file is for timing test. First, a random K-bit (256, 512, 1024, or 2048) big number is generated. Then the big number becomes the same input for all the three algorithms. The three algorithms run one after another to see the difference of the timing.

2 Implementation
1) Basic Implementation (source file: basic_impl.cpp)
Basic implementation uses standard multiplication, modular operation, and binary exponentiation method. Computation of the remainder uses modified Restoring Division Algorithm. The difference is to subtract n from Ri-1 only if Ri-1 >= n in Step 4 of Restoring Division Algorithm.

2) Exponentiation Optimization (source file: exp_opt.cpp)
In the exponentiation optimization, m-ary method is used. The function prototype is 
Bignum mod_exp_mary(const Bignum& exp, const Bignum& n);
Typically, m is 8. This value can be changed to 4, 16, 32 for 4-ary, 16-ary and 32-ary. To change this, modify the following line, and then re-compile it.
#define M_ARY     8

3) Modular Multiplication Optimization (source file: mult_opt.cpp)
In this optimization, instead of using standard modular multiplication, Blakley’s method (shift-add) is used. Shift-add method interleaves the multiplication and shift-subtract of division.
