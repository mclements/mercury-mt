%% Title: test_mt.m
%% Description: Test for the Mercury library for Mersenne Twister random number generator
%% Licence: GPL-3
%% Author: Mark Clements
%% Date: 2023-02-01
%% Version: 0.1.1

:- module test_mt.

:- interface.
:- import_module io.
:- pred main(io::di, io::uo) is det.

:- implementation.

:- import_module mt, float, list, int, uint32, string, random.

:- pred runif(int::in, list(float)::out, random::in, random::out).
runif(N, List, !Rng) :- 
    map_foldl((pred(_I::in,U::out,R0::in,R::out) is det :- runif(U, R0, R)),
	      1..N, List, !Rng).

main(!IO) :-
    R0 = seed(12345u32),
    N = 624,
    runif(2*N,Us,R0,_),
    print_line("Expected: {0.72090390, 0.07196114}", !IO),
    list.det_index1(Us,1,U1),
    list.det_index1(Us,2*N,U2),
    format("Observed: {%10.8f, %10.8f}\n", [f(U1), f(U2)], !IO).

/*

R -q -e "set.seed(12345); runif(624*2)[c(1,624*2)]"
> [1] 0.72090390 0.07196114

R -q -e "set.seed(12345); microsimulation::unsigned(.Random.seed); runif(1); microsimulation::unsigned(.Random.seed)"

*/
