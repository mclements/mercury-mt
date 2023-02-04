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

%% :- pred runif(int::in, list(float)::in, list(float)::out, R::di, R::uo) is det <= random(R).
%% runif(N, !List, !Record) :- runif_loop(1, N, [], List, !Record)
%%     %% map_foldl((pred(_::in,U::out,R0::di,R::uo) is det :- uniform_float_in_01(U, R0, R)),
%%     %% 	      1..N, List, !Record).

%% runif_loop(J, N, !List, !Record) :-
%%     if J=N then !:Record=!.Record, !:List=!.List
%%     else
%%     uniform_float_in_01(U,!Record),
%%     runif_loop(J+1, N, [U|List], !:List, !Record).

main(!IO) :-
    R0 = seed(12345u32),
    print_line(R0, !IO),
    uniform_float_in_01(U, R0, R1),
    print_line(U, !IO),
    print_line(R1, !IO).
	%% N = 624,
	%% runif(2*N,Us,R,_),
	%% print_line("Expected: {0.72090390, 0.07196114}", !IO).
	%% list.det_index0(Us,0,U1),
        %% list.det_index0(Us,2*N-1,U2),
	%% format("Observed: {%10.8f, %10.8f}\n", [f(U1), f(U2)], !IO).

/*

## *Not* comparable to the following code: random(R) uses 64 bit for random draws cf. 32 bit.
R -q -e "set.seed(12345); runif(624*2)[c(1,624*2)]"
> [1] 0.72090390 0.07196114

R -q -e "set.seed(12345); microsimulation::unsigned(.Random.seed); runif(1); microsimulation::unsigned(.Random.seed)"

*/
