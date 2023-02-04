%% Title: mercury-mt
%% Description: Mercury library for Mersenne Twister random number generator
%% Licence: GPL-3
%% Author: Mark Clements
%% Date: 2023-02-01
%% Version: 0.1.1

:- module mt_unique.

:- interface.
:- import_module random, uint8, uint16, uint32, uint64.

%% Types
:- type params.
:- type ustate.

:- instance urandom(params, ustate).
:- instance urandom_dup(ustate).

%% seed(Seed::in, State::out), constructs an instance based on R's set.seed algorithm
:- pred seed(uint32::in, params::out, ustate::uo) is det.

    % Generate a uniformly distributed pseudo-random unsigned integer
    % of 8, 16, 32 or 64 bits, respectively.
    %
:- pred generate_uint8(params::in, uint8::out,
    ustate::di, ustate::uo) is det.
:- pred generate_uint16(params::in, uint16::out,
    ustate::di, ustate::uo) is det.
:- pred generate_uint32(params::in, uint32::out,
    ustate::di, ustate::uo) is det.
:- pred generate_uint64(params::in, uint64::out,
    ustate::di, ustate::uo) is det.

    % Duplicate a 32-bit SFC state.
    %
:- pred urandom_dup(ustate::di, ustate::uo, ustate::uo) is det.

    % Generate a uniformly distributed pseudo-random unsigned integer
    % of 8, 16, 32 or 64 bits, respectively.
    %
    % As above, but does not require the params argument (which is a dummy
    % type only needed to satisfy the typeclass interface).
    %
:- pred generate_uint8(uint8::out, ustate::di, ustate::uo) is det.
:- pred generate_uint16(uint16::out, ustate::di, ustate::uo) is det.
:- pred generate_uint32(uint32::out, ustate::di, ustate::uo) is det.
:- pred generate_uint64(uint64::out, ustate::di, ustate::uo) is det.

:- pred runif(params::in, float::out, ustate::di, ustate::uo) is det.


:- implementation.

:- import_module array, list, int, float.

%% ===================  Mersenne Twister ==========================
%% Based on http://www.math.keio.ac.jp/~matumoto/emt.html

%%    A C-program for MT19937: Real number version([0,1)-interval)
%%    (1999/10/28)
%%    Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura.
%%    REFERENCE
%%    M. Matsumoto and T. Nishimura,
%%    \"Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
%%    Pseudo-Random Number Generator\",
%%    ACM Transactions on Modeling and Computer Simulation,
%%    Vol. 8, No. 1, January 1998, pp 3--30.

%% Types
:- type params
    --->    params.

:- type ustate
   --->    ustate(int, array(uint32)).

:- instance urandom(params, ustate) where [
    pred(generate_uint8/4) is mt_unique.generate_uint8,
    pred(generate_uint16/4) is mt_unique.generate_uint16,
    pred(generate_uint32/4) is mt_unique.generate_uint32,
    pred(generate_uint64/4) is mt_unique.generate_uint64
].

:- instance urandom_dup(ustate) where [
    pred(urandom_dup/3) is mt_unique.urandom_dup
   ].

/* Period parameters */
:- func n = int.
:- func m = int.
:- func matrix_a = uint32.
:- func upper_mask = uint32.
:- func lower_mask = uint32.
:- func mag01 = array(uint32).
n = 624.
m = 397.
matrix_a = 0x9908b0dfu32.   /* constant vector a */
upper_mask = 0x80000000u32. /* most significant w-r bits */
lower_mask = 0x7fffffffu32. /* least significant r bits */
mag01 = from_list([0x0u32, matrix_a]).

/* Tempering parameters */
:- func tempering_mask_b = uint32.
:- func tempering_mask_c = uint32.
:- func tempering_shift_u(uint32) = uint32.
:- func tempering_shift_s(uint32) = uint32.
:- func tempering_shift_t(uint32) = uint32.
:- func tempering_shift_l(uint32) = uint32.
tempering_mask_b = 0x9d2c5680u32.
tempering_mask_c = 0xefc60000u32.
tempering_shift_u(Y) = (Y >> 11).
tempering_shift_s(Y) = (Y << 7).
tempering_shift_t(Y) = (Y << 15).
tempering_shift_l(Y) = (Y >> 18).

:- pred mt_loop1(int::in, uint32::in, uint32::out)  is det.
mt_loop1(J, !Seed) :-
    if J=51 then !:Seed = !.Seed else mt_loop1(J+1, !.Seed*69069u32 + 1u32, !:Seed).
:- pred mt_loop2(int::in, uint32::in, uint32::out, array(uint32)::in, array(uint32)::out)  is det.
mt_loop2(J, !Seed, !Mt) :-
    if J=n then !:Seed = !.Seed, !:Mt = !.Mt
      else
      !:Seed = !.Seed*69069u32 + 1u32,
      Mt0 = !.Mt,
      !:Mt = Mt0 ^ elem(J) := !.Seed,
      mt_loop2(J+1, !Seed, !Mt).
seed(SeedIn, params, S) :-
    some ([!Mt,!Seed]) (
	!:Seed = SeedIn,
	!:Mt = array.init(n, 0u32),
	mt_loop1(0, !Seed),
	mt_loop2(0, !.Seed, _, !Mt),
	S = unsafe_promise_unique(ustate(n,!.Mt))).

generate_uint8(_, N, !S) :-
    mt_unique.generate_uint8(N, !S).

generate_uint16(_, N, !S) :-
    mt_unique.generate_uint16(N, !S).

generate_uint32(_, N, !S) :-
    mt_unique.generate_uint32(N, !S).

generate_uint64(_, N, !S) :-
    mt_unique.generate_uint64(N, !S).

urandom_dup(S, S1, S2) :-
    S = ustate(N,A),
    Sc = ustate(N,array.copy(A)),
    S1 = unsafe_promise_unique(S),
    S2 = unsafe_promise_unique(Sc).

generate_uint8(N, !S) :-
    mt_unique.generate_uint32(N0, !S),
    N1 = uint32.cast_to_int(N0 >> 24),
    N = uint8.cast_from_int(N1).

generate_uint16(N, !S) :-
    mt_unique.generate_uint32(N0, !S),
    N1 = uint32.cast_to_int(N0 >> 16),
    N = uint16.cast_from_int(N1).

:- pred mt_unif_rand_loop1(int::in, uint32::in, uint32::out,
			   array(uint32)::in, array(uint32)::out) is det.
mt_unif_rand_loop1(K,!Y,!Mt) :-
    (if K=n-m
    then !:Y = !.Y, !:Mt = !.Mt
    else
      !:Y = (!.Mt ^ elem(K) /\ upper_mask) \/ (!.Mt ^ elem(K+1) /\ lower_mask),
      Mt0 = !.Mt,
      !:Mt = Mt0 ^ elem(K) := xor(xor(Mt0 ^ elem(K+m), !.Y >> 1),
				  mag01 ^ elem(cast_to_int(!.Y /\ 0x1u32))),
      mt_unif_rand_loop1(K+1,!Y,!Mt)).

:- pred mt_unif_rand_loop2(int::in, uint32::in, uint32::out,
			   array(uint32)::in, array(uint32)::out) is det.
mt_unif_rand_loop2(K,!Y,!Mt) :-
    (if K=n-1
    then !:Y = !.Y, !:Mt = !.Mt
    else
    !:Y = (!.Mt ^ elem(K) /\ upper_mask) \/ (!.Mt ^elem(K+1) /\ lower_mask),
    Mt0 = !.Mt,
    !:Mt = Mt0 ^ elem(K) := xor(xor(Mt0 ^ elem(K+(m-n)), !.Y >> 1),
				mag01 ^ elem(cast_to_int(!.Y /\ 0x1u32))),
      mt_unif_rand_loop2(K+1,!Y,!Mt)).

:- pred xor(uint32::in,uint32::in,uint32::out) is det.
xor(X,!Y) :- !:Y = xor(!.Y, X).

generate_uint32(N, RS0, RS) :-
    some ([!Mt,!Mti,!Y]) (
        ustate(!:Mti,!:Mt) = RS0,
        !:Y = 0u32,
	(if !.Mti = n /* generate N words at one time */
	then
	  mt_unif_rand_loop1(0,!Y,!Mt),
	  mt_unif_rand_loop2(n-m,!.Y,_,!Mt),
	  !:Y = (!.Mt ^ elem(n-1) /\ upper_mask) \/ (!.Mt ^ elem(0) /\ lower_mask),
	  !:Mt = !.Mt ^ elem(n-1) := xor(xor(!.Mt ^ elem(m-1), !.Y >> 1), mag01 ^ elem(cast_to_int(!.Y /\ 0x1u32))),
	  !:Mti = 0
        else !:Mti = !.Mti),
	!:Y = !.Mt ^ elem(!.Mti),
	!:Mti = !.Mti+1,
	xor(tempering_shift_u(!.Y), !Y),
	xor(tempering_shift_s(!.Y) /\ tempering_mask_b, !Y),
	xor(tempering_shift_t(!.Y) /\ tempering_mask_c, !Y),
	xor(tempering_shift_l(!.Y), !Y),
	RS = unsafe_promise_unique(ustate(!.Mti,!.Mt)),
	N = !.Y).

generate_uint64(N, !S) :-
    mt_unique.generate_uint32(A0, !S),
    mt_unique.generate_uint32(B0, !S),
    A = uint32.cast_to_uint64(A0),
    B = uint32.cast_to_uint64(B0),
    N = A + (B << 32).

runif(P, U, !R) :-
    mt_unique.generate_uint32(N0, !R),
    C = 2.3283064365386963e-10,
    U = float.cast_from_uint32(N0) * C.
