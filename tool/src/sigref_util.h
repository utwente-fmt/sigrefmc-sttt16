/*
 * Copyright 2015 Formal Methods and Tools, University of Twente
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdint.h>
#include <sylvan.h>

#ifndef SIGREF_UTIL_H
#define SIGREF_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Calculate random number */
uint64_t trng();
VOID_TASK_DECL_0(init_trng);

/**
 * Calculate A & B & C
 */
TASK_DECL_3(BDD, three_and, BDD, BDD, BDD);
#define three_and(a, b, c) CALL(three_and, a, b, c)

/**
 * Substitute each s by a t and vice versa
 */
TASK_DECL_1(MTBDD, swap_prime, MTBDD);
#define swap_prime(set) CALL(swap_prime, set)

/**
 * Compute \BigSatCount sets
 */
TASK_DECL_3(double, big_satcount, MTBDD*, size_t, size_t);
#define big_satcount(sets, count, nvars) CALL(big_satcount, sets, count, nvars)

/**
 * Compute \BigUnion sets
 */
TASK_DECL_2(MTBDD, big_union, MTBDD*, size_t)
#define big_union(sets, count) CALL(big_union, sets, count)

#ifdef __cplusplus
}
#endif

#endif
