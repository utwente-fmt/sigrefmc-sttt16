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

#include <sigref_util.h>
#include <sylvan_common.h>

/* Calculate random height */
DECLARE_THREAD_LOCAL(thread_rng, uint64_t);

uint64_t trng()
{
    LOCALIZE_THREAD_LOCAL(thread_rng, uint64_t);
    thread_rng = 2862933555777941757ULL * thread_rng + 3037000493ULL;
    SET_THREAD_LOCAL(thread_rng, thread_rng);
    return thread_rng;
}

VOID_TASK_0(init_trng_par)
{
    SET_THREAD_LOCAL(thread_rng, (((uint64_t)rand()) << 32 | rand()));
}

VOID_TASK_IMPL_0(init_trng)
{
    INIT_THREAD_LOCAL(thread_rng);
    TOGETHER(init_trng_par);
}

TASK_IMPL_3(BDD, three_and, BDD, a, BDD, b, BDD, c)
{
    if (a == sylvan_false || b == sylvan_false || c == sylvan_false) return sylvan_false;

    if (a == sylvan_true) return sylvan_and(b, c);
    if (b == sylvan_true) return sylvan_and(a, c);
    if (c == sylvan_true) return sylvan_and(a, b);

    BDD result;
    if (cache_get(a|(260LL<<42), b, c, &result)) return result;

    sylvan_gc_test();

    BDDVAR a_var = sylvan_var(a);
    BDDVAR b_var = sylvan_var(b);
    BDDVAR c_var = sylvan_var(c);

    BDD var = a_var;
    if (var > b_var) var = b_var;
    if (var > c_var) var = c_var;

    BDD a_low, a_high;
    if (var == a_var) {
        a_low = sylvan_low(a);
        a_high = sylvan_high(a);
    } else {
        a_low = a_high = a;
    }

    BDD b_low, b_high;
    if (var == b_var) {
        b_low = sylvan_low(b);
        b_high = sylvan_high(b);
    } else {
        b_low = b_high = b;
    }

    BDD c_low, c_high;
    if (var == c_var) {
        c_low = sylvan_low(c);
        c_high = sylvan_high(c);
    } else {
        c_low = c_high = c;
    }

    bdd_refs_spawn(SPAWN(three_and, a_low, b_low, c_low));
    BDD high = bdd_refs_push(CALL(three_and, a_high, b_high, c_high));
    BDD low = bdd_refs_sync(SYNC(three_and));
    result = sylvan_makenode(var, low, high);
    bdd_refs_pop(1);

    cache_put(a|(260LL<<42), b, c, result);
    return result;
}

TASK_IMPL_1(MTBDD, swap_prime, MTBDD, set)
{
    if (mtbdd_isleaf(set)) return set;

    // TODO: properly ignore action/block variables
    if (mtbdd_getvar(set) >= 99999) return set;

    MTBDD result;
    if (cache_get(set|(258LL<<42), set, 0, &result)) return result;

    sylvan_gc_test();

    mtbdd_refs_spawn(SPAWN(swap_prime, mtbdd_getlow(set)));
    MTBDD high = mtbdd_refs_push(CALL(swap_prime, mtbdd_gethigh(set)));
    MTBDD low = mtbdd_refs_sync(SYNC(swap_prime));
    result = mtbdd_makenode(sylvan_var(set)^1, low, high);
    mtbdd_refs_pop(1);

    cache_put(set|(258LL<<42), set, 0, result);
    return result;
}

TASK_IMPL_3(double, big_satcount, MTBDD*, dds, size_t, count, size_t, nvars)
{
    if (count == 1) return mtbdd_satcount(*dds, nvars);
    SPAWN(big_satcount, dds, count/2, nvars);
    double result = big_satcount(dds+count/2, count-count/2, nvars);
    return result + SYNC(big_satcount);
}

TASK_IMPL_2(MTBDD, big_union, MTBDD*, sets, size_t, count)
{
    if (count == 1) return *sets;
    mtbdd_refs_spawn(SPAWN(big_union, sets, count/2));
    MTBDD right = mtbdd_refs_push(CALL(big_union, sets+count/2, count-count/2));
    MTBDD left = mtbdd_refs_push(mtbdd_refs_sync(SYNC(big_union)));
    MTBDD result = mtbdd_plus(left, right);
    mtbdd_refs_pop(2);
    return result;
}
