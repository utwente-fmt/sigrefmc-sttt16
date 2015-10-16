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

#include <sigref.h>
#include <blocks.h>

int block_length; // number of block variables
BDD block_variables;

VOID_TASK_IMPL_1(prepare_blocks, int, nvars)
{
    block_length = nvars < 25 ? nvars : 25; // Cap it on 2^25 : 33,554,432 blocks max
    uint32_t block_vars[block_length];
    for (int i=0; i<block_length; i++) block_vars[i] = 2000000+i;
    block_variables = sylvan_set_fromarray(block_vars, block_length);
    sylvan_ref(block_variables);
}

TASK_IMPL_1(BDD, encode_block, uint64_t, b)
{
    // for now, assume max 64 bits for a block....
    uint8_t bl[block_length];
    for (int i=0; i<block_length; i++) {
        bl[i] = b & 1 ? 1 : 0;
        b>>=1;
    }
    return sylvan_cube(block_variables, bl);
}

TASK_IMPL_1(uint64_t, decode_block, BDD, block)
{
    uint64_t result = 0;
    uint64_t mask = 1;
    while (block != sylvan_true) {
        BDD b_low = sylvan_low(block);
        if (b_low == sylvan_false) {
            result |= mask;
            block = sylvan_high(block);
        } else {
            block = b_low;
        }
        mask <<= 1;
    }
    return result;
}
