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

#ifndef SIGREF_H
#define SIGREF_H

/* Configuration */
extern int bisimulation; // branching
extern int leaftype; // 0 = float, 1 = fraction, 2 = gmp
extern int quotient; // representative
extern int verbosity; // default: no excessive node counting
extern int merge_relations; // merge relations to 1 relation
extern int closure; // 0 = fixpoint, 1 = squaring, 2 = recursive
extern int reachable; // 0 = no, 1 = yes

/* Obtain current wallclock time */
extern double t_start;
double wctime();
#define INFO(s, ...) fprintf(stdout, "[% 8.2f] " s "\n", wctime()-t_start, ##__VA_ARGS__)

#endif
