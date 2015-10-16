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

#include <argp.h>
#include <locale.h>
#include <stdio.h>
#include <sys/time.h>

#ifdef HAVE_PROFILER
#include <gperftools/profiler.h>
#endif

#include <bisimulation.h>
#include <parse_bdd.hpp>
#include <parse_xml.hpp>
#include <sigref.h>
#include <sigref_util.h>

/* Configuration */
static char* model_filename = NULL;
#ifdef HAVE_PROFILER
static char* profile_filename = NULL;
#endif
static int workers = 0; // autodetect

int bisimulation = 1; // branching
int leaftype = 0; // 0 = float, 1 = fraction, 2 = gmp
int quotient = 1; // representative
int verbosity = 0; // default: no excessive node counting
int merge_relations = 0; // merge relations to 1 relation
int closure = 0; // 0 = fixpoint, 1 = squaring, 2 = recursive
int reachable = 0; // 0 = no, 1 = yes

/* argp configuration */
static struct argp_option options[] =
{
    {"workers", 'w', "<workers>", 0, "Number of workers (default=0: autodetect)", 0},
    {"bisi", 'b', "<bisimulation>", 0, "Bisimulation (branching=1, strong=2)", 0},
    {"leaf", 'l', "<leaf type>", 0, "Leaf type (\"floating point\" (default), \"fraction\", \"gmp\")", 0},
    {"verbosity", 'v', "<verbosity>", 0, "Verbosity (default=0, more=1, too much=2)", 0},
    {"quotient", 'q', "<algorithm>", 0, "Quotient algorithm (representative=1, block numbers=2)", 0},
    {"merge-relations", 'm', 0, 0, "Merge transition relations into one transition relation", 0},
    {"closure", 'c', "<closure>", 0, "Closure algorithm (\"fixpoint\", \"squaring\" or \"recursive\")", 0},
    {"reachable", 'r', 0, 0, "Limit partition to reachable states", 0},
#ifdef HAVE_PROFILER
    {"profiler", 'p', "<filename>", 0, "Filename for profiling", 0},
#endif
    {0, 0, 0, 0, 0, 0}
};

static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    switch (key) {
    case 'w':
        workers = atoi(arg);
        break;
    case 'b':
        bisimulation = atoi(arg);
        if (bisimulation<1 || bisimulation>2) argp_usage(state);
        break;
    case 'l':
        if (arg[0] == 'f' && arg[1] == 'l') {
            leaftype = 0;
        } else if (arg[0] == 'f' && arg[1] == 'r') {
            leaftype = 1;
        } else if (arg[0] == 'g') {
            leaftype = 2;
        } else {
            argp_usage(state);
        }
        break;
    case 'q':
        quotient = atoi(arg);
        if (quotient<1 || quotient>2) argp_usage(state);
        break;
    case 'v':
        verbosity = atoi(arg);
        if (verbosity<0 || verbosity>2) argp_usage(state);
        break;
    case 'm':
        merge_relations = 1;
        break;
    case 'r':
        reachable = 1;
        break;
    case 'c':
        if (arg[0] == 'f') {
            closure = 0;
        } else if (arg[0] == 's') {
            closure = 1;
            merge_relations = 1;
        } else if (arg[0] == 'r') {
            closure = 2;
            merge_relations = 1;
        } else {
            argp_usage(state);
        }
        break;
#ifdef HAVE_PROFILER
    case 'p':
        profile_filename = arg;
        break;
#endif
    case ARGP_KEY_ARG:
        if (state->arg_num >= 1) argp_usage(state);
        model_filename = arg;
        break;
    case ARGP_KEY_END:
        if (state->arg_num < 1) argp_usage(state);
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, "<model>", 0, 0, 0, 0};

double t_start;

double
wctime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

VOID_TASK_0(gc_start)
{
    INFO("(GC) Starting garbage collection...");
}

VOID_TASK_0(gc_end)
{
    INFO("(GC) Garbage collection done.");
}

VOID_TASK_1(main_lace, void*, arg)
{
    setlocale(LC_NUMERIC, "en_US.utf-8");

    sylvan_init_package(1LL<<26, 1LL<<31, 1LL<<25, 1LL<<30);
    sylvan_init_bdd(3);
    sylvan_init_mtbdd();
    sylvan_gc_add_mark(0, TASK(gc_start));
    sylvan_gc_add_mark(40, TASK(gc_end));

    CALL(init_trng);

    t_start = wctime();

    /* Initialization done, now read model from file */

    sigref::SystemType sysType;
    sigref::LTS lts;
    sigref::CTMC ctmc;
    sigref::IMC imc;

    const char *dot = strrchr(model_filename, '.');
    if (dot) {
        if (strcmp(dot+1, "bdd") == 0) {
            sigref::BddLtsParser parser(model_filename);
            sysType = sigref::lts_type;
            lts = *parser.getLTS();
        } else if ((strcmp(dot+1, "xlts") == 0) || (strcmp(dot+1, "xctmc") == 0) || (strcmp(dot+1, "ximc") == 0) || (strcmp(dot+1, "xml") == 0)) {
            sigref::SystemParser reader(model_filename, 0, leaftype == 0 ? sigref::float_type : sigref::simple_fraction_type);
            sysType = reader.getType();
            if (sysType == sigref::lts_type) {
                lts = *reader.getLTS();
            } else if (sysType == sigref::ctmc_type) {
                ctmc = *reader.getCTMC();
            } else {
                imc = *reader.getIMC();
            }
        } else {
            fprintf(stderr, "Unknown extension '%s'!\n", dot+1);
            return;
        }
    } else {
        fprintf(stderr, "Unknown extension ''!\n");
        return;
    }

    INFO("Finished reading system from %s.", model_filename);

#ifdef HAVE_PROFILER
    if (profile_filename != NULL) ProfilerStart(profile_filename);
#endif

    if (sysType == sigref::lts_type && bisimulation == 1) CALL(min_lts_branching, lts);
    else if (sysType == sigref::lts_type && bisimulation == 2) CALL(min_lts_strong, lts);
    else if (sysType == sigref::ctmc_type) CALL(min_ctmc, ctmc);
    else if (sysType == sigref::imc_type && bisimulation == 1) CALL(min_imc_branching, imc);
    else if (sysType == sigref::imc_type && bisimulation == 2) CALL(min_imc_strong, imc);
    else fprintf(stderr, "Unsupported system type!\n");

#ifdef HAVE_PROFILER
    if (profile_filename != NULL) ProfilerStop();
#endif

    (void)arg;
}

int
main(int argc, char **argv)
{
    argp_parse(&argp, argc, argv, 0, 0, 0);

    // boot program
    lace_init(workers, 1024*1024*16);
    lace_startup(0, TASK(main_lace), NULL);

    return 0;
}
