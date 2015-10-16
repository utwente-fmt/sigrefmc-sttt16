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

#ifndef __SYSTEMS__HPP__
#define __SYSTEMS__HPP__

#include <map>
#include <vector>

#include <sylvan.h>
#include <sylvan_obj.hpp>

using namespace sylvan;

namespace sigref {

class StateSystem {
    friend class SystemParser;
    friend class BddLtsParser;

    Bdd states;
    Bdd initialStates;
    std::vector<Bdd> initialPartition;

    Bdd varS;
    Bdd varT;
    Bdd varA;

public:
    Bdd getStates() const { return states; }
    Bdd getInitialStates() const { return initialStates; }
    std::vector<Bdd> getInitialPartition() const { return initialPartition; }
    Bdd getVarS() const { return varS; }
    Bdd getVarT() const { return varT; }
    Bdd getVarA() const { return varA; }
};

class LTS: public StateSystem
{
    friend class SystemParser;
    friend class BddLtsParser;

    std::vector<std::pair<Bdd,Bdd>> transitions;
    Bdd tau;

public:
    std::vector<std::pair<Bdd,Bdd>> getTransitions() const { return transitions; }
    Bdd getTau() const { return tau; }
};

class CTMC: public StateSystem
{
    friend class SystemParser;

    Mtbdd markov_transitions;

public:
    Mtbdd getMarkovTransitions() const { return markov_transitions; }
};

class IMC: public LTS {
    friend class SystemParser;

    Mtbdd markov_transitions;

public:
    Mtbdd getMarkovTransitions() const { return markov_transitions; }
};

} // end namespace sigref

#endif
