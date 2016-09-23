#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
from subprocess32 import call,TimeoutExpired
import math
import time
import random
import itertools
from tabulate import tabulate

# import framework
from expfw import *


class ExperimentMC(Experiment):
    """
    Base class for the multi-core tool
    """
    def parse_log(self, contents):
        res = {}
        s = re.compile(r'Time for computing the bisimulation relation: ([\d\.,]+)').findall(contents)
        if len(s) != 1: return None
        res['time'] = float(s[0].replace(',',''))
        s = re.compile('Time needed for signature computation: ([\d\.,]+)').findall(contents)
        if len(s) == 1: res['tsig'] = float(s[0].replace(',',''))
        s = re.compile('Time needed for partition refinement: ([\d\.,]+)').findall(contents)
        if len(s) == 1: res['tref'] = float(s[0].replace(',',''))
        s = re.compile('Time for signature computation: ([\d\.,]+)').findall(contents)
        if len(s) == 1: res['tsig'] = float(s[0].replace(',',''))
        s = re.compile('Time for partition refinement: ([\d\.,]+)').findall(contents)
        if len(s) == 1: res['tref'] = float(s[0].replace(',',''))
        s = re.compile('Time for computing the quotient of the transition relation: ([\d\.,]+)').findall(contents)
        if len(s) == 1: res['tquot'] = float(s[0].replace(',',''))
        s = re.compile('New Markov transition relation: ([\d\.,]+) transitions, ([\d\.,]+) MTBDD nodes').findall(contents)
        if len(s) == 1: res['newmarkov'] = int(s[0][1].replace(',',''))
        s = re.compile('New interactive transition relation: ([\d\.,]+) transitions, ([\d\.,]+) MTBDD nodes').findall(contents)
        if len(s) == 1: res['newtrans'] = int(s[0][1].replace(',',''))
        return res


class ExperimentMCQ(ExperimentMC):
    def parse_log(self, contents):
        res = super(ExperimentMCQ, self).parse_log(contents)
        if res is None: return None
        if 'tquot' not in res: return None
        return res


class ExperimentRW(Experiment):
    """
    Base class for Ralf Wimmer's tools

    Only parses total time done by default.
    """
    def parse_log(self, contents):
        res = {}
        s = re.compile(r'Time for computing the bisimulation relation =[\s]*([\d\.]+)').findall(contents)
        if len(s) != 1: return None
        res['time'] = float(s[0].replace(',',''))
        return res


def float_str(f):
    if str(f) == 'nan': return '--'
    else: return str(int(f))


class ExperimentCTMC1ht(ExperimentMC):
    def __init__(self, name, workers, model):
        self.name = "{}-ht-{}".format(name, workers)
        self.call  = ["./sigrefmc_ht", "-w", str(workers), "-l", "fl"] + model


class ExperimentCTMC1ht_fr(ExperimentMC):
    def __init__(self, name, workers, model):
        self.name = "{}-ht-fr-{}".format(name, workers)
        self.call  = ["./sigrefmc_ht", "-w", str(workers), "-l", "fr"] + model

class ExperimentCTMC1(ExperimentMC):
    def __init__(self, name, workers, model):
        self.name = "{}-{}".format(name, workers)
        self.call  = ["./sigrefmc", "-w", str(workers), "-l", "fl"] + model


class ExperimentCTMC1_fr(ExperimentMC):
    def __init__(self, name, workers, model):
        self.name = "{}-fr-{}".format(name, workers)
        self.call  = ["./sigrefmc", "-w", str(workers), "-l", "fr"] + model


class ExperimentCTMC2(Experiment):
    def __init__(self, name, model):
        self.name = "{}-gmp".format(name)
        self.call  = ["./sigref_gmp"] + model

    def parse_log(self, contents):
        res = {}
        s = re.compile(r'Time for refinement: ([\d\.]+)').findall(contents)
        if len(s) != 1: return None
        res['time'] = float(s[0].replace(',',''))
        return res


class CTMCExperiments():
    def __init__(self):
        self.mc_1 = {}
        self.mc_48 = {}
        self.mcht_1 = {}
        self.mcht_48 = {}
        self.mc_fr_1 = {}
        self.mc_fr_48 = {}
        self.mcht_fr_1 = {}
        self.mcht_fr_48 = {}
        self.rw = {}

        for m,a in self.models:
            self.mc_1[m] = ExperimentCTMC1(name=m, workers=1, model=a)
            self.mc_48[m] = ExperimentCTMC1(name=m, workers=48, model=a)
            self.mcht_1[m] = ExperimentCTMC1ht(name=m, workers=1, model=a)
            self.mcht_48[m] = ExperimentCTMC1ht(name=m, workers=48, model=a)
            self.mc_fr_1[m] = ExperimentCTMC1_fr(name=m, workers=1, model=a)
            self.mc_fr_48[m] = ExperimentCTMC1_fr(name=m, workers=48, model=a)
            self.mcht_fr_1[m] = ExperimentCTMC1ht_fr(name=m, workers=1, model=a)
            self.mcht_fr_48[m] = ExperimentCTMC1ht_fr(name=m, workers=48, model=a)
            self.rw[m] = ExperimentCTMC2(name=m, model=a)

    def __iter__(self):
        # only include the fr (gmp) experiments, and Wimmer's tool
        # dicts = ["mc_1", "mc_48", "mcht_1", "mcht_48", "rw"]
        dicts = ["mc_fr_1", "mc_fr_48", "mcht_fr_1", "mcht_fr_48", "rw"]
        return itertools.chain(*(getattr(self, x).values() for x in dicts))

    def analyse_experiment(self, name, results, timeouts):
        r = {}
        # compute (count,average) for all times
        r['n_mc_1'], r['mc_1'] = online_variance([v['time'] for n, v in results if n==self.mc_1[name].name])[0:2]
        r['n_mc_48'], r['mc_48'] = online_variance([v['time'] for n, v in results if n==self.mc_48[name].name])[0:2]
        r['n_mcht_1'], r['mcht_1'] = online_variance([v['time'] for n, v in results if n==self.mcht_1[name].name])[0:2]
        r['n_mcht_48'], r['mcht_48'] = online_variance([v['time'] for n, v in results if n==self.mcht_48[name].name])[0:2]
        r['n_mc_fr_1'], r['mc_fr_1'] = online_variance([v['time'] for n, v in results if n==self.mc_fr_1[name].name])[0:2]
        r['n_mc_fr_48'], r['mc_fr_48'] = online_variance([v['time'] for n, v in results if n==self.mc_fr_48[name].name])[0:2]
        r['n_mcht_fr_1'], r['mcht_fr_1'] = online_variance([v['time'] for n, v in results if n==self.mcht_fr_1[name].name])[0:2]
        r['n_mcht_fr_48'], r['mcht_fr_48'] = online_variance([v['time'] for n, v in results if n==self.mcht_fr_48[name].name])[0:2]
        r['n_rw'], r['rw'] = online_variance([v['time'] for n, v in results if n==self.rw[name].name])[0:2]
        # compute (average) for signature/refine parts from the (gmp) multi-core tool
        r['sig_1'] = online_variance([v['tsig'] for n, v in results if n==self.mc_fr_1[name].name])[1]
        r['sig_48'] = online_variance([v['tsig'] for n, v in results if n==self.mc_fr_48[name].name])[1]
        r['ref_1'] = online_variance([v['tref'] for n, v in results if n==self.mc_fr_1[name].name])[1]
        r['ref_48'] = online_variance([v['tref'] for n, v in results if n==self.mc_fr_48[name].name])[1]
        # compute speedup for the (gmp) multi-core tool
        if r['mc_fr_48'] > 0: r['speedup'] = r['mc_fr_1']/r['mc_fr_48']
        else: r['speedup'] = float('nan')
        # compute speedup relative to Wimmer's tool
        if r['mc_fr_1'] > 0: r['mc_1_rw'] = r['rw']/r['mc_fr_1']
        else: r['mc_1_rw'] = float('nan')
        return r

    def analyse(self, results, timeouts):
        data = {}
        for name, fn in self.models:
            data[name] = self.analyse_experiment(name, results, timeouts)
        self.data = data
        return data

    def report(self, res=None):
        if res is None: res = self.data

        # Report #experiments, times, speedups
        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name,
                          "{}".format(r['n_rw']),
                          "{}".format(r['n_mc_fr_1']),
                          "{}".format(r['n_mc_fr_48']),
                          "{:<6.2f}".format(r['rw']),
                          "{:<6.2f}".format(r['mc_fr_1']),
                          "{:<6.2f}".format(r['mc_fr_48']),
                          "{:<6.2f}".format(r['mc_1_rw']),
                          "{:<6.2f}".format(r['speedup']),
                          "{:<6.2f}".format(r['speedup']*r['mc_1_rw']),
                          ])

        headers = ["Model      ", "N_rw","N_1","N_48","  T_rw", "  T_1", "  T_48", "Speedup", "ParSpeedup", "TotalSpeedup"]
        print(tabulate(table, headers))

        print()

        # Report sig/ref times, sig/ref percentages, speedups
        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name,
                          "{:<6.2f}".format(r['sig_1']),
                          "{:<6.2f}".format(r['ref_1']),
                          "{:<6.2f}".format(r['sig_48']),
                          "{:<6.2f}".format(r['ref_48']),
                          "{:<6.2f}".format(100.0*r['sig_1']/r['mc_fr_1']),
                          "{:<6.2f}".format(100.0*r['ref_1']/r['mc_fr_1']),
                          "{:<6.2f}".format(100.0*r['sig_48']/r['mc_fr_48']),
                          "{:<6.2f}".format(100.0*r['ref_48']/r['mc_fr_48']),
                          "{:<6.2f}".format(r['sig_1']/r['sig_48'] if r['sig_48'] != 0 else float('nan')),
                          "{:<6.2f}".format(r['ref_1']/r['ref_48'] if r['ref_48'] != 0 else float('nan')),
                          ])
        headers = ["Model", "Tsig_1", "Tref_1", "Tsig_48", "Tref_48", "%", "%", "%", "%", "Ssig", "Sref"]
        print(tabulate(table, headers))

    def report_latex(self, out):
        res = self.data

        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name,
                          "{}".format(self.s1[name][0]) if name in self.s1 else "",
                          "{}".format(self.s1[name][1]) if name in self.s1 else "",
                          "{:<6.2f}".format(r['rw']),
                          "{:<6.2f}".format(r['mc_fr_1']),
                          "{:<6.2f}".format(r['mc_fr_48']),
                          "{:<6.2f}".format(r['mc_1_rw']),
                          "{:<6.2f}".format(r['speedup']),
                          "{:<6.2f}".format(r['speedup']*r['mc_1_rw']),
                          ])

        headers = ["Model", "States", "Blocks", "$T_{rw}$", "$T_{1}$", "$T_{48}$", "Speedup", "ParSpeedup", "TotalSpeedup"]
        table = fixnan(table)
        out.write(tabulate(table, headers, tablefmt='latex_booktabs'))

    models = [
        ("cycling-2", ["models/cycling-2.xctmc"]),
        ("cycling-3", ["models/cycling-3.xctmc"]),
        ("cycling-4", ["models/cycling-4.xctmc"]),
        ("cycling-5", ["models/cycling-5.xctmc"]),
        ("kanban-3", ["models/kanban-3.xctmc"]),
        ("kanban-4", ["models/kanban-4.xctmc"]),
        ("fgf", ["models/fgf.xctmc"]),
        #("multiproc-2-2", ["models/multiproc-2-2.xctmc"]),
        #("multiproc-2-3", ["models/multiproc-2-3.xctmc"]),
        #("multiproc-2-4", ["models/multiproc-2-4.xctmc"]),
        #("multiproc-3-1", ["models/multiproc-3-1.xctmc"]),
        #("multiproc-3-2", ["models/multiproc-3-2.xctmc"]),
        ("p2p-3-5", ["models/p2p-3-5.xctmc"]),
        ("p2p-4-4", ["models/p2p-4-4.xctmc"]),
        ("p2p-4-5", ["models/p2p-4-5.xctmc"]),
        ("p2p-4-6", ["models/p2p-4-6.xctmc"]),
        ("p2p-5-4", ["models/p2p-5-4.xctmc"]),
        ("p2p-5-5", ["models/p2p-5-5.xctmc"]),
        ("p2p-5-6", ["models/p2p-5-6.xctmc"]),
        ("p2p-6-5", ["models/p2p-6-5.xctmc"]),
        ("p2p-7-5", ["models/p2p-7-5.xctmc"]),
        ("polling-02", ["models/polling-02.xctmc"]),
        ("polling-03", ["models/polling-03.xctmc"]),
        ("polling-04", ["models/polling-04.xctmc"]),
        ("polling-05", ["models/polling-05.xctmc"]),
        ("polling-06", ["models/polling-06.xctmc"]),
        ("polling-07", ["models/polling-07.xctmc"]),
        ("polling-08", ["models/polling-08.xctmc"]),
        ("polling-09", ["models/polling-09.xctmc"]),
        ("polling-10", ["models/polling-10.xctmc"]),
        ("polling-11", ["models/polling-11.xctmc"]),
        ("polling-12", ["models/polling-12.xctmc"]),
        ("polling-13", ["models/polling-13.xctmc"]),
        ("polling-14", ["models/polling-14.xctmc"]),
        ("polling-15", ["models/polling-15.xctmc"]),
        ("polling-16", ["models/polling-16.xctmc"]),
        ("polling-17", ["models/polling-17.xctmc"]),
        ("polling-18", ["models/polling-18.xctmc"]),
        ("robot-015", ["models/robot-015.xctmc"]),
        ("robot-020", ["models/robot-020.xctmc"]),
        ("robot-021", ["models/robot-021.xctmc"]),
        ("robot-022", ["models/robot-022.xctmc"]),
        ("robot-023", ["models/robot-023.xctmc"]),
        ("robot-024", ["models/robot-024.xctmc"]),
        ("robot-025", ["models/robot-025.xctmc"]),
        ("robot-026", ["models/robot-026.xctmc"]),
        ("robot-027", ["models/robot-027.xctmc"]),
        ("robot-028", ["models/robot-028.xctmc"]),
        ("robot-029", ["models/robot-029.xctmc"]),
        ("robot-030", ["models/robot-030.xctmc"]),
    ]

    s1 = {
        "cycling-2": (4666, 3511),
        "cycling-3": (57667, 40659),
        "cycling-4": (431101, 282943),
        "cycling-5": (2326666, 1424914),
        "fgf": (80616, 38639),
        "kanban-3": (58400, 58400),
        "kanban-4": (454475, 454475),
        "p2p-5-4": (1048576, 105),
        "p2p-5-5": (33554432, 196),
        "p2p-5-6": (1073741824, 336),
        "p2p-6-5": (1073741824, 266),
        "p2p-7-5": (34359738368, 336),
        "polling-10": (15360, 1536),
        "polling-11": (33792, 3072),
        "polling-12": (73728, 6144),
        "polling-13": (159744, 12288),
        "polling-14": (344064, 24576),
        "polling-15": (737280, 49152),
        "polling-16": (1572864, 98304),
        "polling-17": (3342336, 196608),
        "polling-18": (7077888, 393216),
        "robot-015": (13020, 12810),
        "robot-020": (31160, 30780),
        "robot-025": (61200, 60600),
        "robot-030": (106140, 105270),
    }


class ExperimentCTMC_blocks1(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-blocks1-{}".format(name, workers)
        self.call = ["./sigrefmc", "-w", str(workers), "-q", "block-s1"] + model


class ExperimentCTMC_block(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-block-{}".format(name, workers)
        self.call = ["./sigrefmc", "-w", str(workers), "-q", "block"] + model


class ExperimentCTMC_pick(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-pick-{}".format(name, workers)
        self.call = ["./sigrefmc", "-w", str(workers), "-q", "pick"] + model


class CTMCQExperiments(CTMCExperiments):
    def __init__(self):
        # block standard 1,48
        self.bs1 = {}
        self.bs48 = {}
        # block custom 1,48
        self.b1 = {}
        self.b48 = {}
        # pick 1,48
        self.p1 = {}
        self.p48 = {}

        for m,a in self.models:
            self.bs1[m]  = ExperimentCTMC_blocks1(m, 1, a)
            self.bs48[m] = ExperimentCTMC_blocks1(m, 48, a)
            self.b1[m]   = ExperimentCTMC_block(m, 1, a)
            self.b48[m]  = ExperimentCTMC_block(m, 48, a)
            self.p1[m]   = ExperimentCTMC_pick(m, 1, a)
            self.p48[m]  = ExperimentCTMC_pick(m, 48, a)

    def __iter__(self):
        dicts = ["bs1", "bs48", "b1", "b48", "p1", "p48"]
        return itertools.chain(*(getattr(self, x).values() for x in dicts))

    def analyse_experiment(self, name, results, timeouts):
        r = {}
        # compute (count,average) for all times
        dicts = ["bs1", "bs48", "b1", "b48", "p1", "p48"]
        for d in dicts:
            results_for_d = [v for n, v in results if n == getattr(self, d)[name].name]
            r['n_'+d], r[d] = online_variance([v['tquot'] for v in results_for_d])[0:2]
        # get new Markov size from single-worker instances AND CHECK THEY ARE SAME.
        results_for_b1 = [v for n, v in results if n == self.b1[name].name or n == self.bs1[name].name]
        count, r['bnodes'], check = online_variance([v['newmarkov'] for v in results_for_b1])
        if count > 1 and check > 0:
            print("Warning: variance in block encoding for {}!".format(name))
        results_for_p1 = [v for n, v in results if n == self.p1[name].name]
        count, r['pnodes'], check = online_variance([v['newmarkov'] for v in results_for_p1])
        if count > 1 and check > 0:
            print("Warning: variance in pick encoding for {}!".format(name))
        return r

    def analyse(self, results, timeouts):
        data = {}
        for name, fn in self.models:
            data[name] = self.analyse_experiment(name, results, timeouts)
        self.data = data
        return data

    def report(self, res=None):
        if res is None: res = self.data

        # Report #experiments
        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name,
                          "{}".format(r['n_bs1']),
                          "{}".format(r['n_bs48']),
                          "{}".format(r['n_b1']),
                          "{}".format(r['n_b48']),
                          "{}".format(r['n_p1']),
                          "{}".format(r['n_p48'])
                          ])
        headers = ["Model      ", "#bs1", "#bs48", "#b1", "#b48", "#p1", "#p48"]
        print(tabulate(table, headers))

        print()

        # Report times and MTBDD sizes
        table = []
        for name in sorted(res.keys()):
            r = res[name]

            if r['bs48'] > 0: sbs48 = r['bs1']/r['bs48']
            else: sbs48 = float('nan')
            if r['b48'] > 0: sb48 = r['b1']/r['b48']
            else: sb48 = float('nan')
            if r['p48'] > 0: sp48 = r['p1']/r['p48']
            else: sp48 = float('nan')
            if r['pnodes'] > 0: f = float(r['bnodes'])/r['pnodes']
            else: f = float('nan')

            table.append([name,
                          "{:<6.2f}".format(r['bs1']),
                          "{:<6.2f}".format(r['bs48']),
                          "{:<6.2f}".format(sbs48),
                          "{:<6.2f}".format(r['b1']),
                          "{:<6.2f}".format(r['b48']),
                          "{:<6.2f}".format(sb48),
                          "{:<6.2f}".format(r['p1']),
                          "{:<6.2f}".format(r['p48']),
                          "{:<6.2f}".format(sp48),
                          "{}".format(float_str(r['bnodes'])),
                          "{}".format(float_str(r['pnodes'])),
                          "{:<6.2f}".format(f),
                          ])

        headers = ["Model      ", "T_bs1", "T_bs48", "Sbs", "T_b1", "T_b48", "Sb", "T_p1", "T_p48", "Sp", "B", "P", "f"]
        print(tabulate(table, headers))


class ExperimentLTS_s(ExperimentMC):
    def __init__(self, name, workers, model):
        self.name = "{}-s-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "2", "-w", str(workers)] + model


class ExperimentLTS_b(ExperimentMC):
    def __init__(self, name, workers, model):
        self.name = "{}-b-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "1", "-w", str(workers)] + model


class ExperimentLTS2_s(ExperimentRW):
    def __init__(self, name, model):
        self.name = "{}-15-s".format(name)
        self.call = ["./sigref", "--bisi=2", "--infile={}".format(*model), "--verbosity=1"]


class ExperimentLTS2_b(ExperimentRW):
    def __init__(self, name, model):
        self.name = "{}-15-b".format(name)
        self.call = ["./sigref", "--bisi=1", "--infile={}".format(*model), "--verbosity=1"]


class LTSExperiments:
    def __init__(self):
        self.mcb_1 = {}
        self.mcb_48 = {}
        self.mcs_1 = {}
        self.mcs_48 = {}
        self.rws = {}
        self.rwb = {}

        for m,a in self.models:
            self.mcs_1[m] = ExperimentLTS_s(name=m, workers=1, model=a)
            self.mcs_48[m] = ExperimentLTS_s(name=m, workers=48, model=a)
            self.rws[m] = ExperimentLTS2_s(name=m, model=a)
            self.mcb_1[m] = ExperimentLTS_b(name=m, workers=1, model=a)
            self.mcb_48[m] = ExperimentLTS_b(name=m, workers=48, model=a)
            self.rwb[m] = ExperimentLTS2_b(name=m, model=a)

    def __iter__(self):
        # these models are too big for our tool (well, unless I change the maximum number allowed blocks)
        s1 = [v for k,v in self.mcs_1.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        s2 = [v for k,v in self.mcs_48.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        s3 = [v for k,v in self.rws.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        b1 = self.mcb_1.values()
        b2 = self.mcb_48.values()
        b3 = self.rwb.values()
        return itertools.chain(s1, s2, s3, b1, b2, b3)

    models = [
        ("kanban01", ["models/kanban01.xlts"]),
        ("kanban02", ["models/kanban02.xlts"]),
        ("kanban03", ["models/kanban03.xlts"]),
        ("kanban04", ["models/kanban04.xlts"]),
        ("kanban05", ["models/kanban05.xlts"]),
        ("kanban06", ["models/kanban06.xlts"]),
        ("kanban07", ["models/kanban07.xlts"]),
        ("kanban08", ["models/kanban08.xlts"]),
        ("kanban09", ["models/kanban09.xlts"]),
    ]

    s1 = {
        "kanban01-s": (256, 148),
        "kanban02-s": (63772, 5725),
        "kanban03-s": (1024240, 85356),
        "kanban04-s": (16020316, 778485),
        "kanban05-s": (16772032, 5033631),
        "kanban06-s": (264515056, 25293849),
        "kanban07-s": (268430272, 0),
        "kanban08-s": (4224876912, 0),
        "kanban01-b": (256, 24),
        "kanban02-b": (63772, 206),
        "kanban03-b": (1024240, 872),
        "kanban04-b": (16020316, 2785),
        "kanban05-b": (16772032, 7366),
        "kanban06-b": (264515056, 17010),
        "kanban07-b": (268430272, 35456),
        "kanban08-b": (4224876912, 68217),
    }

    def analyse(self, results, timeouts):
        data = {}
        for name in [name for name, fn in self.models]:
            r = {}
            r['n_mc_1'], r['mc_1'] = online_variance([v['time'] for n, v in results if n==self.mcs_1[name].name])[0:2]
            r['n_mc_48'], r['mc_48'] = online_variance([v['time'] for n, v in results if n==self.mcs_48[name].name])[0:2]
            r['n_rw'], r['rw'] = online_variance([v['time'] for n, v in results if n==self.rws[name].name])[0:2]
            r['sig_1'] = online_variance([v['tsig'] for n, v in results if n==self.mcs_1[name].name])[1]
            r['sig_48'] = online_variance([v['tsig'] for n, v in results if n==self.mcs_48[name].name])[1]
            r['ref_1'] = online_variance([v['tref'] for n, v in results if n==self.mcs_1[name].name])[1]
            r['ref_48'] = online_variance([v['tref'] for n, v in results if n==self.mcs_48[name].name])[1]
            if r['mc_48'] > 0: r['speedup'] = r['mc_1']/r['mc_48']
            else: r['speedup'] = float('nan')
            if r['mc_1'] > 0: r['mc_1_rw'] = r['rw']/r['mc_1']
            else: r['mc_1_rw'] = float('nan')
            data[name+"-s"] = r
            r = {}
            r['n_mc_1'], r['mc_1'] = online_variance([v['time'] for n, v in results if n==self.mcb_1[name].name])[0:2]
            r['n_mc_48'], r['mc_48'] = online_variance([v['time'] for n, v in results if n==self.mcb_48[name].name])[0:2]
            r['n_rw'], r['rw'] = online_variance([v['time'] for n, v in results if n==self.rwb[name].name])[0:2]
            r['sig_1'] = online_variance([v['tsig'] for n, v in results if n==self.mcb_1[name].name])[1]
            r['sig_48'] = online_variance([v['tsig'] for n, v in results if n==self.mcb_48[name].name])[1]
            r['ref_1'] = online_variance([v['tref'] for n, v in results if n==self.mcb_1[name].name])[1]
            r['ref_48'] = online_variance([v['tref'] for n, v in results if n==self.mcb_48[name].name])[1]
            if r['mc_48'] > 0: r['speedup'] = r['mc_1']/r['mc_48']
            else: r['speedup'] = float('nan')
            if r['mc_1'] > 0: r['mc_1_rw'] = r['rw']/r['mc_1']
            else: r['mc_1_rw'] = float('nan')
            data[name+"-b"] = r
        self.data = data
        return data

    def report(self, res=None):
        if res is None: res = self.data

        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name,
                          "{}".format(r['n_rw']),
                          "{}".format(r['n_mc_1']),
                          "{}".format(r['n_mc_48']),
                          "{:<6.2f}".format(r['rw']),
                          "{:<6.2f}".format(r['mc_1']),
                          "{:<6.2f}".format(r['mc_48']),
                          "{:<6.2f}".format(r['mc_1_rw']),
                          "{:<6.2f}".format(r['speedup']),
                          "{:<6.2f}".format(r['speedup']*r['mc_1_rw']),
                          ])

        headers = ["Model      ", "N_rw","N_1","N_48","  T_rw", "  T_1", "  T_48", "Speedup", "ParSpeedup", "TotalSpeedup"]
        print(tabulate(table, headers))

        print()

        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name,
                          "{:<6.2f}".format(r['sig_1']),
                          "{:<6.2f}".format(r['ref_1']),
                          "{:<6.2f}".format(r['sig_48']),
                          "{:<6.2f}".format(r['ref_48']),
                          "{:<6.2f}".format(100.0*r['sig_1']/r['mc_1']),
                          "{:<6.2f}".format(100.0*r['ref_1']/r['mc_1']),
                          "{:<6.2f}".format(100.0*r['sig_48']/r['mc_48']),
                          "{:<6.2f}".format(100.0*r['ref_48']/r['mc_48']),
                          "{:<6.2f}".format(r['sig_1']/r['sig_48'] if r['sig_48'] != 0 else float('nan')),
                          "{:<6.2f}".format(r['ref_1']/r['ref_48'] if r['ref_48'] != 0 else float('nan')),
                          ])
        headers = ["Model", "Tsig_1", "Tref_1", "Tsig_48", "Tref_48", "%", "%", "%", "%", "Ssig", "Sref"]
        print(tabulate(table, headers))

    def report_latex(self, out):
        res = self.data

        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name,
                          "{}".format(self.s1[name][0]) if name in self.s1 else "",
                          "{}".format(self.s1[name][1]) if name in self.s1 else "",
                          "{:<6.2f}".format(r['rw']),
                          "{:<6.2f}".format(r['mc_1']),
                          "{:<6.2f}".format(r['mc_48']),
                          "{:<6.2f}".format(r['mc_1_rw']),
                          "{:<6.2f}".format(r['speedup']),
                          "{:<6.2f}".format(r['speedup']*r['mc_1_rw']),
                          ])

        headers = ["Model", "States", "Blocks", "$T_{rw}$", "$T_{1}$", "$T_{48}$", "Speedup", "ParSpeedup", "TotalSpeedup"]
        table = fixnan(table)
        out.write(tabulate(table, headers, tablefmt='latex_booktabs'))


class ExperimentLTS_s_blocks1(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-s-blocks1-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "2", "-w", str(workers), "-q", "block-s1"] + model


class ExperimentLTS_s_block(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-s-block-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "2", "-w", str(workers), "-q", "block"] + model


class ExperimentLTS_s_pick(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-s-pick-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "2", "-w", str(workers), "-q", "pick"] + model


class ExperimentLTS_b_blocks1(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-b-blocks1-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "1", "-w", str(workers), "-q", "block-s1"] + model


class ExperimentLTS_b_block(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-b-block-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "1", "-w", str(workers), "-q", "block"] + model


class ExperimentLTS_b_pick(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-b-pick-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "1", "-w", str(workers), "-q", "pick"] + model


class LTSQExperiments(LTSExperiments):
    def __init__(self):
        # block standard 1,48
        self.s_bs1 = {}
        self.s_bs48 = {}
        self.b_bs1 = {}
        self.b_bs48 = {}
        # block custom 1,48
        self.s_b1 = {}
        self.s_b48 = {}
        self.b_b1 = {}
        self.b_b48 = {}
        # pick 1,48
        self.s_p1 = {}
        self.s_p48 = {}
        self.b_p1 = {}
        self.b_p48 = {}

        for m,a in self.models:
            self.s_bs1[m]  = ExperimentLTS_s_blocks1(m, 1, a)
            self.s_bs48[m] = ExperimentLTS_s_blocks1(m, 48, a)
            self.s_b1[m]   = ExperimentLTS_s_block(m, 1, a)
            self.s_b48[m]  = ExperimentLTS_s_block(m, 48, a)
            self.s_p1[m]   = ExperimentLTS_s_pick(m, 1, a)
            self.s_p48[m]  = ExperimentLTS_s_pick(m, 48, a)
            self.b_bs1[m]  = ExperimentLTS_b_blocks1(m, 1, a)
            self.b_bs48[m] = ExperimentLTS_b_blocks1(m, 48, a)
            self.b_b1[m]   = ExperimentLTS_b_block(m, 1, a)
            self.b_b48[m]  = ExperimentLTS_b_block(m, 48, a)
            self.b_p1[m]   = ExperimentLTS_b_pick(m, 1, a)
            self.b_p48[m]  = ExperimentLTS_b_pick(m, 48, a)

    def __iter__(self):
        # these models are too big for our tool (well, unless I change the maximum number allowed blocks)
        s1 = [v for k,v in self.s_bs1.items() if k != "kanban06" and k != "kanban07" and k != "kanban08" and k != "kanban09"]
        s2 = [v for k,v in self.s_bs48.items() if k != "kanban06" and k != "kanban07" and k != "kanban08" and k != "kanban09"]
        s3 = [v for k,v in self.s_b1.items() if k != "kanban06" and k != "kanban07" and k != "kanban08" and k != "kanban09"]
        s4 = [v for k,v in self.s_b48.items() if k != "kanban06" and k != "kanban07" and k != "kanban08" and k != "kanban09"]
        s5 = [v for k,v in self.s_p1.items() if k != "kanban06" and k != "kanban07" and k != "kanban08" and k != "kanban09"]
        s6 = [v for k,v in self.s_p48.items() if k != "kanban06" and k != "kanban07" and k != "kanban08" and k != "kanban09"]
        # s1 = [v for k,v in self.s_bs1.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        # s2 = [v for k,v in self.s_bs48.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        # s3 = [v for k,v in self.s_b1.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        # s4 = [v for k,v in self.s_b48.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        # s5 = [v for k,v in self.s_p1.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        # s6 = [v for k,v in self.s_p48.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        b1 = self.b_bs1.values()
        b2 = self.b_bs48.values()
        b3 = self.b_b1.values()
        b4 = self.b_b48.values()
        b5 = self.b_p1.values()
        b6 = self.b_p48.values()
        return itertools.chain(s1, s2, s3, s4, s5, s6, b1, b2, b3, b4, b5, b6)
        # dicts = ["s_bs1", "s_bs48", "s_b1", "s_b48", "s_p1", "s_p48", "b_bs1", "b_bs48", "b_b1", "b_b48", "b_p1", "b_p48"]
        # return itertools.chain(*(getattr(self, x).values() for x in dicts))

    def analyse_experiment(self, name, results, timeouts):
        r = {}
        # compute (count,average) for all times
        dicts = ["s_bs1", "s_bs48", "s_b1", "s_b48", "s_p1", "s_p48", "b_bs1", "b_bs48", "b_b1", "b_b48", "b_p1", "b_p48"]
        for d in dicts:
            results_for_d = [v for n, v in results if n == getattr(self, d)[name].name]
            r['n_'+d], r[d] = online_variance([v['tquot'] for v in results_for_d])[0:2]
            r['t_'+d] = online_variance([v['time'] for v in results_for_d])[1]
            r['sig_'+d] = online_variance([v['tsig'] for v in results_for_d])[1]
            r['ref_'+d] = online_variance([v['tref'] for v in results_for_d])[1]
        # get new Markov size from single-worker instances AND CHECK THEY ARE SAME.
        results_for_s_b1 = [v for n, v in results if n == self.s_b1[name].name or n == self.s_bs1[name].name]
        count, r['s_bnodes'], check = online_variance([v['newtrans'] for v in results_for_s_b1])
        if count > 1: assert check == 0
        results_for_s_p1 = [v for n, v in results if n == self.s_p1[name].name]
        count, r['s_pnodes'], check = online_variance([v['newtrans'] for v in results_for_s_p1])
        if count > 1: assert check == 0
        results_for_b_b1 = [v for n, v in results if n == self.b_b1[name].name or n == self.b_bs1[name].name]
        count, r['b_bnodes'], check = online_variance([v['newtrans'] for v in results_for_b_b1])
        if count > 1: assert check == 0
        results_for_b_p1 = [v for n, v in results if n == self.b_p1[name].name]
        count, r['b_pnodes'], check = online_variance([v['newtrans'] for v in results_for_b_p1])
        if count > 1: assert check == 0
        # get signature refinement times
        dictsb1 = ["b_bs1", "b_b1", "b_p1"]
        dictsb48 = ["b_bs48", "b_b48", "b_p48"]
        dictss1 = ["s_bs1", "s_b1", "s_p1"]
        dictss48 = ["s_bs48", "s_b48", "s_p48"]
        timesb1 = []
        timesb48 = []
        timess1 = []
        timess48 = []
        for d in dictsb1: timesb1 += [v['time'] for n, v in results if n == getattr(self, d)[name].name]
        for d in dictsb48: timesb48 += [v['time'] for n, v in results if n == getattr(self, d)[name].name]
        for d in dictss1: timess1 += [v['time'] for n, v in results if n == getattr(self, d)[name].name]
        for d in dictss48: timess48 += [v['time'] for n, v in results if n == getattr(self, d)[name].name]
        r['time_b1'] = online_variance(timesb1)[1]
        r['time_b48'] = online_variance(timesb48)[1]
        r['time_s1'] = online_variance(timess1)[1]
        r['time_s48'] = online_variance(timess48)[1]
        return r

    def analyse(self, results, timeouts):
        data = {}
        for name, fn in self.models:
            data[name] = self.analyse_experiment(name, results, timeouts)
        self.data = data
        return data

    def report(self, res=None):
        if res is None: res = self.data

        # Report #experiments
        table = []
        for name in sorted(res.keys()):
            r = res[name]

            if r['time_s48'] > 0: stimes = r['time_s1']/r['time_s48']
            else: stimes = float('nan')

            table.append([name+"-s",
                          "{}".format(r['n_s_bs1']),
                          "{}".format(r['n_s_bs48']),
                          "{}".format(r['n_s_b1']),
                          "{}".format(r['n_s_b48']),
                          "{}".format(r['n_s_p1']),
                          "{}".format(r['n_s_p48']),
                          "{:<6.2f}".format(r['time_s1']),
                          "{:<6.2f}".format(r['time_s48']),
                          "{:<6.2f}".format(stimes),
                          ])
        for name in sorted(res.keys()):
            r = res[name]

            if r['time_b48'] > 0: stimes = r['time_b1']/r['time_b48']
            else: stimes = float('nan')

            table.append([name+"-b",
                          "{}".format(r['n_b_bs1']),
                          "{}".format(r['n_b_bs48']),
                          "{}".format(r['n_b_b1']),
                          "{}".format(r['n_b_b48']),
                          "{}".format(r['n_b_p1']),
                          "{}".format(r['n_b_p48']),
                          "{:<6.2f}".format(r['time_b1']),
                          "{:<6.2f}".format(r['time_b48']),
                          "{:<6.2f}".format(stimes),
                          ])
        headers = ["Model      ", "#bs1", "#bs48", "#b1", "#_b48", "#p1", "#p48", "T1", "T48", "S"]
        print(tabulate(table, headers))

        print()

        # Report times and MTBDD sizes
        table = []
        for name in sorted(res.keys()):
            r = res[name]

            if r['s_bs48'] > 0: s_sbs48 = r['s_bs1']/r['s_bs48']
            else: s_sbs48 = float('nan')
            if r['s_b48'] > 0: s_sb48 = r['s_b1']/r['s_b48']
            else: s_sb48 = float('nan')
            if r['s_p48'] > 0: s_sp48 = r['s_p1']/r['s_p48']
            else: s_sp48 = float('nan')
            if r['s_pnodes'] > 0: s_f = float(r['s_bnodes'])/r['s_pnodes']
            else: s_f = float('nan')

            table.append([name+"-s",
                          "{:<6.2f}".format(r['s_bs1']),
                          "{:<6.2f}".format(r['s_bs48']),
                          "{:<6.2f}".format(s_sbs48),
                          "{:<6.2f}".format(r['s_b1']),
                          "{:<6.2f}".format(r['s_b48']),
                          "{:<6.2f}".format(s_sb48),
                          "{:<6.2f}".format(r['s_p1']),
                          "{:<6.2f}".format(r['s_p48']),
                          "{:<6.2f}".format(s_sp48),
                          "{}".format(float_str(r['s_bnodes'])),
                          "{}".format(float_str(r['s_pnodes'])),
                          "{:<6.2f}".format(s_f),
                          ])

        for name in sorted(res.keys()):
            r = res[name]

            if r['b_bs48'] > 0: b_sbs48 = r['b_bs1']/r['b_bs48']
            else: b_sbs48 = float('nan')
            if r['b_b48'] > 0: b_sb48 = r['b_b1']/r['b_b48']
            else: b_sb48 = float('nan')
            if r['b_p48'] > 0: b_sp48 = r['b_p1']/r['b_p48']
            else: b_sp48 = float('nan')
            if r['b_pnodes'] > 0: b_f = float(r['b_bnodes'])/r['b_pnodes']
            else: b_f = float('nan')

            table.append([name+"-b",
                          "{:<6.2f}".format(r['b_bs1']),
                          "{:<6.2f}".format(r['b_bs48']),
                          "{:<6.2f}".format(b_sbs48),
                          "{:<6.2f}".format(r['b_b1']),
                          "{:<6.2f}".format(r['b_b48']),
                          "{:<6.2f}".format(b_sb48),
                          "{:<6.2f}".format(r['b_p1']),
                          "{:<6.2f}".format(r['b_p48']),
                          "{:<6.2f}".format(b_sp48),
                          "{}".format(float_str(r['b_bnodes'])),
                          "{}".format(float_str(r['b_pnodes'])),
                          "{:<6.2f}".format(b_f),
                          ])

        headers = ["Model      ", "T_bs1", "T_bs48", "Sbs", "T_b1", "T_b48", "Sb", "T_p1", "T_p48", "Sp", "B", "P", "f"]
        print(tabulate(table, headers))


class ExperimentIMC_s(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-s-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc"] + kwargs['model']
        self.call += ["-b", "2"]
        self.call += ["-w", str(kwargs['workers'])]
        self.call += ["-l", "fr"]


class ExperimentIMC_b(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-b-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc"] + kwargs['model']
        self.call += ["-b", "1"]
        self.call += ["-w", str(kwargs['workers'])]
        self.call += ["-l", "fr"]


class ExperimentIMC_s_fl(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-s-fl-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc"] + kwargs['model']
        self.call += ["-b", "2"]
        self.call += ["-w", str(kwargs['workers'])]
        self.call += ["-l", "fl"]


class ExperimentIMC_b_fl(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-b-fl-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc"] + kwargs['model']
        self.call += ["-b", "1"]
        self.call += ["-w", str(kwargs['workers'])]
        self.call += ["-l", "fl"]


class ExperimentIMC2_s(ExperimentRW):
    def __init__(self, **kwargs):
        self.name = "{}-15-s".format(kwargs['name'])
        self.call  = ["./sigref"] + ["--infile={}".format(*kwargs['model'])]
        self.call += ["--bisi=2"]
        self.call += ["--verbosity=1"]


class ExperimentIMC2_b(ExperimentRW):
    def __init__(self, **kwargs):
        self.name = "{}-15-b".format(kwargs['name'])
        self.call  = ["./sigref"] + ["--infile={}".format(*kwargs['model'])]
        self.call += ["--bisi=1"]
        self.call += ["--verbosity=1"]


class IMCExperiments:
    def __init__(self):
        self.mcb_1 = {}
        self.mcb_48 = {}
        self.mcs_1 = {}
        self.mcs_48 = {}
        self.rws = {}
        self.rwb = {}

        for m,a in self.models:
            self.mcs_1[m] = ExperimentIMC_s_fl(name=m, workers=1, model=a)
            self.mcs_48[m] = ExperimentIMC_s_fl(name=m, workers=48, model=a)
            self.mcb_1[m] = ExperimentIMC_b_fl(name=m, workers=1, model=a)
            self.mcb_48[m] = ExperimentIMC_b_fl(name=m, workers=48, model=a)
            self.rws[m] = ExperimentIMC2_s(name=m, model=a)
            self.rwb[m] = ExperimentIMC2_b(name=m, model=a)

    def __iter__(self):
        dicts = ["mcs_1", "mcs_48", "mcb_1", "mcb_48", "rws", "rwb"]
        return itertools.chain(*(getattr(self, x).values() for x in dicts))

    models = [
        ("ftwc01", ["models/ftwc01.ximc"]),
        ("ftwc02", ["models/ftwc02.ximc"]),
        ("ftwc03", ["models/ftwc03.ximc"]),
    ]

    s1 = {
        'ftwc01-s': (2048, 1133),
        'ftwc02-s': (32768, 16797),
        'ftwc01-b': (2048, 430),
        'ftwc02-b': (32786, 3886),
    }

    def analyse(self, results, timeouts):
        data = {}
        for name in [name for name, fn in self.models]:
            r = {}
            r['n_mc_1'], r['mc_1'] = online_variance([v['time'] for n, v in results if n==self.mcs_1[name].name])[0:2]
            r['n_mc_48'], r['mc_48'] = online_variance([v['time'] for n, v in results if n==self.mcs_48[name].name])[0:2]
            r['n_rw'], r['rw'] = online_variance([v['time'] for n, v in results if n==self.rws[name].name])[0:2]
            if r['mc_48'] > 0: r['speedup'] = r['mc_1']/r['mc_48']
            else: r['speedup'] = float('nan')
            if r['mc_1'] > 0: r['mc_1_rw'] = r['rw']/r['mc_1']
            else: r['mc_1_rw'] = float('nan')
            data[name+"-s"] = r
            r = {}
            r['n_mc_1'], r['mc_1'] = online_variance([v['time'] for n, v in results if n==self.mcb_1[name].name])[0:2]
            r['n_mc_48'], r['mc_48'] = online_variance([v['time'] for n, v in results if n==self.mcb_48[name].name])[0:2]
            r['n_rw'], r['rw'] = online_variance([v['time'] for n, v in results if n==self.rwb[name].name])[0:2]
            if r['mc_48'] > 0: r['speedup'] = r['mc_1']/r['mc_48']
            else: r['speedup'] = float('nan')
            if r['mc_1'] > 0: r['mc_1_rw'] = r['rw']/r['mc_1']
            else: r['mc_1_rw'] = float('nan')
            data[name+"-b"] = r
        self.data = data
        return data

    def report(self, res=None):
        if res is None: res = self.data

        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name,
                          "{}".format(r['n_rw']),
                          "{}".format(r['n_mc_1']),
                          "{}".format(r['n_mc_48']),
                          "{:<6.2f}".format(r['rw']),
                          "{:<6.2f}".format(r['mc_1']),
                          "{:<6.2f}".format(r['mc_48']),
                          "{:<6.2f}".format(r['mc_1_rw']),
                          "{:<6.2f}".format(r['speedup']),
                          "{:<6.2f}".format(r['speedup']*r['mc_1_rw']),
                          ])

        headers = ["Model      ", "N_rw","N_1","N_48","  T_rw", "  T_1", "  T_48", "Speedup", "ParSpeedup", "TotalSpeedup"]
        print(tabulate(table, headers))

    def report_latex(self, out):
        res = self.data

        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name,
                          "{}".format(self.s1[name][0]) if name in self.s1 else "",
                          "{}".format(self.s1[name][1]) if name in self.s1 else "",
                          "{:<6.2f}".format(r['rw']),
                          "{:<6.2f}".format(r['mc_1']),
                          "{:<6.2f}".format(r['mc_48']),
                          "{:<6.2f}".format(r['mc_1_rw']),
                          "{:<6.2f}".format(r['speedup']),
                          "{:<6.2f}".format(r['speedup']*r['mc_1_rw']),
                          ])

        headers = ["Model", "States", "Blocks", "$T_{rw}$", "$T_{1}$", "$T_{48}$", "Speedup", "ParSpeedup", "TotalSpeedup"]
        table = fixnan(table)
        out.write(tabulate(table, headers, tablefmt='latex_booktabs'))


class ExperimentIMC_s_blocks1(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-s-blocks1-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "2", "-w", str(workers), "-q", "block-s1"] + model


class ExperimentIMC_s_block(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-s-block-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "2", "-w", str(workers), "-q", "block"] + model


class ExperimentIMC_s_pick(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-s-pick-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "2", "-w", str(workers), "-q", "pick"] + model


class ExperimentIMC_b_blocks1(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-b-blocks1-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "1", "-w", str(workers), "-q", "block-s1"] + model


class ExperimentIMC_b_block(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-b-block-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "1", "-w", str(workers), "-q", "block"] + model


class ExperimentIMC_b_pick(ExperimentMCQ):
    def __init__(self, name, workers, model):
        self.name = "{}-b-pick-{}".format(name, workers)
        self.call = ["./sigrefmc", "-b", "1", "-w", str(workers), "-q", "pick"] + model


class IMCQExperiments(IMCExperiments):
    def __init__(self):
        # block standard 1,48
        self.s_bs1 = {}
        self.s_bs48 = {}
        self.b_bs1 = {}
        self.b_bs48 = {}
        # block custom 1,48
        self.s_b1 = {}
        self.s_b48 = {}
        self.b_b1 = {}
        self.b_b48 = {}
        # pick 1,48
        self.s_p1 = {}
        self.s_p48 = {}
        self.b_p1 = {}
        self.b_p48 = {}

        for m,a in self.models:
            self.s_bs1[m]  = ExperimentIMC_s_blocks1(m, 1, a)
            self.s_bs48[m] = ExperimentIMC_s_blocks1(m, 48, a)
            self.s_b1[m]   = ExperimentIMC_s_block(m, 1, a)
            self.s_b48[m]  = ExperimentIMC_s_block(m, 48, a)
            self.s_p1[m]   = ExperimentIMC_s_pick(m, 1, a)
            self.s_p48[m]  = ExperimentIMC_s_pick(m, 48, a)
            self.b_bs1[m]  = ExperimentIMC_b_blocks1(m, 1, a)
            self.b_bs48[m] = ExperimentIMC_b_blocks1(m, 48, a)
            self.b_b1[m]   = ExperimentIMC_b_block(m, 1, a)
            self.b_b48[m]  = ExperimentIMC_b_block(m, 48, a)
            self.b_p1[m]   = ExperimentIMC_b_pick(m, 1, a)
            self.b_p48[m]  = ExperimentIMC_b_pick(m, 48, a)

    def __iter__(self):
        dicts = ["s_bs1", "s_bs48", "s_b1", "s_b48", "s_p1", "s_p48", "b_bs1", "b_bs48", "b_b1", "b_b48", "b_p1", "b_p48"]
        return itertools.chain(*(getattr(self, x).values() for x in dicts))

    def analyse_experiment(self, name, results, timeouts):
        r = {}
        # compute (count,average) for all times
        dicts = ["s_bs1", "s_bs48", "s_b1", "s_b48", "s_p1", "s_p48", "b_bs1", "b_bs48", "b_b1", "b_b48", "b_p1", "b_p48"]
        for d in dicts:
            results_for_d = [v for n, v in results if n == getattr(self, d)[name].name]
            r['n_'+d], r[d] = online_variance([v['tquot'] for v in results_for_d])[0:2]
        # get new Markov size from single-worker instances AND CHECK THEY ARE SAME.
        results_for_s_b1 = [v for n, v in results if n == self.s_b1[name].name or n == self.s_bs1[name].name]
        count, r['s_mbnodes'], check = online_variance([v['newmarkov'] for v in results_for_s_b1])
        if count > 1 and check > 0:
            print("Warning: variance in block encoding for {}!".format(name))
        results_for_s_p1 = [v for n, v in results if n == self.s_p1[name].name]
        count, r['s_mpnodes'], check = online_variance([v['newmarkov'] for v in results_for_s_p1])
        if count > 1 and check > 0:
            print("Warning: variance in pick encoding for {}!".format(name))
        results_for_b_b1 = [v for n, v in results if n == self.b_b1[name].name or n == self.b_bs1[name].name]
        count, r['b_mbnodes'], check = online_variance([v['newmarkov'] for v in results_for_b_b1])
        if count > 1 and check > 0:
            print("Warning: variance in block encoding for {}!".format(name))
        results_for_b_p1 = [v for n, v in results if n == self.b_p1[name].name]
        count, r['b_mpnodes'], check = online_variance([v['newmarkov'] for v in results_for_b_p1])
        if count > 1 and check > 0:
            print("Warning: variance in pick encoding for {}!".format(name))
         # get new interactive relation size from single-worker instances AND CHECK THEY ARE SAME.
        results_for_s_b1 = [v for n, v in results if n == self.s_b1[name].name or n == self.s_bs1[name].name]
        count, r['s_bnodes'], check = online_variance([v['newtrans'] for v in results_for_s_b1])
        if count > 1: assert check == 0
        results_for_s_p1 = [v for n, v in results if n == self.s_p1[name].name]
        count, r['s_pnodes'], check = online_variance([v['newtrans'] for v in results_for_s_p1])
        if count > 1: assert check == 0
        results_for_b_b1 = [v for n, v in results if n == self.b_b1[name].name or n == self.b_bs1[name].name]
        count, r['b_bnodes'], check = online_variance([v['newtrans'] for v in results_for_b_b1])
        if count > 1: assert check == 0
        results_for_b_p1 = [v for n, v in results if n == self.b_p1[name].name]
        count, r['b_pnodes'], check = online_variance([v['newtrans'] for v in results_for_b_p1])
        if count > 1: assert check == 0
        return r

    def analyse(self, results, timeouts):
        data = {}
        for name, fn in self.models:
            data[name] = self.analyse_experiment(name, results, timeouts)
        self.data = data
        return data

    def report(self, res=None):
        if res is None: res = self.data

        # Report #experiments
        table = []
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name+"-s",
                          "{}".format(r['n_s_bs1']),
                          "{}".format(r['n_s_bs48']),
                          "{}".format(r['n_s_b1']),
                          "{}".format(r['n_s_b48']),
                          "{}".format(r['n_s_p1']),
                          "{}".format(r['n_s_p48'])
                          ])
        for name in sorted(res.keys()):
            r = res[name]
            table.append([name+"-b",
                          "{}".format(r['n_b_bs1']),
                          "{}".format(r['n_b_bs48']),
                          "{}".format(r['n_b_b1']),
                          "{}".format(r['n_b_b48']),
                          "{}".format(r['n_b_p1']),
                          "{}".format(r['n_b_p48'])
                          ])
        headers = ["Model      ", "#bs1", "#bs48", "#b1", "#_b48", "#p1", "#p48"]
        print(tabulate(table, headers))

        print()

        # Report times and MTBDD sizes
        table = []
        for name in sorted(res.keys()):
            r = res[name]

            if r['s_bs48'] > 0: s_sbs48 = r['s_bs1']/r['s_bs48']
            else: s_sbs48 = float('nan')
            if r['s_b48'] > 0: s_sb48 = r['s_b1']/r['s_b48']
            else: s_sb48 = float('nan')
            if r['s_p48'] > 0: s_sp48 = r['s_p1']/r['s_p48']
            else: s_sp48 = float('nan')
            if r['s_pnodes'] > 0: s_f = float(r['s_bnodes'])/r['s_pnodes']
            else: s_f = float('nan')

            table.append([name+"-s",
                          "{:<6.2f}".format(r['s_bs1']),
                          "{:<6.2f}".format(r['s_bs48']),
                          "{:<6.2f}".format(s_sbs48),
                          "{:<6.2f}".format(r['s_b1']),
                          "{:<6.2f}".format(r['s_b48']),
                          "{:<6.2f}".format(s_sb48),
                          "{:<6.2f}".format(r['s_p1']),
                          "{:<6.2f}".format(r['s_p48']),
                          "{:<6.2f}".format(s_sp48),
                          "{}".format(float_str(r['s_bnodes'])),
                          "{}".format(float_str(r['s_pnodes'])),
                          "{:<6.2f}".format(s_f),
                          ])

        for name in sorted(res.keys()):
            r = res[name]

            if r['b_bs48'] > 0: b_sbs48 = r['b_bs1']/r['b_bs48']
            else: b_sbs48 = float('nan')
            if r['b_b48'] > 0: b_sb48 = r['b_b1']/r['b_b48']
            else: b_sb48 = float('nan')
            if r['b_p48'] > 0: b_sp48 = r['b_p1']/r['b_p48']
            else: b_sp48 = float('nan')
            if r['b_pnodes'] > 0: b_f = float(r['b_bnodes'])/r['b_pnodes']
            else: b_f = float('nan')

            table.append([name+"-b",
                          "{:<6.2f}".format(r['b_bs1']),
                          "{:<6.2f}".format(r['b_bs48']),
                          "{:<6.2f}".format(b_sbs48),
                          "{:<6.2f}".format(r['b_b1']),
                          "{:<6.2f}".format(r['b_b48']),
                          "{:<6.2f}".format(b_sb48),
                          "{:<6.2f}".format(r['b_p1']),
                          "{:<6.2f}".format(r['b_p48']),
                          "{:<6.2f}".format(b_sp48),
                          "{}".format(float_str(r['b_bnodes'])),
                          "{}".format(float_str(r['b_pnodes'])),
                          "{:<6.2f}".format(b_f),
                          ])

        headers = ["Model      ", "T_bs1", "T_bs48", "Sbs", "T_b1", "T_b48", "Sb", "T_p1", "T_p48", "Sp", "B", "P", "f"]
        print(tabulate(table, headers))


# signature refinement experiments
ctmc = CTMCExperiments()
lts = LTSExperiments()
imc = IMCExperiments()

sr = ExperimentEngine(outdir='out', timeout=3600)
sr += ctmc
sr += lts
sr += imc

# quotient computation experiments
ctmcq = CTMCQExperiments()
ltsq = LTSQExperiments()
imcq = IMCQExperiments()

q = ExperimentEngine(outdir='out-q-blocks-last', timeout=3600)
# q += ctmcq
q += ltsq
# q += imcq

if __name__ == "__main__":
    # select engine
    engine = sr

    if len(sys.argv) > 1:
        if sys.argv[1] == 'run':
            sr.run_experiments()
        elif sys.argv[1] == 'qrun':
            q.run_experiments()
            pass
        elif sys.argv[1] == 'report':
            n, no, results, timeouts = sr.get_results()
            ctmc.analyse(results, timeouts)
            lts.analyse(results, timeouts)
            imc.analyse(results, timeouts)

            ctmc.report()
            print()
            lts.report()
            print()
            imc.report()

            with open('results_ctmc.tex', 'w') as f:
                ctmc.report_latex(f)

            with open('results_lts.tex', 'w') as f:
                lts.report_latex( f)

            with open('results_imc.tex', 'w') as f:
                imc.report_latex(f)
        elif sys.argv[1] == 'qreport':
            n, no, results, timeouts = q.get_results()
            ctmcq.analyse(results, timeouts)
            ltsq.analyse(results, timeouts)
            imcq.analyse(results, timeouts)

            ctmcq.report()
            print()
            ltsq.report()
            print()
            imcq.report()
