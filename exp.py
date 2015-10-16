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


class Experiment:
    NOTDONE=0
    DONE=1
    TIMEOUT=2

    def parse_logfile(self, handle):
        return None, None
    
    def get_status(self, filename):
        if (os.path.isfile(filename)):
            with open(filename, 'r') as res:
                res = self.parse_logfile(res)
                if res is not None: return Experiment.DONE, res['time'], res

        timeout_filename = "{}.timeout".format(filename)
        if (os.path.isfile(timeout_filename)):
            with open(timeout_filename, 'r') as to:
                return Experiment.TIMEOUT, int(to.read()), None

        return Experiment.NOTDONE, 0, None

    def run_experiment(self, timeout, filename):
        status, value, db = self.get_status(filename)
        if status == Experiment.DONE: return
        if status == Experiment.TIMEOUT and value >= int(timeout): return

        # remove output and timeout files
        if os.path.isfile(filename): os.unlink(filename)
        timeout_filename = "{}.timeout".format(filename)
        if os.path.isfile(timeout_filename): os.unlink(timeout_filename)

        print("Performing {}... ".format(self.name), end='')
        sys.stdout.flush()

        try:
            with open(filename, 'w+') as out:
                call(self.call, stdout=out, stderr=out, timeout=timeout)
        except KeyboardInterrupt:
            os.unlink(filename)
            print("interrupted!")
            sys.exit()
        except OSError:
            os.unlink(filename)
            print("OS failure! (missing executable?)")
            sys.exit()
        except TimeoutExpired:
            with open(timeout_filename, 'w') as to: to.write(str(timeout))
            print("timeout!")
        else:
            status, value, db = self.get_status(filename)
            if status == Experiment.DONE: print("done; {}!".format(value))
            else: print("not done!")
        time.sleep(2)


class ExperimentSet:
    def __init__(self, **kwargs):
        self.experiments = []
        if 'outdir' in kwargs: self.outdir = kwargs['outdir']
        else: self.outdir = 'out'
        if 'timeout' in kwargs: self.timeout = int(kwargs['timeout'])
        else: self.timeout = 1200

    def get_results(self):
        results = []
        timeouts = []
        yes = 0
        no = 0
        timeout = 0

        for i in itertools.count():
            stop = True
            for e in self.experiments:
                status, value, db = e.get_status("{}/{}-{}".format(self.outdir, e.name, i))
                if not status == Experiment.NOTDONE: stop = False
                if status == Experiment.DONE: results.append((e.name, db))
                elif status == Experiment.TIMEOUT: timeouts.append((e.name, value))
                else: no += 1
            if stop: return i, no - len(self.experiments), results, timeouts

    def run_experiments(self):
        if not os.path.exists(self.outdir): os.makedirs(self.outdir)
        for i in itertools.count():
            n, not_done, results, timeouts = self.get_results()
            print("In {} repetitions, {} succesful, {} timeouts, {} not done.".format(n, len(results), len(timeouts), not_done))

            print("Running iteration {}.".format(i))
            random.shuffle(self.experiments)
            for e in self.experiments: e.run_experiment(self.timeout, "{}/{}-{}".format(self.outdir, e.name, i))

    def append(self, experiment):
        self.experiments.append(experiment)


def online_variance(data):
    n = 0
    mean = 0
    M2 = 0

    for x in data:
        n = n + 1
        delta = x - mean
        mean = mean + delta/n
        M2 = M2 + delta*(x - mean)

    if n < 1: return n, float('nan'), float('nan')
    if n < 2: return n, mean, float('nan')

    variance = M2/(n - 1)
    return n, mean, variance


def fixnan(table):
    return [['--' if s.strip() == 'nan' else s for s in row] for row in table]


class ExperimentMC(Experiment):
    def parse_logfile(self, handle):
        res = {}
        contents = handle.read()
        s = re.compile(r'Time for computing the bisimulation relation: ([\d\.,]+)').findall(contents)
        if len(s) != 1: return None
        res['time'] = float(s[0].replace(',',''))
        s = re.compile('Time needed for signature computation: ([\d\.,]+)').findall(contents)
        if len(s) == 1: res['tsig'] = float(s[0].replace(',',''))
        s = re.compile('Time needed for partition refinement: ([\d\.,]+)').findall(contents)
        if len(s) == 1: res['tref'] = float(s[0].replace(',',''))
        return res


class ExperimentRW(Experiment):
    def parse_logfile(self, handle):
        res = {}
        s = re.compile(r'Time for computing the bisimulation relation =[\s]*([\d\.]+)').findall(handle.read())
        if len(s) != 1: return None
        res['time'] = float(s[0].replace(',',''))
        return res


class ExperimentCTMC1ht(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-ht-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc_ht"] + kwargs['model']
        self.call += ["-w", str(kwargs['workers'])]
        self.call += ["-l", "fl"]


class ExperimentCTMC1ht_fr(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-ht-fr-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc_ht"] + kwargs['model']
        self.call += ["-w", str(kwargs['workers'])]
        self.call += ["-l", "fr"]

class ExperimentCTMC1(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc"] + kwargs['model']
        self.call += ["-w", str(kwargs['workers'])]
        self.call += ["-l", "fl"]


class ExperimentCTMC1_fr(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-fr-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc"] + kwargs['model']
        self.call += ["-w", str(kwargs['workers'])]
        self.call += ["-l", "fr"]


class ExperimentCTMC2(Experiment):
    def __init__(self, **kwargs):
        self.name = kwargs['name'] + "-gmp"
        self.call  = ["./sigref_gmp"] + kwargs['model']

    def parse_logfile(self, handle):
        res = {}
        s = re.compile(r'Time for refinement: ([\d\.]+)').findall(handle.read())
        if len(s) != 1: return None
        res['time'] = float(s[0].replace(',',''))
        return res


class CTMCExperiments():
    def __init__(self):
        self.mc_1 = {}
        self.mc_48 = {}
        self.mcht_1 = {}
        self.mcht_48 = {}
        self.rw = {}

        for m,a in self.models:
            #self.mc_1[m] = ExperimentCTMC1(name=m, workers=1, model=a)
            #self.mc_48[m] = ExperimentCTMC1(name=m, workers=48, model=a)
            #self.mcht_1[m] = ExperimentCTMC1ht(name=m, workers=1, model=a)
            #self.mcht_48[m] = ExperimentCTMC1ht(name=m, workers=48, model=a)
            self.mc_1[m] = ExperimentCTMC1_fr(name=m, workers=1, model=a)
            self.mc_48[m] = ExperimentCTMC1_fr(name=m, workers=48, model=a)
            #self.mcht_1[m] = ExperimentCTMC1ht_fr(name=m, workers=1, model=a)
            #self.mcht_48[m] = ExperimentCTMC1ht_fr(name=m, workers=48, model=a)
            self.rw[m] = ExperimentCTMC2(name=m, model=a)

    def add_to_set(self, experiments):
        experiments.experiments += self.mc_1.values()
        experiments.experiments += self.mc_48.values()
        experiments.experiments += self.mcht_1.values()
        experiments.experiments += self.mcht_48.values()
        experiments.experiments += self.rw.values()

    def analyse(self, results, timeouts):
        res = {}
        for name in [name for name, fn in self.models]:
            r = {}
            #r['n_mcht_1'], r['mcht_1'] = online_variance([v['time'] for n, v in results if n==self.mcht_1[name].name])[0:2]
            #r['n_mcht_48'], r['mcht_48'] = online_variance([v['time'] for n, v in results if n==self.mcht_48[name].name])[0:2]
            r['n_mc_1'], r['mc_1'] = online_variance([v['time'] for n, v in results if n==self.mc_1[name].name])[0:2]
            r['n_mc_48'], r['mc_48'] = online_variance([v['time'] for n, v in results if n==self.mc_48[name].name])[0:2]
            r['n_rw'], r['rw'] = online_variance([v['time'] for n, v in results if n==self.rw[name].name])[0:2]
            r['sig_1'] = online_variance([v['tsig'] for n, v in results if n==self.mc_1[name].name])[1]
            r['sig_48'] = online_variance([v['tsig'] for n, v in results if n==self.mc_48[name].name])[1]
            r['ref_1'] = online_variance([v['tref'] for n, v in results if n==self.mc_1[name].name])[1]
            r['ref_48'] = online_variance([v['tref'] for n, v in results if n==self.mc_48[name].name])[1]
            if r['mc_48'] > 0: r['speedup'] = r['mc_1']/r['mc_48']
            else: r['speedup'] = float('nan')
            if r['mc_1'] > 0: r['mc_1_rw'] = r['rw']/r['mc_1']
            else: r['mc_1_rw'] = float('nan')
            res[name] = r
        return res

    def report(self, results, timeouts):
        res = self.analyse(results, timeouts)
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
        headers = ["Model", "Tsig_1", "Tref_1", "Tsig_48", "Tref_48", "", "", "", "", "Ssig", "Sref"]
        print(tabulate(table, headers))


    def report_latex(self, results, timeouts, out):
        res = self.analyse(results, timeouts)
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


class ExperimentLTS_s(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-s-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc"] + kwargs['model']
        self.call += ["-b", "2"]
        self.call += ["-w", str(kwargs['workers'])]


class ExperimentLTS_b(ExperimentMC):
    def __init__(self, **kwargs):
        self.name = "{}-b-{}".format(kwargs['name'], kwargs['workers'])
        self.call  = ["./sigrefmc"] + kwargs['model']
        self.call += ["-b", "1"]
        self.call += ["-w", str(kwargs['workers'])]


class ExperimentLTS2_s(ExperimentRW):
    def __init__(self, **kwargs):
        self.name = "{}-15-s".format(kwargs['name'])
        self.call  = ["./sigref"] + ["--infile={}".format(*kwargs['model'])]
        self.call += ["--bisi=2"]
        self.call += ["--verbosity=1"]


class ExperimentLTS2_b(ExperimentRW):
    def __init__(self, **kwargs):
        self.name = "{}-15-b".format(kwargs['name'])
        self.call  = ["./sigref"] + ["--infile={}".format(*kwargs['model'])]
        self.call += ["--bisi=1"]
        self.call += ["--verbosity=1"]


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

    def add_to_set(self, experiments):
        # these models are too big for our tool (well, unless I change the maximum number allowed blocks)
        experiments.experiments += [v for k,v in self.mcs_1.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        experiments.experiments += [v for k,v in self.mcs_48.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        experiments.experiments += [v for k,v in self.rws.items() if k != "kanban07" and k != "kanban08" and k != "kanban09"]
        experiments.experiments += self.mcb_1.values()
        experiments.experiments += self.mcb_48.values()
        experiments.experiments += self.rwb.values()

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
        res = {}
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
            res[name+"-s"] = r
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
            res[name+"-b"] = r
        return res

    def report(self, results, timeouts):
        res = self.analyse(results, timeouts)
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
        headers = ["Model", "Tsig_1", "Tref_1", "Tsig_48", "Tref_48", "", "", "", "", "Ssig", "Sref"]
        print(tabulate(table, headers))



    def report_latex(self, results, timeouts, out):
        res = self.analyse(results, timeouts)
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

    def add_to_set(self, experiments):
        experiments.experiments += self.mcs_1.values()
        experiments.experiments += self.mcs_48.values()
        experiments.experiments += self.mcb_1.values()
        experiments.experiments += self.mcb_48.values()
        experiments.experiments += self.rws.values()
        experiments.experiments += self.rwb.values()

    models = [
        ("ftwc01", ["models/ftwc01.ximc"]),
        ("ftwc02", ["models/ftwc02.ximc"]),
        #("ftwc03", ["models/ftwc03.ximc"]),
    ]

    s1 = {
        'ftwc01-s': (2048, 1133),
        'ftwc02-s': (32768, 16797),
        'ftwc01-b': (2048, 430),
        'ftwc02-b': (32786, 3886),
    }

    def analyse(self, results, timeouts):
        res = {}
        for name in [name for name, fn in self.models]:
            r = {}
            r['n_mc_1'], r['mc_1'] = online_variance([v['time'] for n, v in results if n==self.mcs_1[name].name])[0:2]
            r['n_mc_48'], r['mc_48'] = online_variance([v['time'] for n, v in results if n==self.mcs_48[name].name])[0:2]
            r['n_rw'], r['rw'] = online_variance([v['time'] for n, v in results if n==self.rws[name].name])[0:2]
            if r['mc_48'] > 0: r['speedup'] = r['mc_1']/r['mc_48']
            else: r['speedup'] = float('nan')
            if r['mc_1'] > 0: r['mc_1_rw'] = r['rw']/r['mc_1']
            else: r['mc_1_rw'] = float('nan')
            res[name+"-s"] = r
            r = {}
            r['n_mc_1'], r['mc_1'] = online_variance([v['time'] for n, v in results if n==self.mcb_1[name].name])[0:2]
            r['n_mc_48'], r['mc_48'] = online_variance([v['time'] for n, v in results if n==self.mcb_48[name].name])[0:2]
            r['n_rw'], r['rw'] = online_variance([v['time'] for n, v in results if n==self.rwb[name].name])[0:2]
            if r['mc_48'] > 0: r['speedup'] = r['mc_1']/r['mc_48']
            else: r['speedup'] = float('nan')
            if r['mc_1'] > 0: r['mc_1_rw'] = r['rw']/r['mc_1']
            else: r['mc_1_rw'] = float('nan')
            res[name+"-b"] = r
        return res

    def report(self, results, timeouts):
        res = self.analyse(results, timeouts)
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


    def report_latex(self, results, timeouts, out):
        res = self.analyse(results, timeouts)
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


if __name__ == "__main__":
    engine = ExperimentSet(outdir='out', timeout=3600)
    ctmc = CTMCExperiments()
    lts = LTSExperiments()
    imc = IMCExperiments()
    ctmc.add_to_set(engine)
    lts.add_to_set(engine)
    imc.add_to_set(engine)

    if len(sys.argv) > 1:
        if sys.argv[1] == 'run':
            engine.run_experiments()
        elif sys.argv[1] == 'report':
            n, no, results, timeouts = engine.get_results()
            ctmc.report(results, timeouts)
            print()
            lts.report(results, timeouts)
            print()
            imc.report(results, timeouts)

            with open('results_ctmc.tex', 'w') as f:
                ctmc.report_latex(results, timeouts, f)

            with open('results_lts.tex', 'w') as f:
                lts.report_latex(results, timeouts, f)

            with open('results_imc.tex', 'w') as f:
                imc.report_latex(results, timeouts, f)
