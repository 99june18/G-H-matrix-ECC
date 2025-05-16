#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fault-sim for RS(18,16) SSC-DEC code
Outputs CE/DUE/SDC rates for SBE, DBE, TBE, SSE, DSE.
"""

import argparse, random, numpy as np, galois
from collections import Counter

# ---------- load H (symbol) ----------
def load_symbol_H(path="H_symbol.txt"):
    rows = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            rows.append([int(x, 16) for x in ln.split()])
    if len(rows) != 2 or len(rows[0]) != 18:
        raise ValueError("H_symbol.txt must be 2 × 18 hex matrix")
    return np.array(rows, dtype=np.uint32)

H_sym_int = load_symbol_H()                       # 2 × 18 int
GF = galois.GF(2**16)
H_sym = GF(H_sym_int)                             # FieldArray

# ---------- pre-compute syndrome look-ups ----------
SYM_BITS = 16
basis16 = [GF(1 << b) for b in range(SYM_BITS)]

single_tbl = {}      # {(s0,s1): col_id}
double_tbl = {}      # {(s0,s1): 'dummy'}

# 1-symbol syndromes
for col in range(18):
    h = H_sym[:, col]
    for a in GF.Range(1, 1<<16):
        s = (int(h[0]*a), int(h[1]*a))
        single_tbl[s] = col

# 2-bit syndromes (same 또는 cross symbol)
for col in range(18):
    h_new = H_sym[:, col]
    # same-symbol 2-bit
    for b1 in range(SYM_BITS):
        for b2 in range(b1+1, SYM_BITS):
            s = h_new*basis16[b1] + h_new*basis16[b2]
            double_tbl[(int(s[0]), int(s[1]))] = True
    # cross-symbol
    for prev in range(col):            # 이미 테이블에 있으므로 중복 방지
        h_prev = H_sym[:, prev]
        for b1 in range(SYM_BITS):
            for b2 in range(SYM_BITS):
                s = h_new*basis16[b1] + h_prev*basis16[b2]
                double_tbl[(int(s[0]), int(s[1]))] = True

def classify_error(symbol_vec):
    """
    symbol_vec : length-18 list of GF elements (error pattern)
    returns CE / DUE / SDC
    """
    s = H_sym @ GF(symbol_vec).T
    st = (int(s[0]), int(s[1]))
    if st == (0, 0):
        return "SDC"
    if st in single_tbl or st in double_tbl:
        return "CE"
    return "DUE"

# ---------- error pattern generators ----------
def sbe():
    """single bit error"""
    sym = [GF.Zero()] * 18
    col = random.randrange(18)
    bit = random.randrange(16)
    sym[col] = GF(1 << bit)
    return sym

def dbe():
    """double bit error (distinct bits)"""
    while True:
        col1, col2 = random.randrange(18), random.randrange(18)
        bit1, bit2 = random.randrange(16), random.randrange(16)
        if col1 != col2 or bit1 != bit2:
            break
    sym = [GF.Zero()] * 18
    sym[col1] += GF(1 << bit1)
    sym[col2] += GF(1 << bit2)
    return sym

def tbe():
    """triple bit error"""
    sym = [GF.Zero()] * 18
    picks = set()
    while len(picks) < 3:
        picks.add((random.randrange(18), random.randrange(16)))
    for col, bit in picks:
        sym[col] += GF(1 << bit)
    return sym

def sse():
    """single-symbol arbitrary 16-bit value"""
    sym = [GF.Zero()] * 18
    col = random.randrange(18)
    val = random.randrange(1, 1<<16)
    sym[col] = GF(val)
    return sym

def dse():
    """two symbols arbitrary values"""
    sym = [GF.Zero()] * 18
    cols = random.sample(range(18), 2)
    sym[cols[0]] = GF(random.randrange(1, 1<<16))
    sym[cols[1]] = GF(random.randrange(1, 1<<16))
    return sym

pattern_funcs = {
    "SBE": sbe,
    "DBE": dbe,
    "TBE": tbe,
    "SSE": sse,
    "DSE": dse,
}

# ---------- Monte-Carlo simulation ----------
def run_sim(trials):
    results = {p: Counter() for p in pattern_funcs}
    for p_name, gen in pattern_funcs.items():
        for _ in range(trials):
            vec = gen()
            outcome = classify_error(vec)
            results[p_name][outcome] += 1
    return results

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=10000,
                    help="trials per pattern (default 10 000)")
    args = ap.parse_args()

    res = run_sim(args.trials)
    print(f"Trials per pattern: {args.trials}\n")
    print(" Pattern |    CE     |    DUE    |    SDC")
    print("---------+-----------+-----------+-----------")
    for p in ["SBE","DBE","TBE","SSE","DSE"]:
        tot = sum(res[p].values())
        ce  = res[p]["CE"]  / tot
        due = res[p]["DUE"] / tot
        sdc = res[p]["SDC"] / tot
        print(f" {p:7s}| {ce:9.6f} | {due:9.6f} | {sdc:9.6f}")


# 실행 python Fault_sim_RS.py --trials 100000
# git pull test