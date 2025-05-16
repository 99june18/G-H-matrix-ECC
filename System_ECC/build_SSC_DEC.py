#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build an SSC-DEC H matrix for RS(18,16) over GF(2^16).
Outputs: H_indices.txt, H_symbol.txt, H_binary_32x288.txt
"""

import itertools, random, sys, re
import galois
import numpy as np
import datetime

# -------- parameters --------
SYM_BITS = 16
N_BITS, K_BITS = 288, 256
N_SYM = N_BITS // SYM_BITS            # 18
# MAX_TRIALS_PER_POLY = 2000            # early-exit
RAND_SEED = 2025
random.seed(RAND_SEED)
np.random.seed(RAND_SEED)

# -------- polynomial loader (ASCII only) --------
def load_primitive_polynomials(path="GF_2^16__primitive_polynomial.txt"):
    hex_ok   = re.compile(r"^0x[0-9a-f]+$", re.I)
    power_rx = re.compile(r"[xd]\^?(\d+)", re.I)
    polys = []
    with open(path, "r") as fh:
        for ln, raw in enumerate(fh, 1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            if hex_ok.match(line):
                polys.append(int(line, 16))
                continue
            if re.search(r"[xd]", line, re.I):
                line = re.sub(r"(\+|^)([xd])(\s*($|\+))",
                              r"\1\g<2>1\3", line, flags=re.I)
                degs = [int(m.group(1)) for m in power_rx.finditer(line)]
                if not degs:
                    raise ValueError("Line {} bad poly {}".format(ln, raw))
                val = 1
                for d in degs:
                    val |= 1 << d
                polys.append(val)
                continue
            raise ValueError("Line {} unknown format {}".format(ln, raw))
    return polys

# -------- GF helpers --------
def make_field(poly_int):
    poly  = galois.Poly.Int(poly_int)
    GF    = galois.GF(1 << SYM_BITS, irreducible_poly=poly)
    return GF, GF.primitive_element

def rs_parity_check(GF, alpha, n_sym=(1 << SYM_BITS) - 1):
    idx  = np.arange(n_sym, dtype=np.uint32)
    row0 = alpha ** idx
    row1 = alpha ** (2 * idx)
    return np.vstack([row0, row1])

# vector 1..65535 (once)
A_RANGE = np.arange(1, 1 << SYM_BITS, dtype=np.uint32)

# -------- syndrome routines --------
def all_single_symbol_syndromes(h_col):
    # vectorised multiply (2 × 65 535)
    prod = h_col[:, None] * A_RANGE            # FieldArray (2, 65k)
    # cast to numpy uint32 then to Python int -> hashable
    p0 = prod[0].astype(np.uint32, copy=False).tolist()
    p1 = prod[1].astype(np.uint32, copy=False).tolist()
    return set(zip(p0, p1))

def all_double_bit_syndromes(H, new_id, selected):
    GF     = type(H[0, 0])
    basis  = [GF(1 << b) for b in range(SYM_BITS)]
    h_new  = H[:, new_id]
    out    = set()
    # same symbol
    for b1, b2 in itertools.combinations(basis, 2):
        s = h_new * (b1 + b2)
        out.add((int(s[0]), int(s[1])))
    # cross symbol
    for prev in selected:
        h_prev = H[:, prev]
        for b1 in basis:
            for b2 in basis:
                s = h_new * b1 + h_prev * b2
                out.add((int(s[0]), int(s[1])))
    return out

def try_construct(H):
    available = list(range(H.shape[1]))
    random.shuffle(available)
    chosen, synd = [], set()
    trials_since_pick = 0
    for col in available:
        s1 = all_single_symbol_syndromes(H[:, col])
        if (0, 0) in s1 or synd & s1:
            trials_since_pick += 1
           
        s2 = all_double_bit_syndromes(H, col, chosen)
        if (0, 0) in s2 or synd & s2 or s1 & s2:
            trials_since_pick += 1
            
        # accept
        trials_since_pick = 0
        synd |= s1 | s2
        chosen.append(col)
        if len(chosen) == N_SYM:
            return chosen
    return None

# -------- main loop --------
if __name__ == "__main__":
    polys = load_primitive_polynomials()
    print("Loaded", len(polys), "polynomials")

    for poly in polys:
        print("-- poly 0x%X" % poly)
        GF, alpha = make_field(poly)
        H_full = rs_parity_check(GF, alpha)

        cols = try_construct(H_full)
        if cols is None:
            print("   fail → next poly")
            continue

        print("SUCCESS", cols)
        H_sel = H_full[:, cols]

        # ---------- 파일에 누적 기록 ----------
        import datetime
        tag = "# poly 0x%X  •  %s" % (
            poly, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

        # 1) 열 인덱스
        with open("H_indices.txt", "a") as f:
            f.write("\n%s\n%s\n" % (tag, " ".join(map(str, cols))))

        # 2) 2×18 심볼 행렬
        with open("H_symbol.txt", "a") as f:
            f.write("\n%s\n" % tag)
            for row in H_sel:
                f.write(" ".join("%04X" % int(x) for x in row) + "\n")

        # 3) 32×288 비트 행렬
        with open("H_binary_32x288.txt", "a") as f:
            f.write("\n%s\n" % tag)
            for r in range(H_sel.shape[0]):              # parity rows 0,1
                for plane in reversed(range(SYM_BITS)):  # bit-plane 15…0
                    row_bits = []
                    for c in range(H_sel.shape[1]):      # 18 symbols
                        val = int(H_sel[r, c])
                        for b in reversed(range(SYM_BITS)):  # 15…0
                            bit_val = (val >> b) & 1 if b == plane else 0
                            row_bits.append(str(bit_val))
                    f.write(" ".join(row_bits) + "\n")

        # 한 번 성공하면 바로 종료하려면 다음 줄 유지,
        # 여러 다항식 결과를 계속 누적하려면 주석 처리.
        sys.exit(0)

    print("All polynomials processed.")


    # build: pip install galois(한 번만)
    # python build_SSC_DEC.py
    # 만약 위 커맨드로 업데이트 후 python -V로 확인했는데 여전히 구 버젼일 경우,  alias python=python3.10  alias pip=pip3.10 실행 시 업데이트됨
    # 사용 함수 때문에 python 2.x 버전에서는 실행 시 에러 발생