#!/usr/bin/env python3
"""Phase 2: detection power of Syn2b's two channels, measured against exact truth.

For each (event type, event length, substitution load, enzyme panel) cell we build
a genome carrying K non-overlapping events of that exact length at random
positions, digest it, and compare with the reference. Truth is recorded at
construction time, so detection is checked, not assumed.

Two questions per cell:
  junction channel   -- was the event's boundary pair reported as junctions?
  orientation channel-- did the inverted extent come back unbiased?
"""
import os, sys, json, random, subprocess, csv

BIN  = "/Users/shihuang/Downloads/Syn2b/target/release/Syn2b"
REF  = "/Users/shihuang/Downloads/Syn2bANI/prototype/mg1655.fasta"
WORK = os.path.dirname(os.path.abspath(__file__)) + "/work"
COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def rc(s): return s.translate(COMP)[::-1]

def load_ref():
    parts = [l.strip() for l in open(REF) if not l.startswith(">")]
    return "".join(parts)

def write_fa(path, name, s):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(s), 80):
            f.write(s[i:i+80] + "\n")

def mutate(s, rate, rng):
    """Uniform substitutions at the given per-site rate."""
    if rate <= 0:
        return s
    n = len(s)
    k = int(n * rate)
    b = bytearray(s, "ascii")
    alt = {ord('A'): b"CGT", ord('C'): b"AGT", ord('G'): b"ACT", ord('T'): b"ACG"}
    for _ in range(k):
        i = rng.randrange(n)
        o = alt.get(b[i])
        if o:
            b[i] = o[rng.randrange(3)]
    return b.decode("ascii")

def place(n_events, length, L, rng, margin=60000, spacing=50000):
    """n_events non-overlapping [start, end) intervals of exactly `length`."""
    need = n_events * (length + spacing)
    if need > L - 2 * margin:
        return None
    slots = []
    cursor = margin
    room = (L - 2*margin - n_events*length - (n_events-1)*spacing)
    gaps = sorted(rng.randrange(room + 1) for _ in range(n_events)) if room > 0 else [0]*n_events
    for k in range(n_events):
        st = cursor + gaps[k] + k*(length + spacing) - (gaps[k-1] if k else 0)*0
        slots.append((st, st + length))
        cursor = margin
    # recompute cleanly: cumulative layout with random slack
    slots = []
    pos = margin
    slack = room
    for k in range(n_events):
        take = rng.randrange(slack + 1) if slack > 0 else 0
        pos += take
        slack -= take
        slots.append((pos, pos + length))
        pos += length + spacing
    return slots

def apply_inversions(s, intervals):
    b = list(s)
    for st, en in intervals:
        b[st:en] = list(rc("".join(b[st:en])))
    return "".join(b)

def apply_translocations(s, intervals, rng):
    """Excise every segment and reinsert them elsewhere, orientation preserved.

    Done as one block permutation rather than a sequence of edits: applying
    translocations one at a time shifts the coordinates of every interval to the
    right of each insertion, so later excisions cut the wrong bases. That bug
    made detection power *fall* with event size, which is the opposite of what a
    resolution limit looks like.
    """
    iv = sorted(intervals)
    segs, rest, keep_map = [], [], []
    cur = 0
    for st, en in iv:
        rest.append(s[cur:st])
        keep_map.append((cur, st))
        segs.append(s[st:en])
        cur = en
    rest.append(s[cur:])
    keep_map.append((cur, len(s)))
    rest_s = "".join(rest)

    # Insertion points inside the retained sequence, kept clear of each other
    # and of the seams left by the excisions.
    seams = []
    acc = 0
    for a, b in keep_map:
        acc += b - a
        seams.append(acc)
    guard = max(20000, 2 * (iv[0][1] - iv[0][0]))
    picks = []
    tries = 0
    while len(picks) < len(segs) and tries < 200000:
        tries += 1
        c = rng.randrange(guard, len(rest_s) - guard)
        if all(abs(c - x) > guard for x in seams) and all(abs(c - p) > guard for p in picks):
            picks.append(c)
    picks.sort()
    if len(picks) < len(segs):
        segs = segs[:len(picks)]
        iv = iv[:len(picks)]

    order = list(range(len(picks)))
    rng.shuffle(order)
    out, prev = [], 0
    for k, c in enumerate(picks):
        out.append(rest_s[prev:c])
        out.append(segs[order[k]])
        prev = c
    out.append(rest_s[prev:])

    # An insertion at offset c in the retained sequence sits between two
    # reference landmarks; map it back so the acceptor junction is scored as
    # truth rather than counted against the method as a false positive.
    def to_ref(c):
        acc = 0
        for a, b in keep_map:
            if acc + (b - a) >= c:
                return a + (c - acc)
            acc += b - a
        return keep_map[-1][1]

    return "".join(out), [(st, en) for st, en in iv], [to_ref(c) for c in picks]

def digest(fa, tgt, panel):
    subprocess.run([BIN, "digest", "-i", fa, "-o", tgt, "-e", panel],
                   check=True, capture_output=True)

def synteny(tgt_a, tgt_b, out):
    pair = out.replace(".csv", ".pair.tgt")   # the CLI dispatches on extension
    with open(pair, "w") as f:
        f.write(open(tgt_a).read())
        f.write(open(tgt_b).read())
    subprocess.run([BIN, "synteny", "-i", pair, "-o", out],
                   check=True, capture_output=True)
    rows = [r for r in csv.DictReader(l for l in open(out) if not l.startswith("#"))]
    juncs = []
    jp = out.replace(".csv", "") + ".junctions.tsv"
    if os.path.exists(jp):
        for line in open(jp).readlines()[1:]:
            p = line.split("\t")
            if len(p) >= 3:
                juncs.append(int(p[2]))
    return rows[0], sorted(juncs)

def detected(truth_bounds, juncs, spacing):
    """Greedy nearest matching of junctions to true boundaries.

    A junction is reported at the position of the *left* landmark of the broken
    adjacency, so it lands somewhere in [boundary - gap, boundary] where the gap
    is exponentially distributed with mean equal to the landmark spacing. Any
    fixed window therefore throws away a predictable share of true detections
    (a 3x-spacing window misses e^-3 = 5% of them), which reads as a method
    failure and is not one. Matching to the nearest unclaimed boundary removes
    the window from the answer; the cap only guards against absurd pairings.

    Returns (events detected on both boundaries, junctions matching nothing).
    """
    cap = 20 * spacing
    bounds = []
    for i, (st, en) in enumerate(truth_bounds):
        bounds.append((st, i)); bounds.append((en, i))
    pairs = sorted(
        ((abs(j - b), j, b, i) for j in juncs for b, i in bounds if abs(j - b) <= cap)
    )
    used_j, used_b, hits = set(), set(), {}
    for _, j, b, i in pairs:
        if j in used_j or (b, i) in used_b:
            continue
        used_j.add(j); used_b.add((b, i))
        hits[i] = hits.get(i, 0) + 1
    both = sum(1 for i in hits if hits[i] == 2)
    return both, len(juncs) - len(used_j)

def main():
    rng = random.Random(20260829)
    ref = load_ref(); L = len(ref)
    panels = {"BcgI": "BcgI", "panel4": "BcgI,AlfI,AloI,FalI"}
    write_fa(f"{WORK}/ref.fasta", "ref", ref)
    ref_tgt = {}
    for pname, spec in panels.items():
        ref_tgt[pname] = f"{WORK}/ref.{pname}.tgt"
        digest(f"{WORK}/ref.fasta", ref_tgt[pname], spec)

    LENGTHS = [500, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000, 256000]
    TARGET  = 40          # events accumulated per cell
    results = []
    out = open(f"{WORK}/../power.tsv", "w")
    out.write("type\tpanel\tsub\tlength\tevents\tdetected\tpower\t"
              "junctions\texpected_junctions\tfalse_pos\ttrue_frac\tobs_frac\tshared\n")

    def run_cell(etype, pname, sub, length, seedbase):
        spec = panels[pname]
        got_events = got_det = got_j = got_fp = 0
        exp_j = 0
        tf_sum = of_sum = 0.0; n_gen = 0; shared = 0
        seed = seedbase
        while got_events < TARGET:
            rng2 = random.Random(seed); seed += 1
            cap = max(1, int(0.30 * L / length))      # stay clear of the 0.5 saturation
            per = min(TARGET - got_events, 40, cap)
            iv = place(per, length, L, rng2)
            while iv is None and per > 1:
                per = per // 2
                iv = place(per, length, L, rng2)
            if iv is None:
                break
            if etype == "inversion":
                g = apply_inversions(ref, iv)
                bounds = list(iv)
                extra = []
                exp_per = 2
            else:
                g, moved, acceptors = apply_translocations(ref, iv, rng2)
                bounds = list(moved)
                # the acceptor site is the translocation's third junction; pair
                # it with itself so the matcher can claim it
                extra = [(a, a) for a in acceptors]
                exp_per = 3
            g = mutate(g, sub, rng2)
            tag = f"{etype}.{pname}.s{sub}.L{length}.{seed}"
            fa  = f"{WORK}/{tag}.fasta"; tg = f"{WORK}/{tag}.tgt"
            write_fa(fa, tag.replace(".", "_"), g)
            digest(fa, tg, spec)
            row, juncs = synteny(ref_tgt[pname], tg, f"{WORK}/{tag}.csv")
            # A junction can only land on a surviving landmark, so the window
            # has to scale with landmark spacing: at 5% divergence BcgI keeps
            # 18% of its tags and the nearest landmark to a true boundary can be
            # tens of kb away. A fixed window scores that as a miss and reads as
            # a method failure when it is a measurement artifact.
            spacing = L / max(1, int(row["shared_tags"]))
            d, fp = detected(bounds + extra, juncs, spacing)
            got_det += d
            got_fp += fp
            got_events += len(bounds)
            got_j += int(row["breakpoints"]) if row["breakpoints"] != "NA" else 0
            exp_j += per * exp_per
            if etype == "inversion":
                tf_sum += per * length / L
                of_sum += float(row["inverted_fraction"])
                n_gen += 1
            shared = int(row["shared_tags"])
            for f in (fa, tg, f"{WORK}/{tag}.csv", f"{WORK}/{tag}.pair.tgt",
                      f"{WORK}/{tag}.junctions.tsv"):
                if os.path.exists(f): os.remove(f)
        power = got_det / got_events if got_events else 0.0
        line = (f"{etype}\t{pname}\t{sub}\t{length}\t{got_events}\t{got_det}\t{power:.4f}\t"
                f"{got_j}\t{exp_j}\t{got_fp}\t{tf_sum:.5f}\t{of_sum:.5f}\t{shared}\n")
        out.write(line); out.flush()
        print(line.strip(), file=sys.stderr)

    # Pass A: power vs length, no substitutions
    for etype in ("inversion", "translocation"):
        for pname in ("BcgI", "panel4"):
            for length in LENGTHS:
                run_cell(etype, pname, 0.0, length, hash((etype, pname, length)) % 10**6)

    # Pass B: substitution robustness at four lengths
    for pname in ("BcgI", "panel4"):
        for sub in (0.005, 0.01, 0.02, 0.05):
            for length in (2000, 8000, 32000, 128000):
                run_cell("inversion", pname, sub, length, hash((pname, sub, length)) % 10**6)

    out.close()

main()
