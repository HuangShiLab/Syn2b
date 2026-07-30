#!/usr/bin/env bash
# Prove that the structural synteny metric separates substitution load from
# structural change, and that the previous metric did not.
#
# The property being tested is an invariance, not an accuracy: substitutions
# change tag PRESENCE, structure changes tag ORDER, and a structural metric must
# respond only to the second. So the substitution row must be flat at 1.0 with
# zero breakpoints, at every divergence, while the inversion row must show a
# constant response.
#
# Every genome here is derived from one reference by counting out substitutions,
# so the divergence is exact by construction, and the inversion is placed by us,
# so its presence is known rather than inferred.
#
# Requires: python3 + numpy, and a built syn2b. Downloads ~5 MB.
set -euo pipefail

OUT="${1:-invariance}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${SYN2B:-$REPO/target/release/syn2b}"
ENZYME="${ENZYME:-BcgI}"

if [ ! -x "$BIN" ]; then
    echo "error: no binary at $BIN — run: (cd $REPO && cargo build --release)" >&2
    exit 1
fi
newer=$(find "$REPO/src" "$REPO/Cargo.toml" -newer "$BIN" -print -quit 2>/dev/null || true)
if [ -n "$newer" ]; then
    echo "error: $BIN is older than the source tree ($newer)." >&2
    echo "       run: (cd $REPO && cargo build --release)" >&2
    exit 1
fi

mkdir -p "$OUT"
cd "$OUT"

# E. coli K-12 MG1655, ENA U00096.3. A pinned public accession so the numbers
# below are reproducible rather than approximately reproducible.
if [ ! -s ref.fa ]; then
    curl -sSL "https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.3?download=true" -o ref.fa
    bases=$(grep -v '^>' ref.fa | tr -d '\n\r ' | wc -c | tr -d ' ')
    if [ "$bases" -ne 4641652 ]; then
        echo "error: got $bases bases, expected 4641652 — incomplete download" >&2
        exit 1
    fi
fi

python3 - <<'PY'
import numpy as np

LEVELS = [1.0, 0.999, 0.995, 0.99, 0.98, 0.95]
INVERSION_BP = 400_000
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)

seq = []
for line in open("ref.fa"):
    if not line.startswith(">"):
        seq.append(line.strip())
seq = "".join(seq).upper()
seq = "".join(c for c in seq if c in "ACGT")
u8 = np.frombuffer(seq.encode(), dtype=np.uint8).copy()

lut = np.arange(256, dtype=np.uint8)
for a, b in zip(b"ACGT", b"TGCA"):
    lut[a] = b


def write(path, name, arr):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        s = arr.tobytes().decode()
        for i in range(0, len(s), 80):
            fh.write(s[i : i + 80] + "\n")


write("ref_clean.fa", "ref", u8)

for ani in LEVELS:
    for invert in (False, True):
        # Same seed for both arms so the substitution pattern is identical and
        # the only difference is the inversion.
        rng = np.random.default_rng(4242)
        out = u8.copy()
        n_sub = int(round((1.0 - ani) * out.size))
        if n_sub:
            pos = rng.choice(out.size, size=n_sub, replace=False)
            idx = np.searchsorted(BASES, out[pos])
            idx = np.where(idx >= 4, 0, idx)
            out[pos] = BASES[(idx + rng.integers(1, 4, size=n_sub, dtype=np.int64)) % 4]
        if invert:
            lo = out.size // 3
            hi = min(lo + INVERSION_BP, out.size)
            out[lo:hi] = lut[out[lo:hi][::-1]]
        tag = f"{'inv' if invert else 'sub'}_{ani:.4f}"
        write(f"{tag}.fa", tag, out)
print(f"wrote {2 * len(LEVELS)} genomes: substitutions only, and the same "
      f"substitutions plus a {INVERSION_BP // 1000} kb inversion")
PY

"$BIN" digest -i ref_clean.fa -o ref.tgt -e "$ENZYME" >/dev/null

printf "\n%-10s %-6s %11s %19s %11s %8s %s\n" \
    popANI arm structural breakpoint_density legacy breakpts shared
printf '%s\n' "------------------------------------------------------------------------------------"
for a in 1.0000 0.9990 0.9950 0.9900 0.9800 0.9500; do
    for k in sub inv; do
        "$BIN" digest -i "${k}_${a}.fa" -o q.tgt -e "$ENZYME" >/dev/null
        cat ref.tgt q.tgt > pair.tgt
        "$BIN" synteny -i pair.tgt -o pair.csv >/dev/null
        row=$(grep -v '^#' pair.csv | grep -v '^genome_A' | head -1)
        printf "%-10s %-6s %11s %19s %11s %8s %s\n" \
            "$(python3 -c "print(f'{$a*100:.2f}%')")" "$k" \
            "$(echo "$row" | cut -d, -f3)" "$(echo "$row" | cut -d, -f5)" \
            "$(echo "$row" | cut -d, -f8)" "$(echo "$row" | cut -d, -f4)" \
            "$(echo "$row" | cut -d, -f6)"
    done
done

cat <<'EOF'

Expected, and what each column proves:

  structural         1.0000 on every `sub` row. Substitutions remove tags but do
                     not reorder the survivors, so a structural metric must not
                     move. The `inv` rows sit near 0.82 regardless of divergence.
  breakpoint_density 0.000 on every `sub` row, ~0.19-0.20 on every `inv` row —
                     roughly a 50x separation, and normalised by shared-tag count
                     so it stays comparable as tags are lost.
  legacy             The previous metric, for contrast: it falls from 1.0000 to
                     ~0.011 across the `sub` rows with no structural variation
                     present at all, because it is the tag-survival curve
                     (0.99^32 = 0.725 per tag, squared for an adjacency).
EOF
