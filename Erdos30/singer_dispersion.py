"""
Singer Difference Set Channel Dispersion Computation
Erdős Problem #30 — Sidon Set Density

Computes V_Singer(q) for the first 20 prime powers.

CONSTRUCTION (Singer 1938):
  For a prime power q = p^k, build GF(p^{3k}).
  Let α be a primitive element of GF(p^{3k}).
  Set N = q²+q+1.
  The Singer difference set is:
      A = { n ∈ {0,...,N-1} : Tr_{GF(q^3)/GF(q)}(α^n) = 0 }
  This gives a (N, q+1, 1)-perfect difference set in Z_N.

CHANNEL MODEL:
  Sidon MAC: (a,b) → a+b mod N, where (a,b) ~ Unif over symmetric pairs {(a,b): a ≤ b, a,b ∈ A}.

  For unordered symmetric pairs:
    p(s | a,b) = 1  [deterministic channel]
    p(s) = r_sym(s) / |{(a,b): a ≤ b}|
    where r_sym(s) = |{(a,b): a ≤ b, a+b ≡ s (mod N), a,b ∈ A}|

  Information density: i(a,b;s) = log(p(s|a,b)/p(s)) = log(M/r_sym(s))
  where M = |{(a,b): a ≤ b}| = |A|(|A|+1)/2.

  Channel dispersion: V = Var_{(a,b)~Unif(sym pairs)} [ i(a,b; a+b) ]
  V = 0  iff  r_sym(s) is constant over all s in the sumset.

  For a Sidon set (all sums a+b with a < b distinct, and 2a distinct from all a+b):
    r_sym(s) ∈ {0, 1} for all s → r_sym is constant (=1) on the sumset.
    → i(a,b;s) = log(M) is constant → V = 0.
"""

import math
import sys
from collections import Counter
from itertools import product as iproduct


# ─── Primality and prime-power detection ──────────────────────────────────────

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True

def prime_power_base(q):
    """Return (p, k) s.t. q = p^k (prime p), or None."""
    if q < 2: return None
    for p in range(2, q + 1):
        if not is_prime(p): continue
        k, pk = 0, 1
        while pk <= q:
            pk *= p; k += 1
        pk //= p; k -= 1
        if pk == q:
            return (p, k)
    return None


# ─── GF(p^n) polynomial arithmetic ───────────────────────────────────────────

def _poly_divmod(a, b, p):
    """Polynomial division: returns (q, r) s.t. a = b*q + r over GF(p).
    Polynomials stored as lists, constant term first."""
    a = list(a)
    b = list(b)
    while b and b[-1] == 0: b.pop()
    deg_b = len(b) - 1
    while a and a[-1] == 0: a.pop()
    if len(a) < len(b):
        return [0], a
    inv_lead = pow(b[-1], p - 2, p)  # modular inverse (b leading coeff)
    quotient = []
    while len(a) >= len(b):
        coef = a[-1] * inv_lead % p
        quotient.append(coef)
        for i in range(len(b)):
            a[len(a) - len(b) + i] = (a[len(a) - len(b) + i] - coef * b[i]) % p
        a.pop()
    quotient.reverse()
    while a and a[-1] == 0: a.pop()
    return quotient, a if a else [0]


def _poly_mul(a, b, p):
    result = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] = (result[i + j] + ai * bj) % p
    return result


def _poly_mod(a, b, p):
    return _poly_divmod(a, b, p)[1]


def _poly_pow_mod(base, exp, mod, p):
    result = [1]
    b = list(base)
    while exp > 0:
        if exp & 1:
            result = _poly_mod(_poly_mul(result, b, p), mod, p)
        b = _poly_mod(_poly_mul(b, b, p), mod, p)
        exp >>= 1
    return result


def _poly_gcd(a, b, p):
    a, b = list(a), list(b)
    while True:
        while b and b[-1] == 0: b.pop()
        if not b or (len(b) == 1 and b[0] == 0): break
        _, a = _poly_divmod(a, b, p)
        a, b = b, a
    while a and a[-1] == 0: a.pop()
    if not a: return [0]
    inv = pow(a[-1], p - 2, p)
    return [(c * inv) % p for c in a]


def find_irreducible(p, n):
    """
    Find a monic irreducible polynomial of degree n over GF(p).
    Uses Rabin's irreducibility test.
    Returns coefficient list (constant term first), length n+1.
    """
    def is_irred(f):
        # Rabin test: check for each prime r|n that gcd(f, x^(p^(n/r)) - x) = 1
        # and x^(p^n) ≡ x (mod f)
        deg_f = len(f) - 1
        prime_factors = set()
        tmp = n
        for d in range(2, tmp + 1):
            if tmp % d == 0:
                prime_factors.add(d)
                while tmp % d == 0: tmp //= d
        if tmp > 1: prime_factors.add(tmp)

        for r in prime_factors:
            exp = p ** (n // r)
            xpow = _poly_pow_mod([0, 1], exp, f, p)
            # subtract x
            diff = list(xpow) + [0] * max(0, 2 - len(xpow))
            while len(diff) < 2: diff.append(0)
            diff[1] = (diff[1] - 1) % p
            while diff and diff[-1] == 0: diff.pop()
            if not diff: diff = [0]
            g = _poly_gcd(f, diff, p)
            if len(g) > 1:  # degree ≥ 1 → not irreducible
                return False

        # Check x^(p^n) ≡ x (mod f)
        xpow = _poly_pow_mod([0, 1], p ** n, f, p)
        diff = list(xpow) + [0] * max(0, 2 - len(xpow))
        while len(diff) < 2: diff.append(0)
        diff[1] = (diff[1] - 1) % p
        rem = _poly_mod(diff, f, p)
        while rem and rem[-1] == 0: rem.pop()
        return not rem or rem == [0]

    # Enumerate monic degree-n polynomials (constant term ≠ 0 for n > 1)
    for coeffs in iproduct(range(p), repeat=n):
        if n > 1 and coeffs[0] == 0:
            continue  # x divides polynomial → reducible
        f = list(coeffs) + [1]  # monic, degree n
        if is_irred(f):
            return f
    raise ValueError(f"No monic irreducible of degree {n} over GF({p}) found")


class GFpn:
    """
    GF(p^n) via Z_p[x] / (irred poly of degree n).
    Elements are tuples of n coefficients (constant term first).
    """

    def __init__(self, p, n, irred=None):
        self.p = p
        self.n = n
        if irred is None:
            irred = find_irreducible(p, n)
        self.irred = list(irred)

    def _reduce(self, poly):
        """Reduce to canonical form: tuple of n coefficients."""
        p, n = self.p, self.n
        a = list(poly)
        irred = self.irred
        # Remove leading zeros
        while len(a) > n:
            while a and a[-1] == 0:
                a.pop()
            if len(a) <= n: break
            coef = a[-1]
            shift = len(a) - 1 - n
            for i, c in enumerate(irred):
                a[shift + i] = (a[shift + i] - coef * c) % p
            a.pop()
        while len(a) < n: a.append(0)
        return tuple(a[:n])

    def zero(self): return tuple([0] * self.n)
    def one(self): r = [0]*self.n; r[0] = 1; return tuple(r)

    def add(self, a, b): return tuple((x+y)%self.p for x,y in zip(a,b))
    def sub(self, a, b): return tuple((x-y)%self.p for x,y in zip(a,b))

    def mul(self, a, b):
        p, n = self.p, self.n
        result = [0] * (2*n - 1)
        for i, ai in enumerate(a):
            for j, bj in enumerate(b):
                result[i+j] = (result[i+j] + ai*bj) % p
        return self._reduce(result)

    def pow(self, a, exp):
        if exp == 0: return self.one()
        result = self.one()
        base = a
        while exp > 0:
            if exp & 1: result = self.mul(result, base)
            base = self.mul(base, base)
            exp >>= 1
        return result

    def all_elements(self):
        for c in iproduct(range(self.p), repeat=self.n):
            yield c

    def find_primitive(self):
        """Find a primitive element (generator of GF(p^n)*)."""
        order = self.p**self.n - 1
        factors = set()
        tmp = order
        for d in range(2, int(tmp**0.5)+1):
            if tmp % d == 0:
                factors.add(d)
                while tmp % d == 0: tmp //= d
        if tmp > 1: factors.add(tmp)

        for elem in self.all_elements():
            if elem == self.zero(): continue
            ok = True
            for f in factors:
                if self.pow(elem, order//f) == self.one():
                    ok = False; break
            if ok: return elem
        return None


# ─── Trace map ────────────────────────────────────────────────────────────────

def trace_gfq3_over_gfq(F, x, q):
    """
    Absolute trace Tr_{GF(q^3)/GF(q)}(x) = x + x^q + x^(q^2).
    F is GF(p^{3k}) where q = p^k.
    Returns an element of F (which should lie in GF(q) ⊂ F).
    """
    xq  = F.pow(x, q)
    xq2 = F.pow(xq, q)
    return F.add(F.add(x, xq), xq2)


# ─── Singer set construction ──────────────────────────────────────────────────

def construct_singer_set(q):
    """
    Construct the Singer (q²+q+1, q+1, 1)-difference set.
    Returns (A, N, error_msg).

    Algorithm:
      1. Build GF(p^{3k}) where q = p^k.
      2. Find primitive element α of GF(p^{3k}).
      3. Singer set = { n ∈ {0,...,N-1} : Tr_{q^3/q}(α^n) = 0 }
         Note: index n ∈ Z_N = Z_{q²+q+1}, NOT Z_{q^3-1}.
         The trace condition is periodic mod N since α^N ∈ GF(q)* and
         the trace is GF(q)-linear.
    """
    pp = prime_power_base(q)
    if pp is None:
        return None, q*q+q+1, f"q={q} is not a prime power"
    p, k = pp
    N = q*q + q + 1

    n3 = 3*k  # GF(q^3) = GF(p^{3k})
    try:
        F = GFpn(p, n3)
    except Exception as e:
        return None, N, f"GF({p}^{n3}) failed: {e}"

    alpha = F.find_primitive()
    if alpha is None:
        return None, N, "No primitive element found"

    # Precompute α^n for n = 0, ..., N-1
    alpha_pows = [F.one()]
    cur = F.one()
    for _ in range(N - 1):
        cur = F.mul(cur, alpha)
        alpha_pows.append(cur)

    zero = F.zero()
    A = [n for n, a in enumerate(alpha_pows)
         if trace_gfq3_over_gfq(F, a, q) == zero]

    if len(A) != q + 1:
        return None, N, f"Trace kernel has {len(A)} elements, expected q+1={q+1}"

    return A, N, None


# ─── Difference set verification ─────────────────────────────────────────────

def verify_perfect_diff_set(A, N):
    """
    Check that A is a (N, |A|, 1)-perfect difference set:
    every nonzero element of Z_N appears exactly once as a-b (a,b ∈ A, a≠b).
    Returns (is_perfect, n_diffs_covered, max_mult).
    """
    diff_count = Counter()
    for a in A:
        for b in A:
            if a != b:
                diff_count[(a - b) % N] += 1
    vals = list(diff_count.values())
    if not vals:
        return False, 0, 0
    return (all(v == 1 for v in vals) and len(diff_count) == N - 1), len(diff_count), max(vals)


# ─── Channel dispersion ───────────────────────────────────────────────────────

def compute_dispersion(A, N):
    """
    Compute channel dispersion V for the Sidon MAC (a,b) → a+b mod N.

    Input distribution: SYMMETRIC (unordered) pairs (a,b) with a ≤ b, a,b ∈ A.
    This matches the task specification "a ≤ b" in the pair sampling.

    Channel: deterministic, p(s | a, b) = 1 if s = (a+b) mod N, else 0.

    Output distribution: p(s) = r_sym(s) / M
      where r_sym(s) = |{(a,b): a≤b, a,b∈A, (a+b)%N = s}|
      and M = total number of symmetric pairs = |A|(|A|+1)/2.

    Information density: i(a,b) = log(p(s|a,b)/p(s)) = log(M / r_sym(s))

    Dispersion: V = Var_{(a,b)~Unif(sym pairs)} [ i(a,b; a+b) ]
      V = 0  iff  r_sym(s) is constant on the sumset support.

    For a Sidon set:
      All sums a+b (a < b) are distinct, so r_sym(s) ≤ 1 for off-diagonal.
      Diagonal sums 2a ∈ A are also distinct (different a give different 2a mod N
      unless N is even and two elements sum to 0 mod N/2, but for Singer sets
      with N = q²+q+1 odd, this is fine).
      → r_sym(s) = 1 for all s in the sumset → V = 0.
    """
    size = len(A)
    M = size * (size + 1) // 2  # number of symmetric pairs

    r_sym = Counter()
    sym_pairs = []
    for i, a in enumerate(A):
        for b in A[i:]:  # a ≤ b (using sorted A)
            s = (a + b) % N
            r_sym[s] += 1
            sym_pairs.append((a, b, s))

    log_M = math.log(M)
    densities = [log_M - math.log(r_sym[s]) for _, _, s in sym_pairs]

    mean_d = sum(densities) / len(densities)
    variance = sum((d - mean_d)**2 for d in densities) / len(densities)

    max_r = max(r_sym.values())
    min_r = min(r_sym.values())
    sumset_size = len(r_sym)
    is_sidon = (max_r == 1)  # All representation numbers = 1

    return dict(
        variance=variance,
        mean_density=mean_d,
        sumset_size=sumset_size,
        max_r=max_r,
        min_r=min_r,
        is_sidon=is_sidon,
        M=M,
    )


# ─── Main ─────────────────────────────────────────────────────────────────────

PRIME_POWERS_20 = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41]


def main():
    print("=" * 110)
    print("Singer Difference Set Channel Dispersion  —  Erdős #30 (Sidon set density)")
    print("=" * 110)
    hdr = (f"{'q':>4}  {'N':>6}  {'|A|':>4}  {'|A+A|':>6}  "
           f"{'V_Singer':>14}  {'max(r_sym)':>10}  {'min(r_sym)':>10}  "
           f"{'PerfDS':>7}  {'Sidon':>6}")
    print(hdr)
    print("-" * 110)

    results = []
    for q in PRIME_POWERS_20:
        N_exp = q*q + q + 1
        print(f"  building q={q}...", end="", flush=True)
        A, N, err = construct_singer_set(q)
        if A is None:
            print(f"\r{q:>4}  {N_exp:>6}  {'—':>4}  {'—':>6}  {'—':>14}  {'—':>10}  {'—':>10}  {'—':>7}  {'—':>6}  SKIP: {err}")
            results.append({'q': q, 'N': N_exp, 'error': err})
            continue

        # A must be sorted for symmetric pair enumeration
        A.sort()

        is_perf, n_diffs, max_diff_r = verify_perfect_diff_set(A, N)
        stats = compute_dispersion(A, N)
        V = stats['variance']
        V_str = "0.000000e+00" if V < 1e-14 else f"{V:.6e}"
        sidon_str = "YES" if stats['is_sidon'] else "NO"
        perf_str = "YES" if is_perf else f"NO(max={max_diff_r})"

        print(f"\r{q:>4}  {N:>6}  {len(A):>4}  {stats['sumset_size']:>6}  "
              f"{V_str:>14}  {stats['max_r']:>10}  {stats['min_r']:>10}  "
              f"{perf_str:>7}  {sidon_str:>6}")

        results.append({
            'q': q, 'N': N, 'A_size': len(A),
            'sumset_size': stats['sumset_size'],
            'V': V,
            'max_r': stats['max_r'], 'min_r': stats['min_r'],
            'is_sidon': stats['is_sidon'],
            'is_perfect': is_perf,
            'M': stats['M'],
        })

    print("=" * 110)
    print()

    computed = [r for r in results if 'V' in r]
    skipped  = [r for r in results if 'error' in r]

    print("SUMMARY")
    print("-" * 70)
    print(f"  Sets computed   : {len(computed)}")
    print(f"  Sets skipped    : {len(skipped)}")
    if skipped:
        for r in skipped:
            print(f"    q={r['q']}: {r['error']}")
    print()

    if computed:
        vs = [r['V'] for r in computed]
        all_zero     = all(v < 1e-12 for v in vs)
        perfect_cnt  = sum(1 for r in computed if r.get('is_perfect', False))
        sidon_cnt    = sum(1 for r in computed if r.get('is_sidon', False))

        print("CHANNEL DISPERSION RESULTS")
        print("-" * 70)
        print(f"  Max  |V| = {max(vs):.6e}")
        print(f"  Mean |V| = {sum(vs)/len(vs):.6e}")
        print(f"  Min  |V| = {min(vs):.6e}")
        print()
        if all_zero:
            print("  *** ALL V = 0  ←  HYPOTHESIS CONFIRMED ***")
            print("  Singer difference sets are dispersion-free codes.")
        else:
            print("  [!] Some V > 0. Details:")
            for r in computed:
                if r['V'] > 1e-12:
                    print(f"    q={r['q']}: V={r['V']:.4e}, max_r={r['max_r']}, "
                          f"is_perfect={r['is_perfect']}, is_sidon={r['is_sidon']}")
        print()
        print(f"  Perfect difference sets : {perfect_cnt}/{len(computed)}")
        print(f"  Sidon (all r_sym = 1)   : {sidon_cnt}/{len(computed)}")

    print()
    print("THEORETICAL ANALYSIS")
    print("-" * 70)
    print("  Sidon MAC channel: (a,b) -> a+b mod N, (a,b) ~ Unif(symmetric pairs)")
    print()
    print("  For a Sidon set A in Z_N:")
    print("    All pairwise sums {a+b : a < b, a,b ∈ A} are distinct.")
    print("    Each diagonal sum {2a : a ∈ A} is also distinct (N odd).")
    print("    Therefore: r_sym(s) ∈ {0,1} for all s,")
    print("    and r_sym(s) = 1 for all s in the sumset A+A.")
    print()
    print("  Information density: i(a,b) = log(M / r_sym(a+b)) = log(M) = constant.")
    print("  ⟹  Var[i(a,b; a+b)] = 0  for all Sidon sets.  QED")
    print()
    print("  Singer difference sets are Sidon sets because they are perfect")
    print("  difference sets: each nonzero difference appears exactly once,")
    print("  which implies all sums are distinct (Sidon property follows from")
    print("  the perfect difference set property in Z_N with N odd).")


if __name__ == "__main__":
    main()
