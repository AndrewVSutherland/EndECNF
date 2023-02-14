# Wrapper around Sage's is_cm_j_invariant() function to test a (monic
# integer irreducible) polynomial for being an HCP and finding its
# discriminant.


def CMDiscriminant_orig(H):
    """
    Return the CM discriminant of an HCP, or 0, by calling Sage's
    is_cm_j_invariant() function directly, using the 'old' method
    which is described in section 3.2 of the paper (where it was
    called the 'first method'.
    """
    (flag, df) = is_cm_j_invariant(NumberField(H, 'j').gen(), method='old')
    if flag:
        d, f = df
        return d*f**2
    return 0

# For our profiling we instead use a precomputed table of all
# discriminants with class number h for all h<=1000, which we first
# read in.

def parse_one(L):
    """
    Parse one line from the data file.

    Each line has the form

    1:[3,4,7,8,11,12,16,19,27,28,43,67,163]

    giving a complete list of all discriminants for each class number.
    """
    h, Dlist = L.split(":")
    h = ZZ(h)
    Dlist = Dlist[1:-1]
    return h, [-ZZ(D) for D in Dlist.split(",")]

# Read the data file into a dict with keys class numbers h, values the
# list of discriminants of class number h

Ddict = dict(read_data("cmdiscs1000.txt", parse_one))

def CMDiscriminant(H):
    """
    Return the CM discriminant of an HCP, or 0, by comparing the input
    polynomial with `H_D` for each discriminant `D` of class number
    equal to the degree, in turn, as in Sage's is_cm_j_invariant()
    function using the 'old' method, but using the precomputed lsit of
    discriminants for each class number.
    """
    Dlist = Ddict[H.degree()]
    for D in Dlist:
        hilbert_class_polynomial.clear_cache()
        HD = hilbert_class_polynomial(D)
        if HD == H:
            return D
    return 0

AbelianD1 = [-3, -4, -7, -8, -11, -12, -16, -19, -27, -28, -43, -67, -163]
AbelianD2 = [ -15, -20, -24, -32, -35, -36, -40, -48, -51, -52, -60, -64, -72, -75, -88, -91, -99, -100, -112, -115, -123, -147, -148, -187, -232, -235, -267, -403, -427]
AbelianD4 = [ -84, -96, -120, -132, -160, -168, -180, -192, -195, -228, -240, -280, -288, -312, -315, -340, -352, -372, -408, -435, -448, -483, -520, -532, -555, -595, -627, -708, -715, -760, -795, -928, -1012, -1435 ]
AbelianD8 = [ -420, -480, -660, -672, -840, -960, -1092, -1120, -1155, -1248, -1320, -1380, -1428, -1540, -1632, -1848, -1995, -2080, -3003, -3040, -3315 ]
AbelianD16 = [ -3360, -5280, -5460, -7392 ]
AbelianD = AbelianD1 + AbelianD2 + AbelianD4 + AbelianD8 + AbelianD16
assert len(AbelianD) == 101

DD1000 = { i: [Ddict[i][0],Ddict[i][len(Ddict[i])//2],Ddict[i][-1]] for i in range(1,1001) }

def TestCMDiscriminant(B=1000):
    start = cputime()
    x = polygen(ZZ)
    DD = [-d for d in srange(3,B+1) if d%4 in [0,3]]
    print("Checking H_D(x),H_D(x+1),H_D(x)+1 for |D| <= {}".format(B))
    for D in DD:
        H = hilbert_class_polynomial(D)
        hilbert_class_polynomial.clear_cache()
        assert CMDiscriminant(H) == D
        hilbert_class_polynomial.clear_cache()
        assert CMDiscriminant(H(x+1)) == 0
        hilbert_class_polynomial.clear_cache()
        assert CMDiscriminant(H+1) == 0
        print(".", end="")
        sys.stdout.flush()
    t0 = cputime()
    print("\nVerified discriminants |D| <= {} in {:.3f}s".format(B, t0-start))
    print("Checking H_D(x),H_D(x+1),H_D(x)+1 for D with class group C2^n")
    for D in AbelianD:
        hilbert_class_polynomial.clear_cache()
        H = hilbert_class_polynomial(D)
        assert CMDiscriminant(H) == D
        assert CMDiscriminant(H(x+1)) == 0
        assert CMDiscriminant(H+1) == 0
    t1 = cputime()
    print("Verified abelian discriminants in {:.3f}s".format(t1-t0))

    print("Checking H_D(x),H_D(x)+1 for D of class number 5,10,15,20,25,50,75,100,...,1000.")
    DD = [DD1000[i][1] for i in [5,10,15,20,25,30,35,40,45,50,75]] + DD1000[100] + [DD1000[i][1] for i in range(125,1000,25)]
    for d in DD:
        D = ZZ(d)
        H = hilbert_class_polynomial(D)
        hilbert_class_polynomial.clear_cache()
        t1 = cputime()
        assert CMDiscriminant(H) == D
        t2 = cputime()
        t1 = t2-t1
        hilbert_class_polynomial.clear_cache()
        assert CMDiscriminant(H+1) == 0
        t2 = cputime()-t2
        print("h({})={}: {:.3f}s & {:.3f}s".format(D, D.class_number(),t1,t2))

    print("All tests completed successfully, total time {:.3f}s".format(cputime()-start))
