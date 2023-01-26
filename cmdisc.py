# Sage implementations of Algorithms 1 and 2

from sage.libs.pari.convert_sage import gen_to_sage
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial

def AlgorithmOne(H, check_monic_irreducible=False, verbose=False):
    """ Returns True if H is a Hilbert class polynomial, False otherwise. """

    # optional check that input is monic and irreducible
    if check_monic_irreducible:
        if not f.is_monic() or not f.is_irreducible():
            return False

    h = H.degree()
    rh = RR(h)
    lh = ((2.15*rh.log() + 4*((rh+2).log().log()+2).log() + 26.19)**2).floor()
    T = polygen(ZZ)
    X,Y = polygens(ZZ,'XY',2)
    H = H.change_ring(ZZ)(T)
    Fq = GF(101)
    Tq = polygen(Fq);
    Xq,Yq = polygens(Fq,['Xq', 'Yq'],2)
    Hmodq = H.change_ring(Fq)
    HYq = H(Yq)
    l = ZZ(2)
    while True:
        l = l.next_prime()
        if l > lh:
            print("Exceeded GRH bound l(h)={}, terminating".format(lh))
            return 0
        Hl = H.change_ring(GF(l))
        if not Hl.is_squarefree():
            continue
        if supersingular_j_polynomial(l).divides(Hl):
            continue
        if verbose:
            print("Using l={}".format(l))
        phi = gen_to_sage(pari.polmodular(l),{'x':X,'y':Y})
        res = phi([Xq,Yq]).resultant(HYq,Yq)([Tq,0])
        quo, rem = res.quo_rem(Hmodq)
        if rem or not Hmodq.divides(quo):
            return False
        res = phi.resultant(H(Y),Y)([T,0])
        quo, rem = res.quo_rem(H)
        if rem or not H.divides(quo):
            return False
        return True

def AlgorithmTwo(H, check_monic_irreducible=False, verbose=False):
    """ Returns D if H is the Hilbert class polynomial H_D, 0 otherwise. """

    # optional check that input is monic and irreducible
    if check_monic_irreducible:
        if not H.is_monic() or not H.is_irreducible():
            return False

    h = H.degree()
    h2list = [d for d in h.divisors() if (d-h)%2==0 and d.prime_to_m_part(2)==1]
    pmin = 33 * (h**2 * (RR(h+2).log().log()+2)**2).ceil()
    # Guarantees 4*p > |D| for fundamental D under GRH
    p = pmin-1
    n = 0
    while True:
        p = next_prime(p)  # Remark 11 is not applied here (this is asymptotically suboptimal)
        n += 1
        Hp = H.change_ring(GF(p))
        z = Hp.parent().quotient(Hp).gen()
        r = z**p-z
        d = r.lift().gcd(Hp).degree()  # number of roots mod p
        #assert d==len(Hp.roots())
        if d==0:
            continue
        if not Hp.is_squarefree():
            continue
        if d<h and d not in h2list:
            return 0
        jp = Hp.any_root(degree=-1, assume_squarefree=True)
        E = EllipticCurve(j=jp)
        if E.is_supersingular():
            continue
        if verbose:
            print("Found a suitable prime p = {} after testing {} primes\n".format(p, n))
        D = E.frobenius_polynomial().discriminant()
        D0 = D.squarefree_part()
        if D0%4 in [2,3]:
            D0 *= 4
        f = ZZ(D//D0).isqrt()
        for d in f.divisors():
            D1 = D0*d**2
            if h != D1.class_number():
                continue
            hilbert_class_polynomial.clear_cache()
            if H == hilbert_class_polynomial(D1):
                return D1
        return 0

D1000 = [-d for d in srange(3,1001) if d%4 in [0,3]]
assert len(D1000) == 500

AbelianD1 = [-3, -4, -7, -8, -11, -12, -16, -19, -27, -28, -43, -67, -163]
AbelianD2 = [ -15, -20, -24, -32, -35, -36, -40, -48, -51, -52, -60, -64, -72, -75, -88, -91, -99, -100, -112, -115, -123, -147, -148, -187, -232, -235, -267, -403, -427]
AbelianD4 = [ -84, -96, -120, -132, -160, -168, -180, -192, -195, -228, -240, -280, -288, -312, -315, -340, -352, -372, -408, -435, -448, -483, -520, -532, -555, -595, -627, -708, -715, -760, -795, -928, -1012, -1435 ]
AbelianD8 = [ -420, -480, -660, -672, -840, -960, -1092, -1120, -1155, -1248, -1320, -1380, -1428, -1540, -1632, -1848, -1995, -2080, -3003, -3040, -3315 ]
AbelianD16 = [ -3360, -5280, -5460, -7392 ]
AbelianD = AbelianD1 + AbelianD2 + AbelianD4 + AbelianD8 + AbelianD16
assert len(AbelianD) == 101

# [least,median,largest] discriminants of class number 5,10,15,...,100
hD100 = [[-47,-571,-2683],[-119,-2299,-13843],[-239,-6571,-34483],[-455,-9124,-58843],[-479,-25747,-93307],[-671,-21592,-137083],[-1031,-42499,-210907],[-1271,-34180,-274003],
         [-1319,-60748,-308323],[-1799,-64203,-389467],[-4463,-101419,-452083],[-2159,-83176,-662803],[-3527,-138883,-703123],[-3239,-138979,-821683],[-4703,-157051,-916507],
         [-5183,-133620,-1165483],[-4079,-252988,-1285747],[-5951,-204619,-1548523],[-6959,-251443,-1659067],[-7991,-249451,-1856563]]
assert len(hD100) == 20

# [least,median,largest] fundamental discriminants of class number 125,150,175,...,1000
hD1000 = [[-11519,-570139,-2944363],[-11879,-658891,-4163443],[-18719,-1124131,-4742467],[-21311,-984952,-7206763],[-25031,-1911787,-9119203],[-30551,-1881987,-12208267],[-38351,-2682227,-14420323],
           [-34271,-2497819,-19184323],[-38639,-3414379,-19385683],[-47759,-3724003,-24669283],[-53231,-5552779,-26809843],[-67031,-3890859,-32492923],[-84719,-6756619,-35307787],
           [-75599,-6503512,-40868683],[-99839,-8428291,-43078963],[-96599,-6895428,-53698747],[-117911,-9473059,-51065683],[-127151,-9325347,-69758203],[-136751,-11505827,-61595227],
           [-148511,-9248548,-72789403],[-234599,-13993459,-81289723],[-227015,-13278259,-85527187],[-269879,-17452579,-104367187],[-185471,-13972667,-96016243],[-226799,-19566907,-106911667],
           [-218951,-17530243,-113620123],[-265079,-21602923,-134009467],[-233999,-15304228,-129667987],[-316391,-24836347,-130788067],[-266279,-23348968,-135386227],[-250799,-27051163,-149457163],
           [-299519,-22973988,-180619363],[-314159,-30598123,-164446027],[-351911,-28525507,-186089443],[-376631,-35460091,-184039507],[-412079,-25671795,-220287043]]
assert len(hD1000) == 36

def CMProfile(alg, discs, exact=True, detail=1):
    start = cputime()
    x = polygen(ZZ)

    if detail == 3:
        print("\\begin{tabular}{rrrrrrr}\n$h$ & $|H|$ & $D$ & $t_{\\rm HCP}$ & $t_{\\rm CM}$ & $t_{\\rm noCM}$\\\\\\toprule")
    else:
        print("Computing H_D and checking H_D(x),H_D(x)+1) for {} discriminants".format(len(discs)))

    for D in discs:
        t0 = cputime()
        hilbert_class_polynomial.clear_cache()
        H = hilbert_class_polynomial(D)
        t0 = cputime()-t0
        t1 = cputime()
        assert (alg(H) == D if exact else alg(H) != 0)
        t1 = cputime()-t1
        # t2 = cputime()
        # assert alg(H(x+1)) == 0
        # t2 = cputime()-t2
        t3 = cputime()
        assert alg(H+1) == 0
        t3 = cputime()-t3
        h = H.degree()
        ht = max(abs(c) for c in H).log().round()
        if detail==1:
            print(".", end="")
            sys.stdout.flush()
        elif detail==2:
            print("h={}, |H|={}, D={}: {:.2f}s, {:.2f}s, {:.2f}s".format(h,ht,D,t0,t1,t3));
        elif detail==3:
            print("{} & {} & {} & {:.2f} & {:.2f} & {:.2f}\\\\".format(h,ht,D,t0,t1,t3))

    if detail == 1:
        print()
    if detail == 3:
        print("\\bottomrule\n\\end{tabular}")

    print("Verified {} discriminants in {:.3f}s".format(len(discs), cputime()-start))

# display usage instruction on loading this file

print("Example usage:")
print("CMProfile(AlgorithmTwo,[r[2] for r in hD100], detail=3)")
print("or")
print("CMProfile(AlgorithmOne,[r[2] for r in hD100], exact=False, detail=2)")
