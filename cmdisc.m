Attach("gl2.m"); // to get the function EndomorphismRingData

// returns 1 if H=H_D for some D and 0 otherwise
function AlgorithmOne(H:Verbose:=false)
    // Uncomment the lines below for general use
    // if not IsMonic(H) then return false; end if;
    // if not IsIrreducible(H) then return false; end if;
    h := Degree(H);
    lh := Floor((2.15*Log(h) + 4*Log(Log(Log(h+2))+2) + 26.19)^2);
    ZZ := Integers();
    RT<T> := PolynomialRing(ZZ);
    RXY<X,Y> := PolynomialRing(ZZ,2);
    H := Evaluate(H,T);
    HY := Evaluate(H,Y);
    Fp := GF(101); // auxiliary finite field, per Remark 10
    RTp<Tp> := PolynomialRing(Fp);
    RXYp<Xp,Yp> := PolynomialRing(Fp,2);
    Hp := Evaluate(H,Tp);
    l := 2; // start with the first odd prime
    while true do
        l := NextPrime(l);
        if l gt lh then printf "Exceeded GRH bound l(h)=%o, terminating\n",lh; return 0; end if;
        Hl := ChangeRing(H,GF(l));
        if GCD(Hl,Derivative(Hl)) ne 1 then continue; end if;
        if IsDivisibleBy(SupersingularPolynomial(l),Hl) then continue; end if;
        if Verbose then print "Using l =",l; return 0; end if;
        phi := ClassicalModularPolynomial(l);
        res := Evaluate(Resultant(Evaluate(phi,[Xp,Yp]),Evaluate(Hp,Yp),Yp),[Tp,0]);
        if not IsDivisibleBy(res,Hp) then return 0; end if;
        res := ExactQuotient(res,Hp);
        if not IsDivisibleBy(res,Hp) then return 0; end if;
        res := Evaluate(Resultant(Evaluate(phi,[X,Y]),HY,Y),[T,0]);
        if not IsDivisibleBy(res,H) then return 0; end if;
        res := ExactQuotient(res,H);
        if not IsDivisibleBy(res,H) then return 0; end if;
        return 1;
    end while;
end function;

// returns D if H=H_D and 0 otherwise
function AlgorithmTwo(H:Verbose:=false)
    // Uncomment the lines below for general use
    // if not IsMonic(H) then return false; end if;
    // if not IsIrreducible(H) then return false; end if;
    h := Degree(H);
    h2list := [d: d in Divisors(h) | d eq 2^Floor(Log(2,d)) and (d-h) mod 2 eq 0];
    pmin := 33*Ceiling(h^2*(Log(Log(h+2))+2)^2); // heuristically 33 is about right (331 is definitely higher than needed)
    p := pmin-1; cnt := 0;
    while true do
        p := NextPrime(p); cnt +:= 1; // Remark 11 is not applied here (this is asymptoically suboptimal)
        Fp := GF(p);
        Hp := ChangeRing(H,Fp);
        R<x> := PolynomialRing(GF(p));
        S<z> := quo<R|Hp>;
        r := z^p-z;
        d := Degree(GCD(R!r,Hp));
        if d eq 0 then continue; end if;
        if GCD(Hp,Derivative(Hp)) ne 1 then continue; end if;
        if d lt h and not d in h2list then return 0; end if;
        j := Roots(Hp)[1][1];
        E := EllipticCurveFromjInvariant(j);
        if IsSupersingular(E) then continue; end if;
        if Verbose then printf "Found a suitable prime p = %o after testing %o primes\n", p, cnt; end if;
        _,_,D := EndomorphismRingData(E); // Remark 12 is not applied here
        if ClassNumber(D) ne Degree(Hp) then return 0; end if;
        return H eq HilbertClassPolynomial(D) select D else 0;
    end while;
end function;

D1000 := [-d:d in [3..1000]|IsDiscriminant(-d)];
assert #D1000 eq 500;

AbelianD1 := [-3, -4, -7, -8, -11, -12, -16, -19, -27, -28, -43, -67, -163];
AbelianD2 := [ -15, -20, -24, -32, -35, -36, -40, -48, -51, -52, -60, -64, -72, -75, -88, -91, -99, -100, -112, -115, -123, -147, -148, -187, -232, -235, -267, -403, -427];
AbelianD4 := [ -84, -96, -120, -132, -160, -168, -180, -192, -195, -228, -240, -280, -288, -312, -315, -340, -352, -372, -408, -435, -448, -483, -520, -532, -555, -595, -627, -708, -715, -760, -795, -928, -1012, -1435 ];
AbelianD8 := [ -420, -480, -660, -672, -840, -960, -1092, -1120, -1155, -1248, -1320, -1380, -1428, -1540, -1632, -1848, -1995, -2080, -3003, -3040, -3315 ];
AbelianD16 := [ -3360, -5280, -5460, -7392 ];
AbelianD := AbelianD1 cat AbelianD2 cat AbelianD4 cat AbelianD8 cat AbelianD16;
assert #AbelianD eq 101;

// [least,median,largest] discriminants of class number 5,10,15,...,100
hD100 := [[-47,-571,-2683],[-119,-2299,-13843],[-239,-6571,-34483],[-455,-9124,-58843],[-479,-25747,-93307],[-671,-21592,-137083],[-1031,-42499,-210907],[-1271,-34180,-274003],
         [-1319,-60748,-308323],[-1799,-64203,-389467],[-4463,-101419,-452083],[-2159,-83176,-662803],[-3527,-138883,-703123],[-3239,-138979,-821683],[-4703,-157051,-916507],
         [-5183,-133620,-1165483],[-4079,-252988,-1285747],[-5951,-204619,-1548523],[-6959,-251443,-1659067],[-7991,-249451,-1856563]];
assert #hD100 eq 20;

// [least,median,largest] fundamental discriminants of class number 125,150,175,...,1000
hD1000 := [[-11519,-570139,-2944363],[-11879,-658891,-4163443],[-18719,-1124131,-4742467],[-21311,-984952,-7206763],[-25031,-1911787,-9119203],[-30551,-1881987,-12208267],[-38351,-2682227,-14420323],
           [-34271,-2497819,-19184323],[-38639,-3414379,-19385683],[-47759,-3724003,-24669283],[-53231,-5552779,-26809843],[-67031,-3890859,-32492923],[-84719,-6756619,-35307787],
           [-75599,-6503512,-40868683],[-99839,-8428291,-43078963],[-96599,-6895428,-53698747],[-117911,-9473059,-51065683],[-127151,-9325347,-69758203],[-136751,-11505827,-61595227],
           [-148511,-9248548,-72789403],[-234599,-13993459,-81289723],[-227015,-13278259,-85527187],[-269879,-17452579,-104367187],[-185471,-13972667,-96016243],[-226799,-19566907,-106911667],
           [-218951,-17530243,-113620123],[-265079,-21602923,-134009467],[-233999,-15304228,-129667987],[-316391,-24836347,-130788067],[-266279,-23348968,-135386227],[-250799,-27051163,-149457163],
           [-299519,-22973988,-180619363],[-314159,-30598123,-164446027],[-351911,-28525507,-186089443],[-376631,-35460091,-184039507],[-412079,-25671795,-220287043]];
assert #hD1000 eq 36;

procedure CMProfile(alg,discs:exact:=true,detail:=1)
    start := Cputime();
    R<x> := PolynomialRing(Integers());
    if detail ne 3 then printf "Computing H_D and checking H_D(x),H_D(x+1),H_D(x)+1) for %o discriminants\n", #discs; end if;
    if detail eq 3 then printf "\\begin{tabular}{rrrrrrr}\n$h$ & $|H|$ & $D$ & $t_{\\rm HCP}$ & $t_{\\rm CM}$ & $t_{\\rm noCM}$\\\\\\toprule\n"; end if; 
    for D in discs do
        t0 := Cputime(); H := HilbertClassPolynomial(D); t0 := Cputime()-t0;
        t1 := Cputime(); assert exact select alg(H) eq D else alg(H) ne 0; t1 := Cputime()-t1;
        // t2 := Cputime(); assert alg(Evaluate(H,x+1)) eq 0; t2 := Cputime()-t2;
        t3 := Cputime(); assert alg(H+1) eq 0; t3 := Cputime()-t3;
        h := Degree(H); ht := Round(Log(Max([Abs(c):c in Coefficients(H)])));
        if detail eq 1 then printf ".";
        elif detail eq 2 then printf "h=%o, |H|=%o, D=%o: %.2os, %.2os, %.2os\n",h,ht,D,t0,t1,t3;
        elif detail eq 3 then printf "%o & %o & %o & %.2o & %.2o & %.2o\\\\\n",h,ht,D,t0,t1,t3;
        end if;
    end for;
    if detail eq 3 then print "\\bottomrule\n\\end{tabular}"; end if;
    printf "Verified %o discriminants in %.3os\n", #discs, Cputime()-start;
end procedure;

print "Example usage: \"CMProfile(AlgorithmTwo,[r[2]:r in hD100:detail:=3);\" or \"CMProfile(AlgorithmOne,[r[2]:r in hD100]:exact:=false,detail:=2);\"";
