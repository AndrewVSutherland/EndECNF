\\ gp version of cmdisc

default(parisizemax,4000000000)
default(nbthreads,1)

{binomial_p(r,s,p)=
 if (r < s, return(Mod(0,p)));
 if (s < r-s, s=r-s);
 return(prod(i=s+1,r,i,Mod(1,p)) / prod(i=1,r-s,i,Mod(1,p)));
}

{SupersingularPolynomial(p)=
\\ Implementation of Finotti's forumla for the supersingular polynomial
  if (p <= 5, return(Mod(1,p)*x));
  r = (p-1) / 2;
  r1 = ceil(r/3);
  r2 = floor(r/2);
  s1 = floor(r/3);
  s2 = ceil(r/2);
  c = Mod(-27,p)/4;
  c1728 = Mod(1728,p);
  f = Mod(-2,p)^r*sum(i=r1,r2,binomial_p(r,i,p)*binomial_p(i,3*i-r,p)*c^i*x^(i-s1)*(x-c1728)^(s2-i),Mod(0,p));
  return(f);
}


{AlgorithmOne(H, check_monic_irreducible=0, verbose=0)=
\\ Returns 1 if H is a Hilbert class polynomial, 0 otherwise.
 my(q=101,h=poldegree(H),lh,Hmodq,l,Hmodl,phi,phimod1,res,quo,rem);

 \\ optional check that input is monic and irreducible
 if(check_monic_irreducible,
        if(pollead(H)!=1,return(0));
        if(polisirreducible(H)==0,return(0)));

 lh = floor((2.15*log(h) + 4*log(log(log(h+2))+2) + 26.19)^2);
 Hmodq = H*Mod(1,q);
 forprime(l=3,, \\ unlimited loop over primes
    if(l>lh, print("Exceeded GRH bound l(h)=",lh,", terminating"); return(0));
    Hmodl = H*Mod(1,l);
    if (poldegree(gcd(Hmodl,deriv(Hmodl))) !=0, next());
    if (SupersingularPolynomial(l)%Hmodl ==0, next());
    if(verbose, print("Using l =",l); return(0));
    phi = polmodular(l);
    phimodq = phi*Mod(1,q);
    res = polresultant(phimodq, Hmodq, y);
    [quo,rem] = divrem(res,Hmodq);
    if (rem != 0, return(0));
    if (quo%Hmodq != 0, return(0));
    res = polresultant(phi,subst(H,x,y),y);
    [quo,rem] = divrem(res,H);
    if (rem != 0, return(0));
    if (quo%H != 0, return(0));
    return(1)
    );
}

{AlgorithmTwo(H, check_monic_irreducible=0, verbose=0) =
\\ returns D if H=H_D and 0 otherwise
 my(h=poldegree(H), d, h2list, n, p, Hp, z, t, nroots, jp, e, D, D0, fp, D1, pol);

 \\ optional check that input is monic and irreducible
 if(check_monic_irreducible,
        if(pollead(H)!=1,return(0));
        if(polisirreducible(H)==0,return(0)));

 h2list = [];
 fordiv(h,d,if(d==2^valuation(d,2) && (d-h)%2==0, h2list = concat(h2list,[d])));
 n = 0;
 pmin = 33 * ceil(h^2 * (log(log(h+2))+2)^2);
 forprime(p=pmin,,
  n +=1;
  Hp = H*Mod(1,p);
  z = Mod(variable(Hp),Hp);
  t = lift(z^p-z);
  nroots = if(t==0,h,poldegree(gcd(t,Hp)));
  if(nroots==0, next());
  if(poldisc(Hp)==0,next());
  if((nroots<h) && vecsearch(h2list,nroots)==0, return(0));
  jp = polrootsmod(Hp)[1];
  e = ellinit(ellfromj(jp));
  if(ellissupersingular(e), next());
  if(verbose, print("Found a suitable prime p = ",p," after testing ",n," primes"));
  t = ellap(e);
  D = t^2-4*p;
  [D0,f] = core(D,1);
  if(D0%4>1, D0 *= 4; f/=2);
  fordiv(f, d, D1 = D0*d^2;
               if(h != qfbclassno(D1), next());
               if(H == polclass(D1), return(D1));
              );
  return(0)
  );
};

D1000 = select(x->x%4<2, vector(1000,i,-i));
if(#D1000 != 500, error());

AbelianD1 = [-3, -4, -7, -8, -11, -12, -16, -19, -27, -28, -43, -67, -163];
AbelianD2 = [ -15, -20, -24, -32, -35, -36, -40, -48, -51, -52, -60, -64, -72, -75, -88, -91, -99, -100, -112, -115, -123, -147, -148, -187, -232, -235, -267, -403, -427];
AbelianD4 = [ -84, -96, -120, -132, -160, -168, -180, -192, -195, -228, -240, -280, -288, -312, -315, -340, -352, -372, -408, -435, -448, -483, -520, -532, -555, -595, -627, -708, -715, -760, -795, -928, -1012, -1435 ];
AbelianD8 = [ -420, -480, -660, -672, -840, -960, -1092, -1120, -1155, -1248, -1320, -1380, -1428, -1540, -1632, -1848, -1995, -2080, -3003, -3040, -3315 ];
AbelianD16 = [ -3360, -5280, -5460, -7392 ];
AbelianD = concat(AbelianD1, concat(AbelianD2, concat(AbelianD4, concat(AbelianD8,AbelianD16))));
if(#AbelianD != 101, error());

\\ [least,median,largest] discriminants of class number 5,10,15,...,100
{
hD100 = [[-47,-571,-2683],[-119,-2299,-13843],[-239,-6571,-34483],[-455,-9124,-58843],[-479,-25747,-93307],[-671,-21592,-137083],[-1031,-42499,-210907],[-1271,-34180,-274003],
         [-1319,-60748,-308323],[-1799,-64203,-389467],[-4463,-101419,-452083],[-2159,-83176,-662803],[-3527,-138883,-703123],[-3239,-138979,-821683],[-4703,-157051,-916507],
         [-5183,-133620,-1165483],[-4079,-252988,-1285747],[-5951,-204619,-1548523],[-6959,-251443,-1659067],[-7991,-249451,-1856563]];
}
if(#hD100 != 20, error());

\\ [least,median,largest] fundamental discriminants of class number 125,150,175,...,1000
{
hD1000 = [[-11519,-570139,-2944363],[-11879,-658891,-4163443],[-18719,-1124131,-4742467],[-21311,-984952,-7206763],[-25031,-1911787,-9119203],[-30551,-1881987,-12208267],[-38351,-2682227,-14420323],
           [-34271,-2497819,-19184323],[-38639,-3414379,-19385683],[-47759,-3724003,-24669283],[-53231,-5552779,-26809843],[-67031,-3890859,-32492923],[-84719,-6756619,-35307787],
           [-75599,-6503512,-40868683],[-99839,-8428291,-43078963],[-96599,-6895428,-53698747],[-117911,-9473059,-51065683],[-127151,-9325347,-69758203],[-136751,-11505827,-61595227],
           [-148511,-9248548,-72789403],[-234599,-13993459,-81289723],[-227015,-13278259,-85527187],[-269879,-17452579,-104367187],[-185471,-13972667,-96016243],[-226799,-19566907,-106911667],
           [-218951,-17530243,-113620123],[-265079,-21602923,-134009467],[-233999,-15304228,-129667987],[-316391,-24836347,-130788067],[-266279,-23348968,-135386227],[-250799,-27051163,-149457163],
           [-299519,-22973988,-180619363],[-314159,-30598123,-164446027],[-351911,-28525507,-186089443],[-376631,-35460091,-184039507],[-412079,-25671795,-220287043]];
}
if(#hD1000 != 36, error());


{CMProfile(alg, discs, exact=1, detail=1) =
    start = getabstime();

    if( detail == 3,
        print("\\begin{tabular}{rrrrrrr}\n$h$ & $|H|$ & $D$ & $t_{\\rm HCP}$ & $t_{\\rm CM}$ & $t_{\\rm noCM}$\\\\\\toprule"),
        printf("Computing H_D and checking H_D(x),H_D(x)+1) for %d discriminants\n", #discs)
      );

    foreach(discs, D,
        t0 = getabstime();
        H = polclass(D);
        t0 = getabstime()-t0;
        t1 = getabstime();
        if(exact, if(alg(H)!=D, error()), if(alg(H)==0, error()));
        t1 = getabstime()-t1;
        \\t2 = getabstime();
        \\if(alg(subst(H,x,x+1))!=0, error());
        \\t2 = getabstime()-t2;
        t3 = getabstime();
        if(alg(H+1)!=0, error());
        t3 = getabstime()-t3;
        h = poldegree(H);
        ht = round(log(vecmax(vector(h+1,i,abs(polcoeff(H,i-1))))));
        if(detail==1, print1("."));
        if(detail==2, printf("h=%d, |H|=%d, D=%d: %.2fs, %.2fs, %.2fs\n",h,ht,D,t0/1000,t1/1000,t3/1000));
        if(detail==3, printf("%d & %d & %d & %.2f & %.2f & %.2f\\\\\n",h,ht,D,t0/1000,t1/1000,t3/1000))
    );
    if(detail == 1, print(), if (detail == 3, print("\\bottomrule\n\\end{tabular}")));

    printf("Verified %d discriminants in %.3fs\n", #discs, (getabstime()-start)/1000)
};

\\display usage instruction on loading this file

print("Example usage:")
print("CMProfile(AlgorithmTwo,[x[2] | x <- hD100], exact=1, detail=3)")
print("or")
print("CMProfile(AlgorithmOne,[x[1] | x <- hD100], exact=0, detail=2)")

