// Histogram bin lookups (Func*) for tE, lens mass, piE, u0, proper motion, baseline magnitude, blend.
#include "run/histograms.h"

int FunctE(lens & l) {
   int gg = -1;

   if (l.tE <= tE_min)       gg = 0;

   else if (l.tE >= tE_max)  gg = GG;

   else {
       for (int i = 1; i <= GG; ++i) {
          if(double((l.tE - l.tEs[i-1]) * (l.tE - l.tEs[i])) < 0.0 or l.tE == l.tEs[i-1]) { gg = i - 1;  break; }
       }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(l.tE > 0.0);
   CHECK(l.tEs[1] > 0.0);
   CHECK(l.tEs[10] > 0.0);

   return(gg);
}

int FuncMl(lens & l) {
   int gg = -1;

   if (l.Ml <= mlMin())      gg = 0;
   else if (l.Ml >= mlMax()) gg = GG;
   else {
      for(int i = 1; i <= GG; ++i) {
          if (double((l.Ml - l.Mls[i-1]) * (l.Ml - l.Mls[i])) < 0.0 or l.Ml == l.Mls[i-1]) { gg = i - 1;  break; }
      }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   // The lower bound follows the active population (mlMin()), not a fixed value.
   CHECK(l.Ml >= mlMin() * 0.999);

   return(gg);
}

int FuncPi(lens & l) {
   int gg = -1;
   double lpi = std::log10(l.pirel);

   if (lpi <= pi_min)      gg = 0;
   else if (lpi >= pi_max) gg = GG;
   else {
      for (int i = 1; i <= GG; ++i) {
          if (double((lpi - l.pis[i-1]) * (lpi - l.pis[i])) < 0.0 or lpi == l.pis[i-1]) { gg = i - 1;  break; }
      }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);

   return(gg);
}

int Funcu0(lens & l){
   int gg = -1;

   if (l.u0 <= u0_min)      gg = 0;
   else if (l.u0 >= u0_max) gg = GG;
   else {
      for (int i = 1; i <= GG; ++i) {
          if (double((l.u0 - l.u0s[i-1]) * (l.u0 - l.u0s[i])) < 0.0 or l.u0 == l.u0s[i-1]) { gg = i - 1;  break; }
      }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(l.u0 >= 0.0);
   CHECK(l.u0 <= u0m);

   return(gg);
}

int FuncMu(lens & l){
   int gg = -1;
   double mur = double(l.murel * year); //mas/years

   if (mur <= mu_min)      gg = 0;
   else if (mur >= mu_max) gg = GG;
   else {
       for (int i = 1; i <= GG; ++i) {
          if (double((mur - l.mus[i-1]) * (mur - l.mus[i])) < 0.0 or mur == l.mus[i-1]) { gg = i - 1;  break; }
       }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(mur >= 0.0);

   return(gg);
}

int FuncMb(lens & l, double gadr) {
   int gg = -1;

   if (gadr <= mb_min)       gg = 0;
   else if (gadr >= mb_max)  gg = GG;
   else {
       for (int i = 1; i <= GG; ++i) {
          if (double((gadr - l.mbs[i-1]) * (gadr - l.mbs[i])) < 0.0 or gadr == l.mbs[i-1]) { gg = i - 1; break; }
       }
   }
   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(gadr >= 0.0);

   return(gg);
}

int FuncFb(lens & l, double blen) {
   int gg = -1;

   if (blen <= fb_min)       gg = 0;
   else if (blen >= fb_max)  gg = GG;
   else {
       for (int i = 1; i <= GG; ++i) {
          if(double((blen - l.fbs[i-1]) * (blen - l.fbs[i])) < 0.0 or blen == l.fbs[i-1]) { gg = i - 1; break;}
       }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(blen >= 0.0);
   CHECK(blen <= 1.0);

   return(gg);
}
