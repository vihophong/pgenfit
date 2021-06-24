## Note:

This benchmark test compare the mlhfitv1 code which was used for my PhD thesis and the current version of pgenfit using identical input parameters and unbinned data. The purpose is to see if we can reproduce the PhD fitting results using pgenfit.

In the pgenfit_pvar0 program, the  unbinfit.cc code was modified to add:

```
pvar[0]->setMin(6.93147e-12);
pvar[0]->setMax(23.3407);
```

so that the input parameters of mlhfitv1 and pgenfit are identical.

you may check the tiny difference if pgenfit without modification of pvar[0] is used