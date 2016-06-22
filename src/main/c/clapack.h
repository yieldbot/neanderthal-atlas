#ifndef CLAPACK_H

int dgesvd_(char *__jobu, char *__jobvt, int *__m,
        int *__n, double *__a, int *__lda,
        double *__s, double *__u, int *__ldu,
        double *__vt, int *__ldvt,
        double *__work, int *__lwork,
        int *__info);

#endif  /* end #ifdef CLAPACK_H */
