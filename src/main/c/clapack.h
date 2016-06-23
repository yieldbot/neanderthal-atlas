#ifndef CLAPACK_H

int dgelsd_(int *__m, int *__n, int *__nrhs,
        double *__a, int *__lda, double *__b,
        int *__ldb, double *__s,
        double *__rcond, int *__rank,
        double *__work, int *__lwork,
        int *__iwork,
        int *__info);

int dgesdd_(char *__jobz, int *__m, int *__n,
        double *__a, int *__lda, double *__s,
        double *__u, int *__ldu, double *__vt,
        int *__ldvt, double *__work,
        int *__lwork, int *__iwork,
        int *__info);

int dgesvd_(char *__jobu, char *__jobvt, int *__m,
        int *__n, double *__a, int *__lda,
        double *__s, double *__u, int *__ldu,
        double *__vt, int *__ldvt,
        double *__work, int *__lwork,
        int *__info);

#endif  /* end #ifdef CLAPACK_H */
