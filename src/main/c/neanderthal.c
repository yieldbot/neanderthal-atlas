#include <stdlib.h>
#include <jni.h>
#include "cblas.h"
#include "uncomplicate_neanderthal_CBLAS.h"


JavaVM *javavm;

JNIEXPORT jint JNICALL JNI_OnLoad (JavaVM *jvm, void *reserved) {
  javavm=jvm;
  return JNI_VERSION_1_2;
}

int cblas_errprn(int ierr, int info, char *form, ...) {
  JNIEnv *env;
  (*javavm)->AttachCurrentThread
    (javavm, (void **)&env, NULL);
  jclass iaexception = (*env)->FindClass
    (env, "java/lang/IllegalArgumentException");

  va_list argptr;
  va_start(argptr, form);

  char *message = malloc(vsnprintf(NULL, 0, form, argptr) + 1);
  vsprintf(message, form, argptr);

  (*env)->ThrowNew(env, iaexception, message);
  va_end(argptr);
  if (ierr < info)
    return(ierr);
  else return(info);
};

void cblas_xerbla(int p, const char *rout, const char *form, ...) {
  // Override exit(-1) of the original cblas_xerbla.
  // In ATLAS, form is empty, so we have to use cblass_errprn
  // to get the error information.
  return;
};

/*
 * ======================================================
 * Level 1 BLAS functions
 * ======================================================
 */

/*
 * ------------------------------------------------------
 * DOT
 * ------------------------------------------------------
 */

JNIEXPORT jfloat JNICALL Java_uncomplicate_neanderthal_CBLAS_sdsdot
(JNIEnv *env, jclass clazz, jint N, jfloat alpha,
 jobject X, jint incX,
 jobject Y, jint incY) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  return cblas_sdsdot(N, alpha, cX, incX, cY, incY);
};


JNIEXPORT jdouble JNICALL Java_uncomplicate_neanderthal_CBLAS_dsdot
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint incX,
 jobject Y, jint incY) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  return cblas_dsdot(N, cX, incX, cY, incY);
};

JNIEXPORT jdouble JNICALL Java_uncomplicate_neanderthal_CBLAS_ddot
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint incX,
 jobject Y, jint incY) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  return cblas_ddot(N, cX, incX, cY, incY);
};

JNIEXPORT jfloat JNICALL Java_uncomplicate_neanderthal_CBLAS_sdot
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint incX,
 jobject Y, jint incY) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  return cblas_sdot(N, cX, incX, cY, incY);
};

/*
 * ------------------------------------------------------
 * NRM2
 * ------------------------------------------------------
 */

JNIEXPORT jfloat JNICALL Java_uncomplicate_neanderthal_CBLAS_snrm2
(JNIEnv *env, jclass clazz, jint N, jobject X, jint incX) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  return cblas_snrm2(N, cX, incX);
};

JNIEXPORT jdouble JNICALL Java_uncomplicate_neanderthal_CBLAS_dnrm2
(JNIEnv *env, jclass clazz, jint N, jobject X, jint incX) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  return cblas_dnrm2(N, cX, incX);
};

/*
 * ------------------------------------------------------
 * ASUM
 * ------------------------------------------------------
 */

JNIEXPORT jfloat JNICALL Java_uncomplicate_neanderthal_CBLAS_sasum
(JNIEnv *env, jclass clazz, jint N, jobject X, jint incX) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  return cblas_sasum(N, cX, incX);
};

JNIEXPORT jdouble JNICALL Java_uncomplicate_neanderthal_CBLAS_dasum
(JNIEnv *env, jclass clazz, jint N, jobject X, jint incX) {
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  return cblas_dasum(N, cX, incX);
};


/*
 * ------------------------------------------------------
 * BLAS PLUS: SUM
 * ------------------------------------------------------
 */

JNIEXPORT jfloat JNICALL Java_uncomplicate_neanderthal_CBLAS_ssum
(JNIEnv *env, jclass clazz, jint N, jobject X, jint incX) {

    float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);

    const int stride = incX;
    const int step = 16;
    const int tail = N % step;
    const int n = N - tail;

    float *end = cX + n * incX;

    float acc = 0.0f;
    float acc0 = 0.0f;
    float acc1 = 0.0f;
    float acc2 = 0.0f;
    float acc3 = 0.0f;
    float acc4 = 0.0f;
    float acc5 = 0.0f;
    float acc6 = 0.0f;
    float acc7 = 0.0f;
    float acc8 = 0.0f;
    float acc9 = 0.0f;
    float acc10 = 0.0f;
    float acc11 = 0.0f;
    float acc12 = 0.0f;
    float acc13 = 0.0f;
    float acc14 = 0.0f;
    float acc15 = 0.0f;

    while (cX != end) {
        acc0 += cX[0];
        acc1 += cX[1];
        acc2 += cX[2];
        acc3 += cX[3];
        acc4 += cX[4];
        acc5 += cX[5];
        acc6 += cX[6];
        acc7 += cX[7];
        acc8 += cX[8];
        acc9 += cX[9];
        acc10 += cX[10];
        acc11 += cX[11];
        acc12 += cX[12];
        acc13 += cX[13];
        acc14 += cX[14];
        acc15 += cX[15];
        cX += 16 * stride;
    }

    for (int i = 0; i < tail; i++) {
        acc += end[i * stride];
    }

    return acc + acc0 + acc1 + acc2 + acc3 + acc4 + acc5 + acc6 + acc7
        + acc8 + acc9 + acc10 + acc11 + acc12 + acc13 + acc14 + acc15;
};

JNIEXPORT jdouble JNICALL Java_uncomplicate_neanderthal_CBLAS_dsum
(JNIEnv *env, jclass clazz, jint N, jobject X, jint incX) {

    double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);

    const int stride = incX;
    const int step = 16;
    const int tail = N % step;
    const int n = N - tail;

    double *end = cX + n * incX;

    double acc = 0.0f;
    double acc0 = 0.0f;
    double acc1 = 0.0f;
    double acc2 = 0.0f;
    double acc3 = 0.0f;
    double acc4 = 0.0f;
    double acc5 = 0.0f;
    double acc6 = 0.0f;
    double acc7 = 0.0f;
    double acc8 = 0.0f;
    double acc9 = 0.0f;
    double acc10 = 0.0f;
    double acc11 = 0.0f;
    double acc12 = 0.0f;
    double acc13 = 0.0f;
    double acc14 = 0.0f;
    double acc15 = 0.0f;

    while (cX != end) {
        acc0 += cX[0];
        acc1 += cX[1];
        acc2 += cX[2];
        acc3 += cX[3];
        acc4 += cX[4];
        acc5 += cX[5];
        acc6 += cX[6];
        acc7 += cX[7];
        acc8 += cX[8];
        acc9 += cX[9];
        acc10 += cX[10];
        acc11 += cX[11];
        acc12 += cX[12];
        acc13 += cX[13];
        acc14 += cX[14];
        acc15 += cX[15];
        cX += 16 * stride;
    }

    for (int i = 0; i < tail; i++) {
        acc += end[i * stride];
    }

    return acc + acc0 + acc1 + acc2 + acc3 + acc4 + acc5 + acc6 + acc7
        + acc8 + acc9 + acc10 + acc11 + acc12 + acc13 + acc14 + acc15;
};


/*
 * ------------------------------------------------------
 * IAMAX
 * ------------------------------------------------------
 */

JNIEXPORT jint JNICALL Java_uncomplicate_neanderthal_CBLAS_isamax
(JNIEnv *env, jclass clazz, jint N, jobject X, jint incX) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  return cblas_isamax(N, cX, incX);
};

JNIEXPORT jint JNICALL Java_uncomplicate_neanderthal_CBLAS_idamax
(JNIEnv *env, jclass clazz, jint N, jobject X, jint incX) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  return cblas_idamax(N, cX, incX);
};

/*
 * ======================================================
 * Level 1 BLAS procedures
 * ======================================================
 */

/*
 * ------------------------------------------------------
 * ROT
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_srot
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint incX,
 jobject Y, jint incY,
 jfloat c, jfloat s) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_srot(N, cX, incX, cY, incY, c, s);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_drot
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint incX,
 jobject Y, jint incY,
 jdouble c, jdouble s) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_drot(N, cX, incX, cY, incY, c, s);
};

/*
 * ------------------------------------------------------
 * ROTG
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_srotg
(JNIEnv *env, jclass clazz, jobject params) {

  float *ca = (float *) (*env)->GetDirectBufferAddress(env, params);
  float *cb = ca + 1;
  float *cc = ca + 2;
  float *cs = ca + 3;
  cblas_srotg(ca, cb, cc, cs);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_drotg
(JNIEnv *env, jclass clazz, jobject params) {

  double *ca = (double *) (*env)->GetDirectBufferAddress(env, params);
  double *cb = ca + 1;
  double *cc = ca + 2;
  double *cs = ca + 3;
  cblas_drotg(ca, cb, cc, cs);
};

/*
 * ------------------------------------------------------
 * ROTMG
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_srotmg
(JNIEnv *env, jclass clazz, jobject args, jobject P) {

  float *cargs = (float *) (*env)->GetDirectBufferAddress(env, args);
  float *cP = (float *) (*env)->GetDirectBufferAddress(env, P);
  cblas_srotmg(cargs, cargs + 1, cargs + 2, cargs[3], cP);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_drotmg
(JNIEnv *env, jclass clazz, jobject args, jobject P) {

  double *cargs = (double *) (*env)->GetDirectBufferAddress(env, args);
  double *cP = (double *) (*env)->GetDirectBufferAddress(env, P);
  cblas_drotmg(cargs, cargs + 1, cargs + 2, cargs[3], cP);
};

/*
 * ------------------------------------------------------
 * ROTM
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_srotm
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint incX,
 jobject Y, jint incY,
 jobject P) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  float *cP = (float *) (*env)->GetDirectBufferAddress(env, P);
  cblas_srotm(N, cX, incX, cY, incY, cP);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_drotm
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint incX,
 jobject Y, jint incY,
 jobject P) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  double *cP = (double *) (*env)->GetDirectBufferAddress(env, P);
  cblas_drotm(N, cX, incX, cY, incY, cP);
};

/*
 * ------------------------------------------------------
 * SWAP
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_sswap
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint offsetX, jint incX,
 jobject Y, jint offsetY, jint incY) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_sswap(N, cX + offsetX, incX, cY + offsetY, incY);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dswap
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint offsetX, jint incX,
 jobject Y, jint offsetY, jint incY) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_dswap(N, cX + offsetX, incX, cY + offsetY, incY);

};

/*
 * ------------------------------------------------------
 * SCAL
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_sscal
(JNIEnv *env, jclass clazz,
 jint N, jfloat alpha,
 jobject X, jint incX) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  cblas_sscal(N, alpha, cX, incX);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dscal
(JNIEnv *env, jclass clazz,
 jint N, jdouble alpha,
 jobject X, jint incX) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  cblas_dscal(N, alpha, cX, incX);
};

/*
 * ------------------------------------------------------
 * COPY
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_scopy
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint offsetX, jint incX,
 jobject Y, jint offsetY, jint incY) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_scopy(N, cX + offsetX, incX, cY + offsetY, incY);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dcopy
(JNIEnv *env, jclass clazz, jint N,
 jobject X, jint offsetX, jint incX,
 jobject Y, jint offsetY, jint incY) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_dcopy(N, cX + offsetX, incX, cY + offsetY, incY);
};

/*
 * ------------------------------------------------------
 * AXPY
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_saxpy
(JNIEnv *env, jclass clazz,
 jint N, jfloat alpha,
 jobject X, jint incX,
 jobject Y, jint incY) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_saxpy(N, alpha, cX, incX, cY, incY);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_daxpy
(JNIEnv *env, jclass clazz,
 jint N, jdouble alpha,
 jobject X, jint incX,
 jobject Y, jint incY) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_daxpy(N, alpha, cX, incX, cY, incY);
};

/*
 * ======================================================
 * Level 2 BLAS procedures
 * ======================================================
 */

/*
 * ------------------------------------------------------
 * GEMV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_sgemv
(JNIEnv * env, jclass clazz,
 jint Order, jint TransA,
 jint M, jint N,
 jfloat alpha,
 jobject A, jint lda,
 jobject X, jint incX,
 jfloat beta,
 jobject Y, jint incY) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_sgemv(Order, TransA, M, N, alpha, cA, lda, cX, incX, beta, cY, incY);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dgemv
(JNIEnv * env, jclass clazz,
 jint Order, jint TransA,
 jint M, jint N,
 jdouble alpha,
 jobject A, jint lda,
 jobject X, jint incX,
 jdouble beta,
 jobject Y, jint incY) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_dgemv(Order, TransA, M, N, alpha, cA, lda, cX, incX, beta, cY, incY);
};

/*
 * ------------------------------------------------------
 * GBMV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_sgbmv
(JNIEnv *env, jclass clazz,
 jint Order, jint TransA,
 jint M, jint N,
 jint KL, jint KU,
 jfloat alpha,
 jobject A, jint lda,
 jobject X, jint incX,
 jfloat beta,
 jobject Y, jint incY) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_sgbmv(Order, TransA, M, N, KL, KU,
              alpha, cA, lda, cX, incX, beta, cY, incY);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dgbmv
(JNIEnv *env, jclass clazz,
 jint Order, jint TransA,
 jint M, jint N,
 jint KL, jint KU,
 jdouble alpha,
 jobject A, jint lda,
 jobject X, jint incX,
 jdouble beta,
 jobject Y, jint incY) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_dgbmv(Order, TransA, M, N, KL, KU,
              alpha, cA, lda, cX, incX, beta, cY, incY);
};

/*
 * ------------------------------------------------------
 * SYMV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_ssymv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jfloat alpha,
 jobject A, jint lda,
 jobject X, jint incX,
 jfloat beta,
 jobject Y, jint incY) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_ssymv(Order, Uplo, N, alpha, cA, lda, cX, incX, beta, cY, incY);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dsymv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jdouble alpha,
 jobject A, jint lda,
 jobject X, jint incX,
 jdouble beta,
 jobject Y, jint incY) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_dsymv(Order, Uplo, N, alpha, cA, lda, cX, incX, beta, cY, incY);
};

/*
 * ------------------------------------------------------
 * SBMV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_ssbmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N, jint K,
 jfloat alpha,
 jobject A, jint lda,
 jobject X, jint incX,
 jfloat beta,
 jobject Y, jint incY) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_ssbmv(Order, Uplo, N, K, alpha, cA, lda, cX, incX, beta, cY, incY);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dsbmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N, jint K,
 jdouble alpha,
 jobject A, jint lda,
 jobject X, jint incX,
 jdouble beta,
 jobject Y, jint incY) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_dsbmv(Order, Uplo, N, K, alpha, cA, lda, cX, incX, beta, cY, incY);
};

/*
 * ------------------------------------------------------
 * SPMV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_sspmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jfloat alpha,
 jobject Ap,
 jobject X, jint incX,
 jfloat beta,
 jobject Y, jint incY) {

  float *cAp = (float *) (*env)->GetDirectBufferAddress(env, Ap);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_sspmv(Order, Uplo, N, alpha, cAp, cX, incX, beta, cY, incY);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dspmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jdouble alpha,
 jobject Ap,
 jobject X, jint incX,
 jdouble beta,
 jobject Y, jint incY) {

  double *cAp = (double *) (*env)->GetDirectBufferAddress(env, Ap);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  cblas_dspmv(Order, Uplo, N, alpha, cAp, cX, incX, beta, cY, incY);
};

/*
 * ------------------------------------------------------
 * TRMV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_strmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA,
 jint N, jfloat alpha,
 jobject A, jint lda,
 jobject X, jint incX) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  cblas_strmv(Order, Uplo, TransA, N, alpha, cA, lda, cX, incX);
};


JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dtrmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA,
 jint N, jdouble alpha,
 jobject A, jint lda,
 jobject X, jint incX) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  cblas_dtrmv(Order, Uplo, TransA, N, alpha, cA, lda, cX, incX);
};

/*
 * ------------------------------------------------------
 * TBMV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_stbmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N, jint K,
 jobject A, jint lda,
 jobject X, jint incX) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  cblas_strmv(Order, Uplo, TransA, Diag, N, cA, lda, cX, incX);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dtbmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N, jint K,
 jobject A, jint lda,
 jobject X, jint incX) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  cblas_dtrmv(Order, Uplo, TransA, Diag, N, cA, lda, cX, incX);
};

/*
 * ------------------------------------------------------
 * TPMV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_stpmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N,
 jobject Ap,
 jobject X, jint incX) {

  float *cAp = (float *) (*env)->GetDirectBufferAddress(env, Ap);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  cblas_stpmv(Order, Uplo, TransA, Diag, N, cAp, cX, incX);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dtpmv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N,
 jobject Ap,
 jobject X, jint incX) {

  double *cAp = (double *) (*env)->GetDirectBufferAddress(env, Ap);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  cblas_dtpmv(Order, Uplo, TransA, Diag, N, cAp, cX, incX);
};

/*
 * ------------------------------------------------------
 * TRSV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_strsv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N,
 jobject A, jint lda,
 jobject X, jint incX) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  cblas_strsv(Order, Uplo, TransA, Diag, N, cA, lda, cX, incX);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dtrsv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N,
 jobject A, jint lda,
 jobject X, jint incX) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  cblas_dtrsv(Order, Uplo, TransA, Diag, N, cA, lda, cX, incX);
};

/*
 * ------------------------------------------------------
 * TBSV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_stbsv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N, jint K,
 jobject A, jint lda,
 jobject X, jint incX) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  cblas_stbsv(Order, Uplo, TransA, Diag, N, K, cA, lda, cX, incX);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dtbsv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N, jint K,
 jobject A, jint lda,
 jobject X, jint incX) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  cblas_dtbsv(Order, Uplo, TransA, Diag, N, K, cA, lda, cX, incX);
};

/*
 * ------------------------------------------------------
 * TPSV
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_stpsv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N,
 jobject Ap,
 jobject X, jint incX) {

  float *cAp = (float *) (*env)->GetDirectBufferAddress(env, Ap);
  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  cblas_stpsv(Order, Uplo, TransA, Diag, N, cAp, cX, incX);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dtpsv
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint TransA, jint Diag,
 jint N,
 jobject Ap,
 jobject X, jint incX) {

  double *cAp = (double *) (*env)->GetDirectBufferAddress(env, Ap);
  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  cblas_dtpsv(Order, Uplo, TransA, Diag, N, cAp, cX, incX);
};

/*
 * ------------------------------------------------------
 * GER
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_sger
(JNIEnv *env, jclass clazz,
 jint Order,
 jint M, jint N,
 jfloat alpha,
 jobject X, jint incX,
 jobject Y, jint incY,
 jobject A, jint lda) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  cblas_sger(Order, M, N, alpha, cX, incX, cY, incY, cA, lda);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dger
(JNIEnv *env, jclass clazz,
 jint Order,
 jint M, jint N,
 jdouble alpha,
 jobject X, jint incX,
 jobject Y, jint incY,
 jobject A, jint lda) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  cblas_dger(Order, M, N, alpha, cX, incX, cY, incY, cA, lda);
};

/*
 * ------------------------------------------------------
 * SYR
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_ssyr
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jfloat alpha,
 jobject X, jint incX,
 jobject A, jint lda) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  cblas_ssyr(Order, Uplo, N, alpha, cX, incX, cA, lda);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dsyr
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jdouble alpha,
 jobject X, jint incX,
 jobject A, jint lda) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  cblas_dsyr(Order, Uplo, N, alpha, cX, incX, cA, lda);
}

/*
 * ------------------------------------------------------
 * SPR
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_sspr
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jfloat alpha,
 jobject X, jint incX,
 jobject Ap) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cAp = (float *) (*env)->GetDirectBufferAddress(env, Ap);
  cblas_sspr(Order, Uplo, N, alpha, cX, incX, cAp);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dspr
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jdouble alpha,
 jobject X, jint incX,
 jobject Ap) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cAp = (double *) (*env)->GetDirectBufferAddress(env, Ap);
  cblas_dspr(Order, Uplo, N, alpha, cX, incX, cAp);
};

/*
 * ------------------------------------------------------
 * SYR2
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_ssyr2
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jfloat alpha,
 jobject X, jint incX,
 jobject Y, jint incY,
 jobject A, jint lda) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  cblas_ssyr2(Order, Uplo, N, alpha, cX, incX, cY, incY, cA, lda);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dsyr2
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jdouble alpha,
 jobject X, jint incX,
 jobject Y, jint incY,
 jobject A, jint lda) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  cblas_dsyr2(Order, Uplo, N, alpha, cX, incX, cY, incY, cA, lda);
};

/*
 * ------------------------------------------------------
 * SPR2
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_sspr2
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jfloat alpha,
 jobject X, jint incX,
 jobject Y, jint incY,
 jobject Ap) {

  float *cX = (float *) (*env)->GetDirectBufferAddress(env, X);
  float *cY = (float *) (*env)->GetDirectBufferAddress(env, Y);
  float *cAp = (float *) (*env)->GetDirectBufferAddress(env, Ap);
  cblas_sspr2(Order, Uplo, N, alpha, cX, incX, cY, incY, cAp);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dspr2
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo,
 jint N,
 jdouble alpha,
 jobject X, jint incX,
 jobject Y, jint incY,
 jobject Ap) {

  double *cX = (double *) (*env)->GetDirectBufferAddress(env, X);
  double *cY = (double *) (*env)->GetDirectBufferAddress(env, Y);
  double *cAp = (double *) (*env)->GetDirectBufferAddress(env, Ap);
  cblas_dspr2(Order, Uplo, N, alpha, cX, incX, cY, incY, cAp);
};

/*
 * ======================================================
 * Level 3 BLAS procedures
 * ======================================================
 */


/*
 * ------------------------------------------------------
 * GEMM
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_sgemm
(JNIEnv *env, jclass clazz,
 jint Order, jint TransA, jint TransB,
 jint M, jint N, jint K,
 jfloat alpha,
 jobject A, jint lda,
 jobject B, jint ldb,
 jfloat beta,
 jobject C, jint ldc) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cB = (float *) (*env)->GetDirectBufferAddress(env, B);
  float *cC = (float *) (*env)->GetDirectBufferAddress(env, C);
  cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, cA, lda, cB, ldb, beta, cC, ldc);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dgemm
(JNIEnv *env, jclass clazz,
 jint Order, jint TransA, jint TransB,
 jint M, jint N, jint K,
 jdouble alpha,
 jobject A, jint lda,
 jobject B, jint ldb,
 jdouble beta,
 jobject C, jint ldc) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cB = (double *) (*env)->GetDirectBufferAddress(env, B);
  double *cC = (double *) (*env)->GetDirectBufferAddress(env, C);
  cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, cA, lda, cB, ldb, beta, cC, ldc);
};

/*
 * ------------------------------------------------------
 * SYMM
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_ssymm
(JNIEnv *env, jclass clazz,
 jint Order, jint Side, jint Uplo,
 jint M, jint N,
 jfloat alpha,
 jobject A, jint lda,
 jobject B, jint ldb,
 jfloat beta,
 jobject C, jint ldc) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cB = (float *) (*env)->GetDirectBufferAddress(env, B);
  float *cC = (float *) (*env)->GetDirectBufferAddress(env, C);
  cblas_ssymm(Order, Side, Uplo, M, N, alpha, cA, lda, cB, ldb, beta, cC, ldc);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dsymm
(JNIEnv *env, jclass clazz,
 jint Order, jint Side, jint Uplo,
 jint M, jint N,
 jdouble alpha,
 jobject A, jint lda,
 jobject B, jint ldb,
 jdouble beta,
 jobject C, jint ldc) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cB = (double *) (*env)->GetDirectBufferAddress(env, B);
  double *cC = (double *) (*env)->GetDirectBufferAddress(env, C);
  cblas_dsymm(Order, Side, Uplo, M, N, alpha, cA, lda, cB, ldb, beta, cC, ldc);
};

/*
 * ------------------------------------------------------
 * SYRK
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_ssyrk
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint Trans,
 jint N, jint K,
 jfloat alpha,
 jobject A, jint lda,
 jfloat beta,
 jobject C, jint ldc) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cC = (float *) (*env)->GetDirectBufferAddress(env, C);
  cblas_ssyrk(Order, Uplo, Trans, N, K, alpha, cA, lda, beta, cC, ldc);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dsyrk
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint Trans,
 jint N, jint K,
 jdouble alpha,
 jobject A, jint lda,
 jfloat beta,
 jobject C, jint ldc) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cC = (double *) (*env)->GetDirectBufferAddress(env, C);
  cblas_dsyrk(Order, Uplo, Trans, N, K, alpha, cA, lda, beta, cC, ldc);
};

/*
 * ------------------------------------------------------
 * SYR2K
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_ssyr2k
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint Trans,
 jint N, jint K,
 jfloat alpha,
 jobject A, jint lda,
 jobject B, jint ldb,
 jfloat beta,
 jobject C, jint ldc) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cB = (float *) (*env)->GetDirectBufferAddress(env, B);
  float *cC = (float *) (*env)->GetDirectBufferAddress(env, C);
  cblas_ssyr2k(Order, Uplo, Trans, N, K, alpha, cA, lda, cB, ldb, beta, cC, ldc);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dsyr2k
(JNIEnv *env, jclass clazz,
 jint Order, jint Uplo, jint Trans,
 jint N, jint K,
 jdouble alpha,
 jobject A, jint lda,
 jobject B, jint ldb,
 jdouble beta,
 jobject C, jint ldc) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cB = (double *) (*env)->GetDirectBufferAddress(env, B);
  double *cC = (double *) (*env)->GetDirectBufferAddress(env, C);
  cblas_dsyr2k(Order, Uplo, Trans, N, K, alpha, cA, lda, cB, ldb, beta, cC, ldc);
};

/*
 * ------------------------------------------------------
 * TRMM
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_strmm
(JNIEnv *env, jclass clazz,
 jint Order, jint Side,
 jint Uplo, jint TransA, jint Diag,
 jint M, jint N,
 jfloat alpha,
 jobject A, jint lda,
 jobject B, jint ldb) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cB = (float *) (*env)->GetDirectBufferAddress(env, B);
  cblas_strmm(Order, Side, Uplo, TransA, Diag, M, N, alpha, cA, lda, cB, ldb);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dtrmm
(JNIEnv *env, jclass clazz,
 jint Order, jint Side,
 jint Uplo, jint TransA, jint Diag,
 jint M, jint N,
 jdouble alpha,
 jobject A, jint lda,
 jobject B, jint ldb) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cB = (double *) (*env)->GetDirectBufferAddress(env, B);
  cblas_dtrmm(Order, Side, Uplo, TransA, Diag, M, N, alpha, cA, lda, cB, ldb);
};

/*
 * ------------------------------------------------------
 * TRSM
 * ------------------------------------------------------
 */

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_strsm
(JNIEnv *env, jclass clazz,
 jint Order, jint Side,
 jint Uplo, jint TransA, jint Diag,
 jint M, jint N,
 jfloat alpha,
 jobject A, jint lda,
 jobject B, jint ldb) {

  float *cA = (float *) (*env)->GetDirectBufferAddress(env, A);
  float *cB = (float *) (*env)->GetDirectBufferAddress(env, B);
  cblas_strsm(Order, Side, Uplo, TransA, Diag, M, N, alpha, cA, lda, cB, ldb);
};

JNIEXPORT void JNICALL Java_uncomplicate_neanderthal_CBLAS_dtrsm
(JNIEnv *env, jclass clazz,
 jint Order, jint Side,
 jint Uplo, jint TransA, jint Diag,
 jint M, jint N,
 jdouble alpha,
 jobject A, jint lda,
 jobject B, jint ldb) {

  double *cA = (double *) (*env)->GetDirectBufferAddress(env, A);
  double *cB = (double *) (*env)->GetDirectBufferAddress(env, B);
  cblas_dtrsm(Order, Side, Uplo, TransA, Diag, M, N, alpha, cA, lda, cB, ldb);
};

/*
 * ------------------------------------------------------
 * LAPACK
 * ------------------------------------------------------
 */

JNIEXPORT jint JNICALL Java_uncomplicate_neanderthal_LAPACK_dgelsd_1
  (JNIEnv *env, jclass clazz,
  jint m, jint n, jint nrhs,
  jobject a, jint lda, jobject b, jint ldb, jobject s,
  jdouble rcond, jobject rank, jobject work, jint lwork, jobject iwork)
{
  double *cA = (double *) (*env)->GetDirectBufferAddress(env, a);
  double *cB = (double *) (*env)->GetDirectBufferAddress(env, b);
  double *cS = (double *) (*env)->GetDirectBufferAddress(env, s);
  int *cRank = (int *) (*env)->GetDirectBufferAddress(env, rank);
  double *cWork = (double *) (*env)->GetDirectBufferAddress(env, work);
  int *cIwork = (int *) (*env)->GetDirectBufferAddress(env, iwork);
  int info;
  dgelsd_(&m, &n, &nrhs, cA, &lda, cB, &ldb, cS, &rcond, cRank, cWork, &lwork, cIwork, &info);
  return info;
}

JNIEXPORT jint JNICALL Java_uncomplicate_neanderthal_LAPACK_dgesdd_1
(JNIEnv *env, jclass clazz, jchar jobz, jint m, jint n, jobject a,
 jint lda, jobject s, jobject u, jint ldu, jobject vt, jint ldvt,
 jobject work, jint lwork, jobject iwork)
{
  double *cA = (double *) (*env)->GetDirectBufferAddress(env, a);
  double *cS = (double *) (*env)->GetDirectBufferAddress(env, s);
  double *cU = (double *) (*env)->GetDirectBufferAddress(env, u);
  double *cVt = (double *) (*env)->GetDirectBufferAddress(env, vt);
  double *cWork = (double *) (*env)->GetDirectBufferAddress(env, work);
  int *cIwork = (int *) (*env)->GetDirectBufferAddress(env, iwork);
  int info;
  dgesdd_(&jobz, &m, &n, cA, &lda, cS, cU, &ldu, cVt, &ldvt, cWork, &lwork, cIwork, &info);
  return info;
}

JNIEXPORT jint JNICALL Java_uncomplicate_neanderthal_LAPACK_dgesvd_1
(JNIEnv *env, jclass clazz,
 jchar jobu, jchar jobvt, jint m, jint n,
 jobject a, jint lda, jobject s, jobject u, jint ldu,
 jobject vt, jint ldvt, jobject work, jint lwork) {
  double *cA = (double *) (*env)->GetDirectBufferAddress(env, a);
  double *cS = (double *) (*env)->GetDirectBufferAddress(env, s);
  double *cU = (double *) (*env)->GetDirectBufferAddress(env, u);
  double *cVt = (double *) (*env)->GetDirectBufferAddress(env, vt);
  double *cWork = (double *) (*env)->GetDirectBufferAddress(env, work);
  int info;
  dgesvd_(&jobu, &jobvt, &m, &n, cA, &lda, cS, cU, &ldu, cVt, &ldvt, cWork, &lwork, &info);
  return info;
}
