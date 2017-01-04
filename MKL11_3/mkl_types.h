/*******************************************************************************
* Copyright 1999-2016 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
! Content:
!      Intel(R) Math Kernel Library (MKL) types definition
!****************************************************************************/

#ifndef _MKL_TYPES_H_
#define _MKL_TYPES_H_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* MKL Complex type for single precision */
#ifndef MKL_Complex8
typedef
struct _MKL_Complex8 {
    float real;
    float imag;
} MKL_Complex8;
#endif

/* MKL Complex type for double precision */
#ifndef MKL_Complex16
typedef
struct _MKL_Complex16 {
    double real;
    double imag;
} MKL_Complex16;
#endif

/* MKL Version type */
typedef
struct {
    int    MajorVersion;
    int    MinorVersion;
    int    UpdateVersion;
    char * ProductStatus;
    char * Build;
    char * Processor;
    char * Platform;
} MKLVersion;

/* MKL integer types for LP64 and ILP64 */
#if (!defined(__INTEL_COMPILER)) & defined(_MSC_VER)
    #define MKL_INT64 __int64
    #define MKL_UINT64 unsigned __int64
#else
    #define MKL_INT64 long long int
    #define MKL_UINT64 unsigned long long int
#endif

#ifdef MKL_ILP64

/* MKL ILP64 integer types */
#ifndef MKL_INT
    #define MKL_INT MKL_INT64
#endif
#ifndef MKL_UINT
    #define MKL_UINT MKL_UINT64
#endif
#define MKL_LONG MKL_INT64

#else

/* MKL LP64 integer types */
#ifndef MKL_INT
    #define MKL_INT int
#endif
#ifndef MKL_UINT
    #define MKL_UINT unsigned int
#endif
#define MKL_LONG long int

#endif

/* MKL threading stuff. MKL Domain names */
#define MKL_DOMAIN_ALL      0
#define MKL_DOMAIN_BLAS     1
#define MKL_DOMAIN_FFT      2
#define MKL_DOMAIN_VML      3
#define MKL_DOMAIN_PARDISO  4

/* MKL CBWR stuff */

/* options */
#define MKL_CBWR_BRANCH 1
#define MKL_CBWR_ALL   ~0

/* common settings */
#define MKL_CBWR_UNSET_ALL 0
#define MKL_CBWR_OFF       0

/* branch specific values */
#define MKL_CBWR_BRANCH_OFF     1
#define MKL_CBWR_AUTO           2
#define MKL_CBWR_COMPATIBLE     3
#define MKL_CBWR_SSE2           4
#define MKL_CBWR_SSSE3          6
#define MKL_CBWR_SSE4_1         7
#define MKL_CBWR_SSE4_2         8
#define MKL_CBWR_AVX            9
#define MKL_CBWR_AVX2          10
#define MKL_CBWR_AVX512_MIC    11
#define MKL_CBWR_AVX512        12

/* error codes */
#define MKL_CBWR_SUCCESS                   0
#define MKL_CBWR_ERR_INVALID_SETTINGS     -1
#define MKL_CBWR_ERR_INVALID_INPUT        -2
#define MKL_CBWR_ERR_UNSUPPORTED_BRANCH   -3
#define MKL_CBWR_ERR_UNKNOWN_BRANCH       -4
#define MKL_CBWR_ERR_MODE_CHANGE_FAILURE  -8

/* Obsolete */
#define MKL_CBWR_SSE3           5

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _MKL_TYPES_H_ */
