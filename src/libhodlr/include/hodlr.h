/**
 * @file   hodlr.h
 * @date   January, 2019
 * @brief  C wrapper for HODLR library
*/

#ifndef HODLR_H
#define HODLR_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/** @brief Calculate \f$\log\mathrm{det}(1-M)\f$ using HODLR approach
 *
 * Compute log det(1-A) for a matrix A of dimension given by dim. The diagonal
 * elements of A are given by diagonal which is an array of dim elements.
 * Arbitrary matrix elements A_ij are given by the callback(i,j,args). The
 * requested numerical precision is given by tolerance.  returns the matrix
 * elements of A.
 *
 * nLeaf is the size (number of rows of the matrix) of the smallest block at
 * the leaf level. The number of levels in the tree is given by
 * n_levels=log_2(dim/nLeaf).
 *
 * Values for sym_psd:
 *  - 0: generic matrix
 *  - 1: matrix is symmetric
 *  - 2: matrix is symmetric and positive definite
 *
 * @param dim           dimension of matrix M
 * @param callback      function that returns matrix elements of M
 * @param args          pointer that is passed as third argument to callback
 * @param diagonal      array with the diagonal elements of M
 * @param nLeaf         nLeaf is the dimension of the smallest block at the leaf level
 * @param tolerance     requested accuracy of result
 * @param sym_spd       specifiy whether matrix is generic (0), symmetric (1) or spd (2)
 * @retval logdet       \f$\log\mathrm{det}(1-M)\f$
 */
EXTERNC double hodlr_logdet_diagonal(int dim, double (*callback)(int,int,void *), void *args, double *diagonal, unsigned int nLeaf, double tolerance, int is_symmetric);

/** @brief Calculate log(det(Id-M)) using HODLR approach
 *
 * See \ref hodlr_logdet_diagonal for more information.
 *
 * @param dim           dimension of matrix M
 * @param callback      function that returns matrix elements of M
 * @param args          pointer that is passed as third argument to callback
 * @param nLeaf         nLeaf is the dimension of the smallest block at the leaf level
 * @param tolerance     requested accuracy of result
 * @param sym_spd       specifiy whether matrix is generic (0), symmetric (1) or spd (2)
 * @retval logdet       \f$\log\mathrm{det}(1-M)\f$
 */
EXTERNC double hodlr_logdet(int dim, double (*callback)(int,int,void *), void *args, unsigned int nLeaf, double tolerance, int is_symmetric);

#endif
