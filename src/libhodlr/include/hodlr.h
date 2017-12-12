#ifndef HODLR_H
#define HODLR_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/* @brief Calculate log(det(Id-M)) using HODLR approach
 *
 * @param dim           dimension of matrix M
 * @param callback      function that returns matrix elements of M
 * @param args          pointer that is passed to callback
 * @param nLeaf         ???
 * @param tolerance     requested accuracy of result
 * @param is_symmetric  matrix is symmetric (1) or not symmetric (0)
 * @retval log(det(Id-M))
 */
EXTERNC double hodlr_logdet(int dim, double (*callback)(int,int,void *), void *args, unsigned int nLeaf, double tolerance, int is_symmetric);

/* @brief Calculate log(det(Id-M)) using HODLR approach
 *
 * @param dim           dimension of matrix M
 * @param callback      function that returns matrix elements of M
 * @param args          pointer that is passed to callback
 * @param diagonal      diagonal elements of M
 * @param nLeaf         ???
 * @param tolerance     requested accuracy of result
 * @param is_symmetric  matrix is symmetric (1) or not symmetric (0)
 * @retval log(det(Id-M))
 */
EXTERNC double hodlr_logdet_diagonal(int dim, double (*callback)(int,int,void *), void *args, double *diagonal, unsigned int nLeaf, double tolerance, int is_symmetric);

#endif
