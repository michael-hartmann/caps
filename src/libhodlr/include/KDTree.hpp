/*!
 \class KDTree
 
 \brief This class is for constructing a KDTree in any dimension.
 
 \note
 
 \author $Sivaram Ambikasaran$
 
 \version
 
 \date $October 22nd, 2013$
 
 Contact: siva.1985@gmail.com
 */

#ifndef __KDTREE_HPP__
#define __KDTREE_HPP__

#include "Eigen/Dense"

using namespace Eigen;

void mergeSortedLists(MatrixXd& list1, MatrixXd& list2, unsigned index, MatrixXd& final_list);
void mergeSort(MatrixXd& locations, unsigned index);
void getKDTreeSorted(MatrixXd& locations, unsigned index);

#endif /* defined(__KDTREE_HPP__) */
