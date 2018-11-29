/*!
 \brief This function is for the construction of KDTree.
 
 \note
 
 \author $Sivaram Ambikasaran$
 
 \version
 
 \date $Date: November 14th, 2013$
 */

#include "KDTree.hpp"
#include "Eigen/Dense"

using namespace Eigen;

void mergeSortedLists(MatrixXd& list1, MatrixXd& list2, unsigned index, MatrixXd& final_list) 
{
    unsigned N1 = list1.rows();
    unsigned N2 = list2.rows();
    unsigned j1 = 0;
    unsigned j2 = 0;
    unsigned j  = 0;

    while (j1 < N1 && j2 < N2) 
    {
        if (list1(j1,index) < list2(j2,index)) 
        {
            final_list.row(j) = list1.row(j1);
            ++j1;
        }
    
        else 
        {
            final_list.row(j)    =   list2.row(j2);
            ++j2;
        }
    
        ++j;
    }

    while (j1<N1) 
    {
        final_list.row(j) = list1.row(j1);
        ++j1;
        ++j;
    }

    while (j2<N2) 
    {
        final_list.row(j) = list2.row(j2);
        ++j2;
        ++j;
    }
}

void mergeSort(MatrixXd& locations, unsigned index) 
{
    unsigned N = locations.rows();
    
    if(N == 1) 
    {
        return;
    }
    
    else 
    {
        // Number of points in the left cluster.
        unsigned N_left = N / 2;
        
        // Number of points in the right cluster.
        unsigned N_right = N - N_left;
        
        // Dimension of the space.
        unsigned dims = locations.cols();
        
        // Left locations.
        MatrixXd left_locations = locations.block(0,0,N_left,dims);
        
        // Right locations.
        MatrixXd right_locations = locations.block(N_left,0,N_right,dims);
        
        // Mergesort for the left.
        mergeSort(left_locations, index);
        
        // Mergesort for the right.
        mergeSort(right_locations, index);
        
        // Merge the sorted left and right lists.
        mergeSortedLists(left_locations, right_locations, index, locations);
    }
}

void getKDTreeSorted(MatrixXd& locations, unsigned index) 
{
    // Get the total number of points
    unsigned N = locations.rows();
    
    if(N == 1) 
    {
        return;
    }
    
    else 
    {
        // Number of points in the left cluster.
        unsigned N_left = N / 2;

        // Number of points in the right cluster.
        unsigned N_right = N - N_left;

        // Dimension of the space.
        unsigned dims = locations.cols();

        // Merge sort on the input locations based on the coordinate index%2.
        mergeSort(locations, index % dims);

        /// Obtain the left and right locations.
        MatrixXd left_locations  = locations.block(0, 0, N_left, dims);
        MatrixXd right_locations = locations.block(N_left, 0, N_right, dims);

        /// Sort the left and right locations based on a KDTree.
        getKDTreeSorted(left_locations, index + 1);
        getKDTreeSorted(right_locations, index + 1);

        /// Output the locations.
        locations.block(0, 0, N_left, dims)       = left_locations;
        locations.block(N_left, 0, N_right, dims) = right_locations;
    }
}
