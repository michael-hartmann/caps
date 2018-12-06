#ifndef __HODLR_Node__
#define __HODLR_Node__

#include <Eigen/Dense>
#include "HODLR_Matrix.hpp"

class HODLR_Node 
{
friend class HODLR_Tree;

// All methods are declared as private since all usage happens from
// the friend class HODLR_Tree:
private:
    HODLR_Node(int node_number, int level_number, int local_number, 
               int n_start, int n_size, double tolerance);
    
    // Storing the information passed to constructor as attribute:
    int node_number, level_number, local_number;
    int n_start, n_size;
    // Storing the start locations and sizes for the children of the node:
    int c_start[2], c_size[2];
    // Tolerance for the computation carried out:
    double tolerance;

    // This stores the matrix directly(i.e at the leaf level)
    MatrixXd K;

    //  Variables and methods needed for HODLR solver
    MatrixXd U[2], V[2];
    MatrixXd U_factor[2], V_factor[2];
    Eigen::PartialPivLU<MatrixXd> K_factor;
    int rank[2];

    // Methods for Leaf Nodes:
    void assembleLeafNode(HODLR_Matrix* A);
    void matmatProductLeaf(MatrixXd x, MatrixXd& b);

    // Methods for Non-leaf Nodes:
    void assembleNonLeafNode(HODLR_Matrix* A, bool is_sym);
    void matmatProductNonLeaf(MatrixXd x, MatrixXd& b);

    // Method to print the parameters of the node(mainly used to debug)
    void printNodeDetails();
};

#endif /*__HODLR_Node__*/
