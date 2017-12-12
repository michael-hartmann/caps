/*!
 \class HODLR_Tree

 \brief This class is the main class for the HODLR tree.

 \note

 \author $Sivaram Ambikasaran$

 \version

 \date $November 8th, 2013$

 Contact: siva.1985@gmail.com
 */

#ifndef __HODLR_TREE_HPP__
#define __HODLR_TREE_HPP__

#include <iostream>
#include <cmath>
#include <vector>
#include "HODLR_Node.hpp"

using std::cout;
using std::endl;
using std::vector;

template <typename MatrixType>
class HODLR_Tree {
private:
    /** Variables of interest for the entire tree */
    friend class HODLR_Node<MatrixType>;
    MatrixType* kernel;
    HODLR_Node<MatrixType>* root;
    int N;
    int nLeaf;
    double lowRankTolerance;
    double determinant;
    VectorXd diagonal;
    int nLevels;
    int nNodes;
    vector< vector<int*> > ranks;
    vector< vector<int> > nCalls_To_Get_Matrix_Entry;
    char s; //  s = 's' means symmetric, else non-symmetric
        /*!
         Creates the tree.
         */
    void create_Tree(HODLR_Node<MatrixType>*& node)
    {
        ++nNodes;

        if(node->nSize <= nLeaf)
            node->isLeaf = true;
        else
        {
            node->isLeaf = false;
            int n0 = node->nSize/2;
            int n1 = node->nSize-n0;
            node->child[0] = new HODLR_Node<MatrixType>(kernel, node->levelNumber+1, 0, 2*node->levelBasedNodeNumber, node->nStart, n0);
            node->child[0]->parent = node;
            node->child[1] = new HODLR_Node<MatrixType>(kernel, node->levelNumber+1, 1, 2*node->levelBasedNodeNumber+1, node->nStart+n0, n1);
            node->child[1]->parent = node;
            create_Tree(node->child[0]);
            create_Tree(node->child[1]);
        }
    };

        /*!
         Assembles the matrix.
         */
    void assemble_Matrix(HODLR_Node<MatrixType>*& node)
    {
        if (node)
        {
            node->assemble_Matrices(lowRankTolerance, diagonal, s);
            assemble_Matrix(node->child[0]);
            assemble_Matrix(node->child[1]);
        }
    };

        /*!
         Obtain the rank of all the off-diagonal blocks.
         */
        void obtain_All_Ranks(HODLR_Node<MatrixType>* node)
        {
            if(node)
            {
                if (ranks.size()<node->levelNumber+1)
                {
                    vector<int*> temp;
                    ranks.push_back(temp);
                    vector<int> temp1;
                    nCalls_To_Get_Matrix_Entry.push_back(temp1);
                }

                int* temp = new int[2];
                temp[0]   = node->nRank[0];
                temp[1]   = node->nRank[1];
                ranks[node->levelNumber].push_back(temp);
                int temp1 = node->nSize*(temp[0]+temp[1])-temp[0]*temp[0]-temp[1]*temp[1];
                nCalls_To_Get_Matrix_Entry[node->levelNumber].push_back(temp1);
                obtain_All_Ranks(node->child[0]);
                obtain_All_Ranks(node->child[1]);
            }
        }

        /*!
         Matrix matrix product.
         */
    void matMatProduct(HODLR_Node<MatrixType>*& node, MatrixXd& x, MatrixXd& b)
    {
        if(node)
        {
            node->matrix_Matrix_Product(x, b);
            matMatProduct(node->child[0], x, b);
            matMatProduct(node->child[1], x, b);
        }
    };

        /*!
         Computes the factor of the HODLR matrix.
         */
    void compute_Factor(HODLR_Node<MatrixType>*& node)
    {
        if(node)
        {
            for(int k = 0; k < 2; k++)
                compute_Factor(node->child[k]);

            node->compute_K(s);
            node->compute_Inverse();
            HODLR_Node<MatrixType>*& mynode = node->parent;
            int number = node->nodeNumber;
            int mStart = node->nStart;
            while (mynode)
            {
                node->apply_Inverse(mynode->Uinverse[number], mStart);
                number = mynode->nodeNumber;
                mStart = mynode->nStart;
                mynode = mynode->parent;
            }
        }
    };

        /*!
         Set the matrices for inversion.
         */
    void set_Matrices_For_Inversion(HODLR_Node<MatrixType>*& node)
    {
        if (node)
        {
            node->set_UV_Inversion();
            set_Matrices_For_Inversion(node->child[0]);
            set_Matrices_For_Inversion(node->child[1]);
        }
    };

        /*!
         Computes the determinant of the HODLR matrix.
         */
    void compute_Determinant(HODLR_Node<MatrixType>*& node)
    {
        if(node)
        {
            for(int k = 0; k < 2; ++k)
            {
                compute_Determinant(node->child[k]);
                node->compute_Determinant();
            }
            determinant = determinant+node->determinant;
        }
    };

public:
        /*!
         Constructor for the HODLR tree.
         */
    HODLR_Tree(MatrixType* kernel, int N, int nLeaf)
    {
        this->kernel  = kernel;
        this->N       = N;
        this->nLeaf   = nLeaf;
        this->nNodes  = 0;
        this->nLevels = log2(nNodes);
        root          = new HODLR_Node<MatrixType>(kernel, 0, 0, 0, 0, N);
        create_Tree(root);
    };

        /*!
         Assembles the matrix.
         */
    void assemble_Matrix(VectorXd& diagonal, double lowRankTolerance, char s)
    {
        this->lowRankTolerance = lowRankTolerance;
        this->diagonal         = diagonal;
        this->s                = s;
        assemble_Matrix(root);
    };

        /*!
         Matrix matrix product.
         */
    void matMatProduct(MatrixXd& x, MatrixXd& b)
    {
        b = MatrixXd::Zero(N, x.cols());
        matMatProduct(root, x, b);
    };

        /*!
         Computes the factor of the HODLR matrix.
         */
    void compute_Factor()
    {
        set_Matrices_For_Inversion(root);
        compute_Factor(root);
    }

        /*!
         Solves the HODLR linear system.
         */
    void solve(MatrixXd& b, MatrixXd& x)
    {
        x = b;
        solve(root, x);
    };

        /*!
         Computes the determinant.
         */
    double compute_Determinant(double& determinant)
    {
        this->determinant = 0;
        compute_Determinant(root);
        determinant = this->determinant;
        return determinant;
    };

    ~HODLR_Tree() {
        delete root;
    }
};

#endif /* defined(__HODLR_TREE_HPP__) */
