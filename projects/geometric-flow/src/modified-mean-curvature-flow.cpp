// Implement member functions for ModifiedMeanCurvatureFlow class.
#include "modified-mean-curvature-flow.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ModifiedMeanCurvatureFlow::ModifiedMeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry, and A (Laplace matrix)
    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: build the Laplace matrix
    this->A = geometry->laplaceMatrix();
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> ModifiedMeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {
    // TODO
    // 从几何体获取负半定拉普拉斯矩阵 (cotangent权重)
    SparseMatrix<double> L_neg = this->A;

    // 计算流算子: A = M - h * L
    // 注意: L_neg 是负半定矩阵，所以 -h * L_neg 是正定贡献
    SparseMatrix<double> A = M - h * L_neg;
    return A;

    return identityMatrix<double>(1); // placeholder
}
