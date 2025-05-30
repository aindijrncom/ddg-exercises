// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"

#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}


/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    // TODO
    // 从几何体获取负半定拉普拉斯矩阵 (cotangent权重)
    SparseMatrix<double> L_neg = geometry->laplaceMatrix();

    // 计算流算子: A = M - h * L
    // 注意: L_neg 是负半定矩阵，所以 -h * L_neg 是正定贡献
    SparseMatrix<double> A = M - h * L_neg;
    return A;
    

    return identityMatrix<double>(1); // placeholder
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    // Note: Update positions via geometry->inputVertexPositions
    // 在平均曲率流模拟中
    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> A = buildFlowOperator(M, h);

    // 2. 准备右侧向量 (M * currentPositions)
    size_t n = mesh->nVertices();
    Eigen::VectorXd rhsX(n), rhsY(n), rhsZ(n);

    // 遍历所有顶点填充右侧向量
    for (Vertex v : mesh->vertices()) {
        int idx = v.getIndex();
        Vector3 pos = geometry->inputVertexPositions[v];
        double area = M.coeff(idx, idx); // 质量矩阵对角元素

        rhsX(idx) = area * pos.x;
        rhsY(idx) = area * pos.y;
        rhsZ(idx) = area * pos.z;
    }

    // 3. 求解线性系统 (对每个坐标分量)
    // 求解线性系统: A x^{n+1} = M x^n
    // Build the solver
    Solver<double> solver(A);
    Eigen::VectorXd newX = solver.solve(rhsX);
    Eigen::VectorXd newY = solver.solve(rhsY);
    Eigen::VectorXd newZ = solver.solve(rhsZ);

    // 4. 更新顶点位置
    for (Vertex v : mesh->vertices()) {
        int idx = v.getIndex();
        geometry->inputVertexPositions[v] = Vector3{newX(idx), newY(idx), newZ(idx)};
    }

    return;

    for (Vertex v : mesh->vertices()) {
        geometry->inputVertexPositions[v] = geometry->inputVertexPositions[v]; // placeholder
    }
}