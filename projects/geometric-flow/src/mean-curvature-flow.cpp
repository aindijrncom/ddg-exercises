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
    // �Ӽ������ȡ���붨������˹���� (cotangentȨ��)
    SparseMatrix<double> L_neg = geometry->laplaceMatrix();

    // ����������: A = M - h * L
    // ע��: L_neg �Ǹ��붨�������� -h * L_neg ����������
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
    // ��ƽ��������ģ����
    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> A = buildFlowOperator(M, h);

    // 2. ׼���Ҳ����� (M * currentPositions)
    size_t n = mesh->nVertices();
    Eigen::VectorXd rhsX(n), rhsY(n), rhsZ(n);

    // �������ж�������Ҳ�����
    for (Vertex v : mesh->vertices()) {
        int idx = v.getIndex();
        Vector3 pos = geometry->inputVertexPositions[v];
        double area = M.coeff(idx, idx); // ��������Խ�Ԫ��

        rhsX(idx) = area * pos.x;
        rhsY(idx) = area * pos.y;
        rhsZ(idx) = area * pos.z;
    }

    // 3. �������ϵͳ (��ÿ���������)
    // �������ϵͳ: A x^{n+1} = M x^n
    // Build the solver
    Solver<double> solver(A);
    Eigen::VectorXd newX = solver.solve(rhsX);
    Eigen::VectorXd newY = solver.solve(rhsY);
    Eigen::VectorXd newZ = solver.solve(rhsZ);

    // 4. ���¶���λ��
    for (Vertex v : mesh->vertices()) {
        int idx = v.getIndex();
        geometry->inputVertexPositions[v] = Vector3{newX(idx), newY(idx), newZ(idx)};
    }

    return;

    for (Vertex v : mesh->vertices()) {
        geometry->inputVertexPositions[v] = geometry->inputVertexPositions[v]; // placeholder
    }
}