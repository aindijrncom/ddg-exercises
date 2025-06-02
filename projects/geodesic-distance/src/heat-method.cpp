// Implement member functions HeatMethod class.
#include "heat-method.h"

#include "geometrycentral/numerical/linear_solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;

    // TODO: Build Laplace and flow matrices.
    // Compute mean edge length
    double h = geometry->meanEdgeLength();

    // Compute time step t = h^2 as suggested in the paper
    double t = h * h;

    // Build mass matrix M (diagonal matrix of Voronoi areas)
    SparseMatrix<double> M = geometry->massMatrix();

    // Build cotangent Laplacian matrix L
    size_t n = mesh->nVertices();
    SparseMatrix<double> L(n, n);
    {
        std::vector<Eigen::Triplet<double>> triplets;
        Eigen::VectorXd diagonalSums = Eigen::VectorXd::Zero(n);
       
        // 遍历所有边计算cotangent权重
        for (Edge e : mesh->edges()) {
            Halfedge he = e.halfedge();
            double w = 0.5 * (geometry->halfedgeCotanWeight(he) + geometry->halfedgeCotanWeight(he.twin()));
       
            Vertex v1 = he.vertex();
            Vertex v2 = he.twin().vertex();
            int i = v1.getIndex();
            int j = v2.getIndex();
       
            // 添加非对角元素 (正值)
            triplets.emplace_back(i, j, w);
            triplets.emplace_back(j, i, w);
       
            // 累加对角元素 (负值)
            diagonalSums[i] += w;
            diagonalSums[j] += w;
        }
       
        // 添加对角元素 (负的权重和)
        for (Vertex v : mesh->vertices()) 
        {
            int i = v.getIndex();
            triplets.emplace_back(i, i, -diagonalSums[i]);
        }

        L.setFromTriplets(triplets.begin(), triplets.end());
        L.makeCompressed();
        this->F = L;
    }
    //SparseMatrix<double> L = geometry->laplaceMatrix();

    // Build A = L for the Poisson equation
    this->A = L;

    // Build F =  M - t * L for the heat equation
    L = geometry->laplaceMatrix();
    this->F = M - t * L;

    // Note: core/geometry.cpp has meanEdgeLength() function
    //this->A = identityMatrix<double>(1); // placeholder
    //this->F = identityMatrix<double>(1); // placeholder
}

    // 计算标量场 u 在每个面上的梯度
Vector3 HeatMethod::faceGradient(const Vector<double>& u, Face f) const {
    // 获取面的顶点
    std::vector<Vertex> vertices;
    for (Vertex v : f.adjacentVertices()) {
        vertices.push_back(v);
    }
    assert(vertices.size() == 3); // 假设为三角形面

    // 获取顶点上的标量值
    double u0 = u[vertices[0].getIndex()];
    double u1 = u[vertices[1].getIndex()];
    double u2 = u[vertices[2].getIndex()];

    // 计算面法向量
    Vector3 N = geometry->faceNormal(f).normalize();


    // 计算边向量（逆时针）
    Vector3 e0 = geometry->inputVertexPositions[vertices[2]] - geometry->inputVertexPositions[vertices[1]];
    Vector3 e1 = geometry->inputVertexPositions[vertices[0]] - geometry->inputVertexPositions[vertices[2]];
    Vector3 e2 = geometry->inputVertexPositions[vertices[1]] - geometry->inputVertexPositions[vertices[0]];

    // 计算梯度
    geometry->requireFaceAreas();
    Vector3 grad_u = (u0 * cross(N, e0) + u1 * cross(N, e1) + u2 * cross(N, e2)) / (2 * geometry->faceAreas[f]);
    // Vector3 grad_u = (u0 * cross(N, e0)) / (2 * geometry->faceAreas[f]);

    return grad_u.normalize();
}

// 计算矢量场 X 在每个顶点上的集成散度
Vector<double> HeatMethod::integratedDivergence(const FaceData<Vector3>& X) const {
    Vector<double> div_X = Vector<double>::Zero(mesh->nVertices());

    for (Face f : mesh->faces()) {
        // 获取面的顶点
        std::vector<Vertex> vertices;
        for (Vertex v : f.adjacentVertices()) {
            vertices.push_back(v);
        }
        assert(vertices.size() == 3); // 假设为三角形面

        // 获取面上的矢量场
        Vector3 X_f = X[f];

        // 计算边向量
        Vector3 e0 = geometry->inputVertexPositions[vertices[2]] - geometry->inputVertexPositions[vertices[1]];
        Vector3 e1 = geometry->inputVertexPositions[vertices[0]] - geometry->inputVertexPositions[vertices[2]];
        Vector3 e2 = geometry->inputVertexPositions[vertices[1]] - geometry->inputVertexPositions[vertices[0]];

        // 计算对边顶点的角度余切
        // 获取顶点位置
        Vector3 p0 = geometry->inputVertexPositions[vertices[0]];
        Vector3 p1 = geometry->inputVertexPositions[vertices[1]];
        Vector3 p2 = geometry->inputVertexPositions[vertices[2]];

        // 计算 cot_alpha0（顶点 v0 的角度余切）
        Vector3 e01 = p1 - p0; // 从 v0 到 v1 的向量
        Vector3 e02 = p2 - p0; // 从 v0 到 v2 的向量
        double cot_alpha0 = dot(e01, e02) / norm(cross(e01, e02));

        // 计算 cot_alpha1（顶点 v1 的角度余切）
        Vector3 e10 = p0 - p1; // 从 v1 到 v0 的向量
        Vector3 e12 = p2 - p1; // 从 v1 到 v2 的向量
        double cot_alpha1 = dot(e10, e12) / norm(cross(e10, e12));

        // 计算 cot_alpha2（顶点 v2 的角度余切）
        Vector3 e20 = p0 - p2; // 从 v2 到 v0 的向量
        Vector3 e21 = p1 - p2; // 从 v2 到 v1 的向量
        double cot_alpha2 = dot(e20, e21) / norm(cross(e20, e21));

        // 计算每个顶点的散度贡献
        div_X[vertices[0].getIndex()] += 0.5 * (cot_alpha1 * dot(e1, X_f) + cot_alpha2 * dot(e2, X_f));
        div_X[vertices[1].getIndex()] += 0.5 * (cot_alpha2 * dot(e2, X_f) + cot_alpha0 * dot(e0, X_f));
        div_X[vertices[2].getIndex()] += 0.5 * (cot_alpha0 * dot(e0, X_f) + cot_alpha1 * dot(e1, X_f));
    }

    return div_X;
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {

    // TODO
    // Compute the normalized negative gradient for each face
    FaceData<Vector3> X(*mesh);
    for (Face f : mesh->faces()) {
        // Compute gradient of u on face f
        Vector3 grad_u = faceGradient(u, f);

        // Compute magnitude
        double mag = norm(grad_u);

        // Avoid division by zero; set X to zero vector if gradient is negligible
        if (mag > 1e-10) {
            X[f] = (-grad_u / mag).normalize();
        } else {
            X[f] = {0, 0, 0};
        }
    }

    return X;

    return FaceData<Vector3>(*mesh, {0, 0, 0}); // placeholder
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {

    // TODO
    // Compute integrated divergence at each vertex
    Vector<double> div_X =integratedDivergence(X);
    return div_X;
    return Vector<double>::Zero(1); // placeholder
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {

    // TODO
    // Step 1: Solve heat equation (M - t L) u = M delta
     Eigen::SimplicialLLT<SparseMatrix<double>> llt(F);
    SparseMatrix<double> M = geometry->massMatrix();
    Vector<double> b = M * delta;
    // 求解 Ax = b
    Vector<double> u = llt.solve(b);

    std::vector<double> val_vec1;
    val_vec1.reserve(u.size()); // 预分配内存以提高性能
    for (Eigen::Index i = 0; i < u.size(); ++i) {
        if (b[i] != 0) {
            int a = 0;
            b[i] += 10;
        }
        val_vec1.push_back(b[i]);
    }
      
    std::vector<double> val_vec;
    val_vec.reserve(u.size()); // 预分配内存以提高性能
    for (Eigen::Index i = 0; i < u.size(); ++i) {
        val_vec.push_back(u(i));
    }

    // Step 2: Compute vector field X = -∇u / |∇u|
    FaceData<Vector3> X = computeVectorField(u);

    // Step 3: Compute divergence ∇.X
    Vector<double> div_X = computeDivergence(X);

    // Step 4: Solve Poisson equation L φ = div_X
    SparseMatrix<double> L = geometry->laplaceMatrix();
    int r = L.rows();
    int c = L.cols();

    Vector<double> phi = solveSquare(L, div_X);

    // Step 5: Shift φ so that the smallest distance is zero
    //double min_phi = phi.minCoeff();
    //phi.array() -= min_phi;
    this->subtractMinimumDistance(phi);
    return phi;

    //Vector<double> phi = Vector<double>::Zero(delta.rows());

    //// Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    //this->subtractMinimumDistance(phi);

    //return phi;
}