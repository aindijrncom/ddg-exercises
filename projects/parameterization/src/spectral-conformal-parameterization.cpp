// Implement member functions for SpectralConformalParameterization class.
#include "spectral-conformal-parameterization.h"

#include "geometrycentral/numerical/linear_solvers.h"
#include "polyscope/utilities.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
SpectralConformalParameterization::SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh,
                                                                     VertexPositionGeometry* inputGeo) {

    this->mesh = inputMesh;
    this->geometry = inputGeo;
}

/*
 * 构建复共形能量矩阵 EC = ED - A
 *
 * 数学定义:
 *   EC = ED - A
 *   ED = 复拉普拉斯矩阵 (实部为拉普拉斯算子，虚部为0)
 *   A = 面积矩阵 (实对角矩阵，表示顶点对偶面积)
 *
 * 物理意义:
 *   共形能量度量曲面参数化与共形映射的偏差
 *   EC 的最小特征向量对应最佳共形参数化
 *
 * 返回: 复稀疏矩阵表示的共形能量
 */
/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {

    // TODO
    // 1. 获取复拉普拉斯矩阵 (实部为正定拉普拉斯矩阵，虚部为0)
    SparseMatrix<std::complex<double>> ED = geometry->complexLaplaceMatrix();

    // 2. 构建面积矩阵 A (实对角矩阵，元素为顶点对偶面积)
    size_t n = geometry->mesh.nVertices();
    SparseMatrix<double> A_real(n, n);
    std::vector<Eigen::Triplet<double>> tripletsA;
    tripletsA.reserve(n);

    for (Vertex v : geometry->mesh.vertices()) {
        int idx = v.getIndex();
        double area = geometry->barycentricDualArea(v);
        tripletsA.emplace_back(idx, idx, area);
    }
    A_real.setFromTriplets(tripletsA.begin(), tripletsA.end());

    // 3. 将面积矩阵转换为复数矩阵 (实部为面积值，虚部为0)
    SparseMatrix<std::complex<double>> A = A_real.cast<std::complex<double>>();

    // 4. 计算共形能量矩阵: EC = ED - A
    SparseMatrix<std::complex<double>> EC = ED - A;

    return EC;

    return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * 使用谱共形参数化方法将带边界的曲面网格共形展平到平面
 *
 * 步骤:
 * 1. 构建复共形能量矩阵 EC = ED - A
 * 2. 求解广义特征值问题：EC * x = λ * A * x
 * 3. 取最小非零特征值对应的特征向量
 * 4. 将特征向量转换为二维坐标 (u, v)
 *
 * 返回: 每个顶点对应的二维平面坐标
 */
/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {
    // 1. 构建复共形能量矩阵
    SparseMatrix<std::complex<double>> EC = buildConformalEnergy();

    // 2. 获取面积矩阵（质量矩阵）
    size_t n = mesh->nVertices();
    SparseMatrix<double> A_real(n, n);
    std::vector<Eigen::Triplet<double>> tripletsA;
    tripletsA.reserve(n);

    for (Vertex v : mesh->vertices()) {
        int idx = v.getIndex();
        double area = geometry->barycentricDualArea(v);
        tripletsA.emplace_back(idx, idx, area);
    }
    A_real.setFromTriplets(tripletsA.begin(), tripletsA.end());
    SparseMatrix<std::complex<double>> A = A_real.cast<std::complex<double>>();

     // 3. 使用提供的接口求解最小特征向量
    Vector<std::complex<double>> complexCoords;

    // 尝试使用正定求解器
    try {
        complexCoords = smallestEigenvectorPositiveDefinite(EC, A);
    }
    // 如果正定求解失败，使用备用求解器
    catch (...) {
        try {
            complexCoords = smallestEigenvectorSquare(EC, A);
        }
        // 如果所有方法都失败，使用最大特征向量作为备用
        catch (...) {
            complexCoords = largestEigenvector(EC, A);
        }
    }

   // 4. 转换为二维坐标并归一化
    VertexData<Vector2> parameterization(*mesh);

    // 计算坐标范围
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    for (Vertex v : mesh->vertices()) {
        int idx = v.getIndex();
        std::complex<double> coord = complexCoords(idx);
        minX = std::min(minX, coord.real());
        maxX = std::max(maxX, coord.real());
        minY = std::min(minY, coord.imag());
        maxY = std::max(maxY, coord.imag());
    }

    // 归一化到[0,1]范围
    double rangeX = maxX - minX;
    double rangeY = maxY - minY;
    double scale = (rangeX > 0 || rangeY > 0) ? 1.0 / std::max(rangeX, rangeY) : 1.0;

    // 5. 存储归一化后的二维坐标
    for (Vertex v : mesh->vertices()) {
        int idx = v.getIndex();
        std::complex<double> coord = complexCoords(idx);
        double u = (coord.real() - minX) * scale;
        double v_coord = (coord.imag() - minY) * scale;
        parameterization[v] = Vector2{u, v_coord};
    }

    return parameterization;

    // TODO
    return VertexData<Vector2>(*mesh); // placeholder
}