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
 * ������������������ EC = ED - A
 *
 * ��ѧ����:
 *   EC = ED - A
 *   ED = ��������˹���� (ʵ��Ϊ������˹���ӣ��鲿Ϊ0)
 *   A = ������� (ʵ�ԽǾ��󣬱�ʾ�����ż���)
 *
 * ��������:
 *   ����������������������빲��ӳ���ƫ��
 *   EC ����С����������Ӧ��ѹ��β�����
 *
 * ����: ��ϡ������ʾ�Ĺ�������
 */
/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {

    // TODO
    // 1. ��ȡ��������˹���� (ʵ��Ϊ����������˹�����鲿Ϊ0)
    SparseMatrix<std::complex<double>> ED = geometry->complexLaplaceMatrix();

    // 2. ����������� A (ʵ�ԽǾ���Ԫ��Ϊ�����ż���)
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

    // 3. ���������ת��Ϊ�������� (ʵ��Ϊ���ֵ���鲿Ϊ0)
    SparseMatrix<std::complex<double>> A = A_real.cast<std::complex<double>>();

    // 4. ���㹲����������: EC = ED - A
    SparseMatrix<std::complex<double>> EC = ED - A;

    return EC;

    return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * ʹ���׹��β��������������߽������������չƽ��ƽ��
 *
 * ����:
 * 1. ������������������ EC = ED - A
 * 2. ����������ֵ���⣺EC * x = �� * A * x
 * 3. ȡ��С��������ֵ��Ӧ����������
 * 4. ����������ת��Ϊ��ά���� (u, v)
 *
 * ����: ÿ�������Ӧ�Ķ�άƽ������
 */
/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {
    // 1. ������������������
    SparseMatrix<std::complex<double>> EC = buildConformalEnergy();

    // 2. ��ȡ���������������
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

     // 3. ʹ���ṩ�Ľӿ������С��������
    Vector<std::complex<double>> complexCoords;

    // ����ʹ�����������
    try {
        complexCoords = smallestEigenvectorPositiveDefinite(EC, A);
    }
    // ����������ʧ�ܣ�ʹ�ñ��������
    catch (...) {
        try {
            complexCoords = smallestEigenvectorSquare(EC, A);
        }
        // ������з�����ʧ�ܣ�ʹ���������������Ϊ����
        catch (...) {
            complexCoords = largestEigenvector(EC, A);
        }
    }

   // 4. ת��Ϊ��ά���겢��һ��
    VertexData<Vector2> parameterization(*mesh);

    // �������귶Χ
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

    // ��һ����[0,1]��Χ
    double rangeX = maxX - minX;
    double rangeY = maxY - minY;
    double scale = (rangeX > 0 || rangeY > 0) ? 1.0 / std::max(rangeX, rangeY) : 1.0;

    // 5. �洢��һ����Ķ�ά����
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