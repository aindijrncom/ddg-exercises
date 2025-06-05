// Implement member functions for HodgeDecomposition class.
#include "hodge-decomposition.h"

#include "geometrycentral/numerical/linear_solvers.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HodgeDecomposition::HodgeDecomposition(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: build DEC operators
    this->hodge1 = geometry->buildHodgeStar1Form(); // placeholder
    this->hodge2 = geometry->buildHodgeStar2Form(); // placeholder
    this->d0 = geometry->buildExteriorDerivative0Form(); // placeholder
    this->d1 = geometry->buildExteriorDerivative1Form(); // placeholder

    // TODO: Build operator inverses.
    // Hint: Use the sparseInverseDiagonal() in utils/src/solvers.cpp to invert sparse diagonal matrices.
    this->hodge1Inv = sparseInverseDiagonal(hodge1); // placeholder
    this->hodge2Inv = sparseInverseDiagonal(hodge2); // placeholder
    this->d0T = d0.transpose();                      // placeholder
    this->d1T = d1.transpose();                      // placeholder

    // TODO: Construct 0-form Laplace matrix.
    // Shift matrix by a small constant (1e-8) to make it positive definite.
    this->A = geometry->laplaceMatrix(); // placeholder

    // TODO: Construct 2-form matrix.
    this->B = d1*hodge1Inv*d1T; // placeholder
}

/*
 * Compute the 0-form potential α by solving the system 𝛿dα = 𝛿ω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The exact component dα of ω.
 */
/*
 * 计算精确分量 dα，通过求解系统 δdα = δω。
 *
 * 输入：定义在网格边上的 1-form omega。
 * 返回：ω 的精确分量 dα。
 */
Vector<double> HodgeDecomposition::computeExactComponent(const Vector<double>& omega) const {

    // 计算 δω = d0T * hodge1 * omega
    Vector<double> deltaOmega = d0T * hodge1 * omega;

    // 解线性系统 A * α = δω，其中 A 是 0-form Laplace 矩阵
    SparseMatrix<double> AA = A;
    Vector<double> alpha = solvePositiveDefinite(AA, deltaOmega);

    // 计算精确分量 dα = d0 * α
    Vector<double> dAlpha = d0 * alpha;

    return dAlpha;

    // TODO
    return Vector<double>::Zero(1); // placeholder
}

/*
 * Compute the 2-form potential β by solving the system d𝛿β = dω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The coexact component 𝛿β of ω.
 */
/*
 * 计算协同精确分量 δβ，通过求解系统 dδβ = dω。
 *
 * 输入：定义在网格边上的 1-form omega。
 * 返回：ω 的协同精确分量 δβ。
 */
Vector<double> HodgeDecomposition::computeCoExactComponent(const Vector<double>& omega) const {
    // 计算 dω = d1 * omega
    Vector<double> dOmega = d1 * omega;

    // 解线性系统 B * β = dω，其中 B 是 2-form 矩阵（hodge2Inv）
    SparseMatrix<double> BB = B;
    Vector<double> beta = solveSquare(BB, dOmega);

    // 计算协同精确分量 δβ = d1T * hodge2Inv * β
    Vector<double> deltaBeta = d1T* (hodge2Inv * beta);

    return deltaBeta;
}

/*
 * Compute the harmonic component γ = ω - dα - 𝛿β of ω.
 *
 * Input: A primal 1-form <omega> on the edges of the input mesh, the exact component <dAlpha> of ω, and the coexact
 * component <deltaBeta> of ω.
 * Returns: The coexact component 𝛿β of ω.
 */
/*
 * 计算调和分量 γ = ω - dα - δβ。
 *
 * 输入：定义在网格边上的 1-form omega，精确分量 dα 和协同精确分量 δβ。
 * 返回：ω 的调和分量 γ。
 */
Vector<double> HodgeDecomposition::computeHarmonicComponent(const Vector<double>& omega, const Vector<double>& dAlpha,
                                                            const Vector<double>& deltaBeta) const {
    // 调和分量是 ω 减去精确分量和协同精确分量
    Vector<double> gamma = omega - dAlpha - deltaBeta;

    return gamma;
}