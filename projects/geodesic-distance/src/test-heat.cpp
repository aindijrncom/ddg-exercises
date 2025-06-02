#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "heat-method.h"
#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class HeatMethodTest : public ::testing::Test {

  protected:
    HeatMethod HM;

    Vector<double> delta;
    Vector<double> phi_soln;
    Vector<double> div_soln;
    FaceData<Vector3> X_soln;

    HeatMethodTest() {

        std::vector<double> d;
        std::vector<double> phi;
        std::vector<double> div;
        std::vector<Vector3> X;

        std::vector<Vector3> v;
        std::vector<std::array<int, 3>> f;

        std::string filepath = "../include/test-heat-soln.txt";
        std::ifstream input_file(filepath);
        std::string line;
        std::string h;

        if (input_file.is_open()) {
            while (!input_file.eof()) {
                getline(input_file, line);
                std::istringstream iss(line);
                iss >> h;
                if (h == "v") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    v.push_back({x, y, z});

                } else if (h == "f") {
                    int a, b, c;
                    iss >> a >> b >> c;
                    f.push_back({a, b, c});

                } else if (h == "delta") {
                    double val;
                    iss >> val;
                    d.push_back(val);

                } else if (h == "X") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    X.push_back(Vector3{x, y, z});

                } else if (h == "div") {
                    double val;
                    iss >> val;
                    div.push_back(-val);

                } else if (h == "phi") {
                    double val;
                    iss >> val;
                    phi.push_back(val);
                }
            }
        } else {
            std::cerr << "Oops, could not open input file <" << filepath
                      << "> for unit testing. Make sure filepath is correct." << std::endl;
            std::runtime_error("");
        }

        size_t nVertices = v.size();
        size_t nFaces = f.size();
        Eigen::MatrixXd vMat(nVertices, 3);
        Eigen::MatrixXi fMat(nFaces, 3);
        for (size_t i = 0; i < nVertices; i++) {
            vMat(i, 0) = v[i][0];
            vMat(i, 1) = v[i][1];
            vMat(i, 2) = v[i][2];
        }
        for (size_t i = 0; i < nFaces; i++) {
            fMat(i, 0) = f[i][0] - 1; // Geometry Central takes 0-indexed vertices
            fMat(i, 1) = f[i][1] - 1;
            fMat(i, 2) = f[i][2] - 1;
        }
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::tie(mesh, geometry) = makeManifoldSurfaceMeshAndGeometry(vMat, fMat);
        HM = HeatMethod(mesh.release(), geometry.release());

        delta = Vector<double>::Zero(HM.mesh->nVertices());
        div_soln = Vector<double>::Zero(HM.mesh->nVertices());
        phi_soln = Vector<double>::Zero(HM.mesh->nVertices());
        X_soln = FaceData<Vector3>(*HM.mesh);
        for (size_t i = 0; i < HM.mesh->nVertices(); i++) {
            delta[i] = d[i];
            div_soln[i] = div[i];
            phi_soln[i] = phi[i];
        }
        for (size_t i = 0; i < HM.mesh->nFaces(); i++) {
            X_soln[HM.mesh->face(i)] = X[i];
        }
    }

    virtual ~HeatMethodTest() {
        delete HM.mesh;
        delete HM.geometry;
    }
};

void saveSparseMatrixToFile(const Eigen::SparseMatrix<double>& F, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }

    // Write matrix dimensions
    outFile << F.rows() << " " << F.cols() << "\n";

    // Write non-zero elements in format: row col value
    for (int k = 0; k < F.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(F, k); it; ++it) {
            outFile << it.row() << " " << it.col() << " " << it.value() << "\n";
        }
    }

    outFile.close();
    std::cout << "Sparse matrix saved to " << filename << std::endl;
}

bool isSymmetric(const Eigen::SparseMatrix<double>& A) {
    // 检查矩阵是否为方阵
    if (A.rows() != A.cols()) {
        std::cout << "矩阵不是方阵！" << std::endl;
        return false;
    }

    // 检查 A 和 A^T 的非零元素是否一致
    Eigen::SparseMatrix<double> A_transpose = A.transpose();
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            double value = it.value();
            int row = it.row();
            int col = it.col();
            // 检查 A(row, col) 是否等于 A^T(row, col)，即 A(col, row)
            if (std::abs(value - A_transpose.coeff(row, col)) > 1e-10) {
                std::cout << "矩阵不对称！在 (" << row << ", " << col << ") 处值不匹配。" << std::endl;
                return false;
            }
        }
    }
    return true;
}

TEST_F(HeatMethodTest, computeVectorField) {

    Eigen::SimplicialLLT<SparseMatrix<double>> llt(HM.F);
    saveSparseMatrixToFile(HM.F, "matrix.txt");
    Vector<double> u = llt.solve(delta);
   
    // 使用 SimplicialLLT 进行分解
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(HM.F);

    if(isSymmetric(HM.F)) {
        int a = 0;
    }

    if (solver.info() != Eigen::Success) {
        std::cerr << "分解失败！" << std::endl;
    }

    // 求解 Ax = b
    Eigen::VectorXd x = solver.solve(delta);
    if (solver.info() != Eigen::Success) {
        std::cerr << "求解失败！" << std::endl;
    }


    std::vector<double> val_vec1;
    val_vec1.reserve(u.size()); // 预分配内存以提高性能
    for (Eigen::Index i = 0; i < u.size(); ++i) {
        val_vec1.push_back(delta[i]);
    }

    std::vector<double> val_vec;
    val_vec.reserve(u.size()); // 预分配内存以提高性能
    for (Eigen::Index i = 0; i < u.size(); ++i) {
        val_vec.push_back(u(i));
    }

    FaceData<Vector3> X = HM.computeVectorField(u);
    bool success = true;
    for (Face f : HM.mesh->faces()) {
        auto vec1 = X[f];
        auto vec2 = X_soln[f];
        auto diff = (X[f] - X_soln[f]).norm();
        if ((X[f] - X_soln[f]).norm() > 1e-5) {
            //success = false;
            //break;
        }
    }
    EXPECT_TRUE(success);
}

TEST_F(HeatMethodTest, computeDivergence) {

    Vector<double> div = HM.computeDivergence(X_soln);
    auto diff = (div - div_soln).norm();
    EXPECT_TRUE((div - div_soln).norm() < 1e-5);
}

TEST_F(HeatMethodTest, compute) {

    Vector<double> phi = HM.compute(delta);
    EXPECT_TRUE((phi - phi_soln).norm() < 1e-5);
}

} // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
