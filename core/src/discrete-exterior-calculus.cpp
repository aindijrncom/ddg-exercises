// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    // 获取网格中的顶点数量
    int nVertices = mesh.getVertexIndices().size();

    // 创建一个稀疏矩阵，初始为零矩阵
    SparseMatrix<double> hodgeStar0(nVertices, nVertices);

    // 设置对角元素为 1（根据约定，每个顶点的面积为 1）
    for (int i = 0; i < nVertices; ++i) {
        hodgeStar0.insert(i, i) = barycentricDualArea(mesh.vertex(i));
    }

    // 完成矩阵构造
    hodgeStar0.makeCompressed();

    return hodgeStar0;

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    // 获取网格的边数
    size_t nEdges = mesh.nEdges();

    // 创建三元组列表存储非零元素（对角线）
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nEdges);

    // 遍历所有边
    for (Edge e : mesh.edges()) {
        double cotSum = 0.0;
        int faceCount = 0;

        // 获取与边关联的两个半边
        Halfedge he1 = e.halfedge();
        Halfedge he2 = he1.twin();

        // 处理第一个三角形（如果存在）
        if (!he1.face().isBoundaryLoop()) {
            // 获取三角形中与边相对的顶点
            Vertex v_opposite = he1.next().next().vertex();

            // 计算边的两个向量
            Vector3 v0 = inputVertexPositions[he1.vertex()];
            Vector3 v1 = inputVertexPositions[he1.twin().vertex()];
            Vector3 v_opp = inputVertexPositions[v_opposite];

            // 计算从相对顶点到边端点的向量
            Vector3 vec1 = v0 - v_opp;
            Vector3 vec2 = v1 - v_opp;

            // 计算余切值：cosθ/sinθ = dot(v1,v2)/|v1×v2|
            double dotProd = dot(vec1, vec2);
            double crossNorm = norm(cross(vec1, vec2));
            double cot = dotProd / crossNorm;

            cotSum += cot;
            faceCount++;
        }

        // 处理第二个三角形（如果存在）
        if (!he2.face().isBoundaryLoop()) {
            // 获取三角形中与边相对的顶点
            Vertex v_opposite = he2.next().next().vertex();

            // 计算边的两个向量
            Vector3 v0 = inputVertexPositions[he2.vertex()];
            Vector3 v1 = inputVertexPositions[he2.twin().vertex()];
            Vector3 v_opp = inputVertexPositions[v_opposite];

            // 计算从相对顶点到边端点的向量
            Vector3 vec1 = v0 - v_opp;
            Vector3 vec2 = v1 - v_opp;

            // 计算余切值
            double dotProd = dot(vec1, vec2);
            double crossNorm = norm(cross(vec1, vec2));
            double cot = dotProd / crossNorm;

            cotSum += cot;
            faceCount++;
        }

        // 计算平均值：如果只有一个相邻三角形，则直接使用该值
        double value = (faceCount > 0) ? cotSum / 2.0 : 0.0;

        // 添加对角线元素
        triplets.emplace_back(e.getIndex(), e.getIndex(), value);
    }

    // 构建并填充稀疏矩阵
    SparseMatrix<double> star1(nEdges, nEdges);
    star1.setFromTriplets(triplets.begin(), triplets.end());
    star1.makeCompressed();

    return star1;
    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

      // 获取网格中的面（三角形）数量
    int nFaces =mesh.nFaces();

    // 创建一个稀疏矩阵，初始为零矩阵
    SparseMatrix<double> hodgeStar2(nFaces, nFaces);

    // 遍历所有面
    for (Face f : mesh.faces()) {
        // 获取面的三个顶点
        std::vector<Vertex> vertices;
        for (Vertex v : f.adjacentVertices()) {
            vertices.push_back(v);
        }

        // 确保面是三角形（应该有 3 个顶点）
        if (vertices.size() != 3) continue;

        // 获取三个顶点的位置 (A, B, C)
        Vector3 A = inputVertexPositions[vertices[0]];
        Vector3 B = inputVertexPositions[vertices[1]];
        Vector3 C = inputVertexPositions[vertices[2]];

        // 计算向量 AB 和 AC
        Vector3 AB = B - A;
        Vector3 AC = C - A;

        // 使用叉积计算三角形的面积
        double area = 0.5 * (cross(AB,AC)).norm();

        // 避免除以零；如果面积接近零则跳过
        if (area < 1e-10) continue;

        // 该面对角元素为 1/面积
        int faceIndex = f.getIndex(); // 假设 Face 提供 getIndex() 方法
        hodgeStar2.insert(faceIndex, faceIndex) = 1.0 / area;
    }

    // 完成矩阵构造
    hodgeStar2.makeCompressed();

    return hodgeStar2;

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // 获取顶点数量和边数量
    int nVertices = mesh.nVertices();
    int nEdges = mesh.nEdges();

    // 创建一个稀疏矩阵，行数为边数，列数为顶点数
    SparseMatrix<double> d0(nEdges, nVertices);

    // 遍历所有边
    for (Edge e : mesh.edges()) {
        // 获取边的索引
        int edgeIndex = e.getIndex(); // 假设 Edge 提供 getIndex() 方法

        // 获取边的起点和终点（通过 halfedge）
        Vertex v0 = e.halfedge().vertex();        // 起点
        Vertex v1 = e.halfedge().next().vertex(); // 终点

        // 获取顶点索引
        int v0Index = v0.getIndex(); // 假设 Vertex 提供 getIndex() 方法
        int v1Index = v1.getIndex(); // 假设 Vertex 提供 getIndex() 方法

        // 外微分矩阵条目：起点 +1，终点 -1
        d0.insert(edgeIndex, v0Index) = 1.0;  // 起点贡献 +1
        d0.insert(edgeIndex, v1Index) = -1.0; // 终点贡献 -1
    }

    // 完成矩阵构造
    d0.makeCompressed();
    return d0;

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // 获取面数量和边数量
    int nFaces = mesh.nFaces();
    int nEdges = mesh.nEdges();

    // 创建一个稀疏矩阵，行数为面数，列数为边数
    SparseMatrix<double> d1(nFaces, nEdges);

    // 遍历所有面
    for (Face f : mesh.faces()) {
        // 获取面的索引
        int faceIndex = f.getIndex(); // 假设 Face 提供 getIndex() 方法

        // 获取面的半边（halfedge），以遍历边界
        Halfedge he = f.halfedge();
        Halfedge heStart = he;

        // 遍历三角形的边界（三个边）
        do {
            // 获取当前半边对应的边
            Edge e = he.edge();
            int edgeIndex = e.getIndex(); // 假设 Edge 提供 getIndex() 方法

            // 确定边的方向是否与面边界方向一致
            // 如果半边的方向与边的方向一致，系数为 +1，否则为 -1
            bool isEdgeDirectionSame = (he == e.halfedge());
            double coefficient = isEdgeDirectionSame ? 1.0 : -1.0;

            // 设置矩阵条目
            d1.insert(faceIndex, edgeIndex) = coefficient;

            // 移动到下一个半边
            he = he.next();
        } while (he != heStart); // 直到回到起始半边
    }

    // 完成矩阵构造
    d1.makeCompressed();
    return d1;

    // TODO
    return identityMatrix<double>(1); // placeholder
}

} // namespace surface
} // namespace geometrycentral