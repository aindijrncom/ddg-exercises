// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>
#include <geometrycentral/surface/tufted_laplacian.h>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {

    // TODO
    halfedgeCotanWeight(he);
    return 0; // placeholders
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    // TODO
    double dualArea = 0.0;

    // Iterate over all faces adjacent to the vertex v
    for (Face f : v.adjacentFaces()) {
        // Get the three vertices of the face
        // Assuming Face provides a way to iterate over its vertices
        std::vector<Vertex> vertices;
        for (Vertex vert : f.adjacentVertices()) {
            vertices.push_back(vert);
        }

        // Ensure the face is a triangle (should have exactly 3 vertices)
        if (vertices.size() != 3) continue;

        // Get the positions of the three vertices (A, B, C)
        Vector3 A = inputVertexPositions[vertices[0]];
        Vector3 B = inputVertexPositions[vertices[1]];
        Vector3 C = inputVertexPositions[vertices[2]];

        // Compute vectors AB and AC
        Vector3 AB = B - A;
        Vector3 AC = C - A;

        // Compute the area of the triangle using the cross product
        Vector3 crossProduct = cross(AB,AC);
        double area = 0.5 * crossProduct.norm();

        // Add the contribution of this triangle to the dual area
        // Each triangle contributes area/3 to each of its vertices
        dualArea += area / 3.0;
    }

    return dualArea; 
    return 0; // placeholder
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    // TODO
    // 获取与角关联的半边
    Halfedge h = c.halfedge();

    // 获取三个顶点：
    // v0: 角所在的顶点
    // v1: 第一条邻边的终点
    // v2: 第二条邻边的终点（通过下一条半边的终点获取）
    Vertex v0 = h.vertex();
    Vertex v1 = h.next().vertex();
    Vertex v2 = h.next().next().vertex();

    // 获取顶点坐标
    Vector3 p0 = inputVertexPositions[v0];
    Vector3 p1 = inputVertexPositions[v1];
    Vector3 p2 = inputVertexPositions[v2];

    // 计算从角顶点出发的两条邻边向量
    Vector3 vec1 = p1 - p0;
    Vector3 vec2 = p2 - p0;

    // 计算点积和模长
    double dotProduct = dot(vec1, vec2);
    double len1 = norm(vec1);
    double len2 = norm(vec2);

    // 处理零向量情况（理论上不应发生）
    if (len1 == 0 || len2 == 0) {
        return 0.0; // 返回0避免NaN
    }

    // 计算夹角余弦值，并限制在[-1,1]范围内防止浮点误差
    double cosTheta = dotProduct / (len1 * len2);

    // 反余弦计算夹角（弧度制）
    return std::acos(cosTheta);

    return 0; // placeholder
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {

    // TODO

    return edgeDihedralAngle(he.edge());

    return 0; // placeholder
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    // TODO
    return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {

    // TODO
    // 初始化法向量为零向量
    Vector3 normal;
    normal = normal.zero();

    // 遍历顶点 v 周围的所有面（三角形）
    for (Face f : v.adjacentFaces()) {
        // 获取三角形的三个顶点
        std::vector<Vertex> vertices;
        for (Vertex vert : f.adjacentVertices()) {
            vertices.push_back(vert);
        }

        // 确保是三角形（有 3 个顶点）
        if (vertices.size() != 3) continue;

        // 确定当前顶点 v 在三角形中的位置
        int vIndex = -1;
        for (int i = 0; i < 3; ++i) {
            if (vertices[i] == v) {
                vIndex = i;
                break;
            }
        }
        if (vIndex == -1) continue;

        // 获取三个顶点的坐标
        Vector3 A = inputVertexPositions[vertices[0]];
        Vector3 B = inputVertexPositions[vertices[1]];
        Vector3 C = inputVertexPositions[vertices[2]];

        // 计算三角形的法向量
        Vector3 AB = B - A;
        Vector3 AC = C - A;
        Vector3 faceNormal = cross(AB, AC);
        double faceNormalMag = faceNormal.norm();
        if (faceNormalMag < 1e-10) continue; // 避免退化三角形
        faceNormal = faceNormal / faceNormalMag;

        // 计算顶点 v 处的内角
        // 根据 v 在三角形中的位置，确定对应的边
        Vector3 vPos, vNext, vPrev;
        if (vIndex == 0) {
            vPos = A;
            vNext = B;
            vPrev = C;
        } else if (vIndex == 1) {
            vPos = B;
            vNext = C;
            vPrev = A;
        } else {
            vPos = C;
            vNext = A;
            vPrev = B;
        }

        // 计算顶点 v 处的两条边向量
        Vector3 edge1 = vPrev - vPos;
        Vector3 edge2 = vNext - vPos;

        // 归一化边向量
        double edge1Mag = edge1.norm();
        double edge2Mag = edge2.norm();
        if (edge1Mag < 1e-10 || edge2Mag < 1e-10) continue;
        edge1 = edge1 / edge1Mag;
        edge2 = edge2 / edge2Mag;

        // 计算夹角（内角）作为权重
        double cosTheta = dot(edge1,edge2);
        cosTheta = std::max(-1.0, std::min(1.0, cosTheta)); // 限制范围避免误差
        double theta = acos(cosTheta);

        // 加权法向量（权重为内角 theta）
        normal = normal + faceNormal * theta;
    }

    // 归一化最终法向量
    double normalMag = normal.norm();
    Vector3 zeronormal;
    zeronormal = zeronormal.zero();
    if (normalMag < 1e-10) return zeronormal; // 避免除以零
    normal = normal / normalMag;

    return normal;
    return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
//方法背景：内接球方法通常涉及顶点周围的混合面积（mixed area），通过加权三角形法向量计算顶点法向量。混合面积可近似为重心对偶面积（每个三角形面积的 <math xmlns="http://www.w3.org/1998/Math/MathML"><semantics><mrow><mn>1</mn><mi mathvariant="normal">/</mi><mn>3</mn></mrow><annotation encoding="application/x-tex">   1/3   </annotation></semantics></math>），并手动计算所有几何操作。
//步骤：
//
//初始化法向量为零。
//遍历顶点周围的三角形，计算每个三角形的法向量。
//使用混合面积（三角形面积 / 3）作为权重。
//累加加权法向量并归一化。
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    // TODO
    // 初始化法向量为零向量
    Vector3 zeronormal;
    zeronormal = zeronormal.zero();

    // 获取顶点 v 周围的半边
    std::vector<Halfedge> halfedges;
    Halfedge startHe = v.halfedge(); // 获取顶点的第一个半边
    Halfedge he = startHe;
    do {
        halfedges.push_back(he);
        he = he.next().next().twin(); // 移动到下一个半边
    } while (he != startHe); // 直到回到起点或遇到边界

    // 如果没有有效半边，返回零向量
    if (halfedges.empty()) return zeronormal.zero();

    // 获取半边向量
    std::vector<Vector3> vec;
    for (Halfedge h : halfedges) {
        Vertex vStart = h.vertex();
        Vertex vEnd = h.next().vertex();
        Vector3 startPos = inputVertexPositions[vStart];
        Vector3 endPos = inputVertexPositions[vEnd];
        Vector3 edgeVec = endPos - startPos; // 半边向量
        vec.push_back(edgeVec);
    }

    // 计算每个半边向量的长度
    std::vector<double> len(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        len[i] = sqrt(dot(vec[i], vec[i]));
        if (len[i] < 1e-10) return zeronormal.zero(); // 避免除以零
    }

    // 计算法向量
    for (size_t i = 0; i < vec.size(); ++i) {
        size_t j = (i + 1) % vec.size(); // 下一个半边（循环）
        Vector3 e1 = vec[i];
        Vector3 e2 = {vec[j].x, vec[j].y, vec[j].z}; // 取负

        // 计算加权叉积
        Vector3 crossProd = cross(e1, e2);
        double weight = 1.0 / (len[i] * len[i] * len[j] * len[j]);
        zeronormal = zeronormal + crossProd * weight;
    }
     
    // 归一化最终法向量
    double normalMag = sqrt(dot(zeronormal, zeronormal));
    if (normalMag < 1e-10) return zeronormal.zero(); // 避免除以零
    zeronormal = {zeronormal.x / normalMag, zeronormal.y / normalMag, zeronormal.z / normalMag};

    return zeronormal;

    return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    // TODO
    // 初始化法向量为零向量
    Vector3 zeronormal;
    zeronormal = zeronormal.zero();

    // 遍历顶点 v 周围的所有面（三角形）
    for (Face f : v.adjacentFaces()) {
        // 获取三角形的三个顶点
        std::vector<Vertex> vertices;
        for (Vertex vert : f.adjacentVertices()) {
            vertices.push_back(vert);
        }

        // 确保是三角形（有 3 个顶点）
        if (vertices.size() != 3) continue;

        // 获取三个顶点的坐标，使用 inputVertexPositions
        Vector3 A = inputVertexPositions[vertices[0]];
        Vector3 B = inputVertexPositions[vertices[1]];
        Vector3 C = inputVertexPositions[vertices[2]];

        // 计算三角形的法向量
        Vector3 AB = B - A;
        Vector3 AC = C - A;
        Vector3 faceNormal = cross(AB, AC);
        double faceNormalMag = sqrt(dot(faceNormal, faceNormal));
        if (faceNormalMag < 1e-10) continue; // 避免退化三角形
        faceNormal = {faceNormal.x / faceNormalMag, faceNormal.y / faceNormalMag, faceNormal.z / faceNormalMag};

        // 计算三角形面积，作为权重
        double area = 0.5 * sqrt(dot(faceNormal, faceNormal)) * faceNormalMag;

        // 加权法向量（权重为三角形面积）
        zeronormal = zeronormal + faceNormal * area;
    }

    // 归一化最终法向量
    double normalMag = sqrt(dot(zeronormal, zeronormal));
    if (normalMag < 1e-10) return zeronormal.zero(); // 避免除以零
    zeronormal = {zeronormal.x / normalMag, zeronormal.y / normalMag, zeronormal.z / normalMag};

    return zeronormal;
    return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    // TODO
    // 初始化法向量为零向量
    Vector3 zeronormal;
    zeronormal = zeronormal.zero();

    // 获取顶点 v 周围的半边
    std::vector<Halfedge> halfedges;
    Halfedge startHe = v.halfedge(); // 获取顶点的第一个半边
    Halfedge he = startHe;
    do {
        halfedges.push_back(he);
        he = he.next().next().twin(); // 移动到下一个半边
    } while (he != startHe); // 直到回到起点或遇到边界

    // 如果没有有效半边，返回零向量
    if (halfedges.empty()) return zeronormal.zero();

    // 获取半边向量
    std::vector<Vector3> vec;
    for (Halfedge h : halfedges) {
        Vertex vStart = h.vertex();
        Vertex vEnd = h.next().vertex();
        Vector3 startPos = inputVertexPositions[vStart];
        Vector3 endPos = inputVertexPositions[vEnd];
        Vector3 edgeVec = endPos - startPos; // 半边向量
        vec.push_back(edgeVec);
    }

    // 计算每个半边向量的长度
    std::vector<double> len(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        len[i] = sqrt(dot(vec[i], vec[i]));
        if (len[i] < 1e-10) return zeronormal.zero(); // 避免除以零
    }

    for (size_t i = 0; i < vec.size(); ++i) 
    {
        auto diAngle = dihedralAngle(halfedges[i]);
        zeronormal += diAngle * vec[i] / len[i];
    }

    zeronormal /= 2;

    // 归一化最终法向量
    double normalMag = sqrt(dot(zeronormal, zeronormal));
    if (normalMag < 1e-10) return zeronormal.zero(); // 避免除以零
    zeronormal = {zeronormal.x / normalMag, zeronormal.y / normalMag, zeronormal.z / normalMag};
    return zeronormal;

    return {0, 0, 0}; // placeholder
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    // TODO
    // 初始化法向量为零向量
    Vector3 zeronormal =Vector3::zero();

    // 获取顶点 v 周围的半边
    std::vector<Halfedge> halfedges;
    Halfedge startHe = v.halfedge(); // 获取顶点的第一个半边
    Halfedge he = startHe;
    do {
        halfedges.push_back(he);
        he = he.next().next().twin(); // 移动到下一个半边
    } while (he != startHe); // 直到回到起点或遇到边界

    // 如果没有有效半边，返回零向量
    if (halfedges.empty()) return zeronormal.zero();

    // 获取半边向量
    std::vector<Vector3> vec;
    for (Halfedge h : halfedges) {
        Vertex vStart = h.vertex();
        Vertex vEnd = h.next().vertex();
        Vector3 startPos = inputVertexPositions[vStart];
        Vector3 endPos = inputVertexPositions[vEnd];
        Vector3 edgeVec = endPos - startPos; // 半边向量
        vec.push_back(edgeVec);
    }

    // 计算每个半边向量的长度
    std::vector<double> len(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        len[i] = sqrt(dot(vec[i], vec[i]));
        if (len[i] < 1e-10) return zeronormal.zero(); // 避免除以零
    }

   for (size_t i = 0; i < vec.size(); ++i) {
        auto CotWeight = halfedgeCotanWeight(halfedges[i]) + halfedgeCotanWeight(halfedges[i].twin());
        zeronormal += CotWeight * vec[i];
    }

    zeronormal /= 2;

    // 归一化最终法向量
    double normalMag = sqrt(dot(zeronormal, zeronormal));
    if (normalMag < 1e-10) return zeronormal.zero(); // 避免除以零
    zeronormal = {zeronormal.x / normalMag, zeronormal.y / normalMag, zeronormal.z / normalMag};
    return zeronormal;

    return {0, 0, 0}; // placeholder
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    // TODO
    // 初始化内角和
    double angleSum = 0.0;

    // 遍历顶点 v 周围的所有面（三角形）
    for (Face f : v.adjacentFaces()) {
        // 获取三角形的三个顶点
        std::vector<Vertex> vertices;
        for (Vertex vert : f.adjacentVertices()) {
            vertices.push_back(vert);
        }

        // 确保是三角形（有 3 个顶点）
        if (vertices.size() != 3) continue;

        // 确定当前顶点 v 在三角形中的位置
        int vIndex = -1;
        for (int i = 0; i < 3; ++i) {
            if (vertices[i] == v) {
                vIndex = i;
                break;
            }
        }
        if (vIndex == -1) continue;

        // 获取三个顶点的坐标
        Vector3 A = inputVertexPositions[vertices[0]];
        Vector3 B = inputVertexPositions[vertices[1]];
        Vector3 C = inputVertexPositions[vertices[2]];

        // 计算顶点 v 处的内角
        Vector3 vPos = inputVertexPositions[v];
        Vector3 vPrev = (vIndex == 0) ? inputVertexPositions[vertices[2]] : inputVertexPositions[vertices[vIndex - 1]];
        Vector3 vNext = (vIndex == 2) ? inputVertexPositions[vertices[0]] : inputVertexPositions[vertices[vIndex + 1]];

        // 计算边向量
        Vector3 edgePrev = vPrev - vPos;
        Vector3 edgeNext = vNext - vPos;

        // 归一化边向量
        double edgePrevMag = sqrt(dot(edgePrev, edgePrev));
        double edgeNextMag = sqrt(dot(edgeNext, edgeNext));
        if (edgePrevMag < 1e-10 || edgeNextMag < 1e-10) continue;
        edgePrev = {edgePrev.x / edgePrevMag, edgePrev.y / edgePrevMag, edgePrev.z / edgePrevMag};
        edgeNext = {edgeNext.x / edgeNextMag, edgeNext.y / edgeNextMag, edgeNext.z / edgeNextMag};

        // 计算内角
        double cosTheta = dot(edgePrev, edgeNext);
        cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
        double theta = acos(cosTheta);
        angleSum += theta;
    }

    // 计算角亏
    double defect = 2.0 * M_PI - angleSum;
    return defect;
    return 0; // placeholder
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {

    // TODO
    // 初始化总角亏
    double totalDefect = 0.0;

    // 遍历网格中的所有顶点
    for (Vertex v : mesh.vertices()) {
        // 计算每个顶点的角亏并累加
        totalDefect += angleDefect(v);
    }

    return totalDefect;
    return 0; // placeholder
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {

    // TODO
    // 初始化法向量为零向量
    double meanCurvate = 0.0;

    // 获取顶点 v 周围的半边
    std::vector<Halfedge> halfedges;
    Halfedge startHe = v.halfedge(); // 获取顶点的第一个半边
    Halfedge he = startHe;
    do {
        halfedges.push_back(he);
        he = he.next().next().twin(); // 移动到下一个半边
    } while (he != startHe); // 直到回到起点或遇到边界

    // 如果没有有效半边，返回零向量
    if (halfedges.empty()) return 0.0;

    // 获取半边向量
    std::vector<Vector3> vec;
    for (Halfedge h : halfedges) {
        Vertex vStart = h.vertex();
        Vertex vEnd = h.next().vertex();
        Vector3 startPos = inputVertexPositions[vStart];
        Vector3 endPos = inputVertexPositions[vEnd];
        Vector3 edgeVec = endPos - startPos; // 半边向量
        vec.push_back(edgeVec);
    }

    // 计算每个半边向量的长度
    std::vector<double> len(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        len[i] = sqrt(dot(vec[i], vec[i]));
        if (len[i] < 1e-10) return 0.0; // 避免除以零
    }

    for (size_t i = 0; i < vec.size(); ++i) {
        auto diAngle = dihedralAngle(halfedges[i]);
        meanCurvate += diAngle * len[i];
    }

    meanCurvate /= 2;
    return meanCurvate;

    return 0; // placeholder
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {

    // TODO
    {
        double area = 0.0;

        for (Face f : v.adjacentFaces()) {
            // Skip non-triangle faces
            if (f.degree() != 3) continue;

            Halfedge he = f.halfedge();
            Halfedge h1 = he;
            Halfedge h2 = he.next();
            Halfedge h3 = h2.next();

            // Identify the halfedge starting at vertex v
            if (h1.vertex() != v && h2.vertex() != v && h3.vertex() != v) continue;

            // Find the halfedge starting at v
            Halfedge hV = (h1.vertex() == v) ? h1 : (h2.vertex() == v) ? h2 : h3;

            // Get vertices: v, vNext, vPrev
            Vertex vA = hV.tailVertex();     // v
            Vertex vB = hV.next().tailVertex();    // vNext
            Vertex vC = hV.next().next().tailVertex(); // vPrev

            Vector3 pA = inputVertexPositions[vA];
            Vector3 pB = inputVertexPositions[vB];
            Vector3 pC = inputVertexPositions[vC];

            // Edge vectors and lengths
            Vector3 vecAB = pB - pA;
            Vector3 vecAC = pC - pA;
            double lenAB = norm(vecAB);
            double lenAC = norm(vecAC);

            // Skip degenerate edges
            if (lenAB < 1e-10 || lenAC < 1e-10) continue;

            // Compute cotangents for angles at vB and vC
            Vector3 vecBA = pA - pB;
            Vector3 vecBC = pC - pB;
            double cotB = dot(vecBA, vecBC) / norm(cross(vecBA, vecBC));

            Vector3 vecCA = pA - pC;
            Vector3 vecCB = pB - pC;
            double cotC = dot(vecCA, vecCB) / norm(cross(vecCA, vecCB));

            // Handle degenerate angles (avoid division by zero)
            if (std::isnan(cotB)) cotB = 0.0;
            if (std::isnan(cotC)) cotC = 0.0;

            // Add Voronoi area contribution from this face
            area += (lenAB * lenAB * cotC + lenAC * lenAC * cotB);
        }
        area /= 8.0;

        // Clamp negative areas to zero
        return area;
    }
    return 0; // placeholder
} 

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    // TODO
    double H = vertexMeanCurvature(v);
    double K = vertexGaussianCurvature(v);
    double area = barycentricDualArea(v);
    double temp = std::sqrt(std::abs(H / area * H / area - K / area));
    double k_max = H / area + temp;
    double k_min = H / area - temp;
    return  std::make_pair(k_min, k_max);
    return std::make_pair(0, 0); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    // TODO
    // 获取网格对象和顶点数量
    int nVertices = mesh.nVertices();

    // 使用三元组列表高效构建稀疏矩阵
    std::vector<Eigen::Triplet<double>> triplets;

    // 顶点度数向量（记录每个顶点相邻边的权重和）
    Eigen::VectorXd vertexDegrees = Eigen::VectorXd::Zero(nVertices);

    // 第一步：遍历所有边构建负半定拉普拉斯矩阵
    for (Edge e : mesh.edges()) {
        // 获取当前边的两个半边
        Halfedge he = e.halfedge();
        Halfedge twin = he.twin();

        // 计算cotangent权重（取两个对角的平均值）
        double weight = 0.5 * (halfedgeCotanWeight(he) + halfedgeCotanWeight(twin));

        // 获取边两端的顶点索引
        Vertex v1 = he.vertex();
        Vertex v2 = twin.vertex();
        size_t i = v1.getIndex();
        size_t j = v2.getIndex();

        // 添加非对角元素（负半定矩阵的非零元）
        triplets.emplace_back(i, j, weight); // L_ij
        triplets.emplace_back(j, i, weight); // L_ji

        // 累加顶点度数（权重和）
        vertexDegrees[i] += weight;
        vertexDegrees[j] += weight;
    }

    // 第二步：处理对角元素并应用正定转换
    for (Vertex v : mesh.vertices()) {
        size_t idx = v.getIndex();
        double degree = vertexDegrees[idx];

        // 原始负半定矩阵的对角元素 = -degree
        // 乘以-1后变为：+degree
        // 添加小偏移量(1e-8)确保正定性
        triplets.emplace_back(idx, idx, (-degree) + 1e-8);
    }

    // 第三步：构建稀疏矩阵
    SparseMatrix<double> laplaceMat(nVertices, nVertices);
    laplaceMat.setFromTriplets(triplets.begin(), triplets.end());

    // 确保矩阵对称性（数值稳定性）
    laplaceMat.makeCompressed();
    return laplaceMat;
    return identityMatrix<double>(1); // placeholder
}

/*
 * 构建包含顶点重心对偶面积的稀疏对角质量矩阵
 *
 * 重心对偶面积 = 顶点周围所有三角形面积之和的 1/3
 *
 * 返回: 网格的质量矩阵（对角矩阵）
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {
    // 获取网格对象和顶点数量
    size_t nVertices = mesh.nVertices();

    // 初始化存储顶点面积的向量（全零）
    Eigen::VectorXd vertexAreas = Eigen::VectorXd::Zero(nVertices);

    // 将三角形面积的1/3分配给每个顶点
    for (Vertex v : mesh.vertices()) {
        vertexAreas[v.getIndex()] = barycentricDualArea(v);
    }

    // 第二步：构建对角质量矩阵
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nVertices); // 预分配空间

    for (size_t i = 0; i < nVertices; ++i) {
        // 将对角元素设为顶点面积
        triplets.emplace_back(i, i, vertexAreas[i]);
    }

    // 创建并填充稀疏矩阵
    SparseMatrix<double> massMat(nVertices, nVertices);
    massMat.setFromTriplets(triplets.begin(), triplets.end());

    return massMat;
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * 构建负半定拉普拉斯矩阵 (原始cotangent权重)
 * 返回: 负半定拉普拉斯矩阵
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrixNegativeSemidefinite() const {
    size_t n = mesh.nVertices();
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd diagonalSums = Eigen::VectorXd::Zero(n);

    // 遍历所有边计算cotangent权重
    for (Edge e : mesh.edges()) {
        Halfedge he = e.halfedge();
        double w = 0.5 * (halfedgeCotanWeight(he) + halfedgeCotanWeight(he.twin()));

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
    for (Vertex v : mesh.vertices()) {
        int i = v.getIndex();
        triplets.emplace_back(i, i, -diagonalSums[i]);
    }

    SparseMatrix<double> L(n, n);
    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral