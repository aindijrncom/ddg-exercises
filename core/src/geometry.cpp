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
        he = he.next().twin(); // 移动到下一个半边
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

    // 初始化角和与混合面积
    double angleSum = 0.0;  // 顶点周围的内角和
    double mixedArea = 0.0; // 混合面积

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
        double faceNormalMag = sqrt(dot(faceNormal, faceNormal));
        if (faceNormalMag < 1e-10) continue; // 避免退化三角形
        faceNormal = {faceNormal.x / faceNormalMag, faceNormal.y / faceNormalMag, faceNormal.z / faceNormalMag};

        // 计算三角形面积
        double area = 0.5 * faceNormalMag;
        mixedArea += area / 3.0; // 重心对偶面积贡献

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

        // 加权法向量（暂时使用面积权重，后续结合高斯曲率调整）
        zeronormal = zeronormal + faceNormal * (area / 3.0);
    }

    // 计算离散高斯曲率
    double K = (2.0 * M_PI - angleSum) / mixedArea;
    if (mixedArea < 1e-10) return zeronormal.zero(); // 避免除以零

    // 修正法向量：高斯曲率影响权重（这里简化，直接使用面积加权结果）
    // 归一化最终法向量
    double normalMag = sqrt(dot(zeronormal, zeronormal));
    if (normalMag < 1e-10) return zeronormal.zero(); // 避免除以零
    zeronormal = {zeronormal.x / normalMag, zeronormal.y / normalMag, zeronormal.z / normalMag};
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
    Vector3 zeronormal;
    zeronormal = zeronormal.zero();

    // 计算顶点 v 的混合面积
    double dualArea = circumcentricDualArea(v);
    if (dualArea < 1e-10) return zeronormal.zero(); // 避免除以零

    // 遍历顶点 v 周围的半边
    Halfedge startHe = v.halfedge();
    Halfedge he = startHe;
    do {
        // 获取当前半边对应的边和顶点
        Vertex vStart = he.vertex();
        Vertex vEnd = he.next().vertex();
        Vector3 xStart = inputVertexPositions[vStart];
        Vector3 xEnd = inputVertexPositions[vEnd];

        // 计算边的向量
        Vector3 edgeVec = xEnd - xStart;

        // 获取相邻面
        Face f1 = he.face();
        Face f2 = he.twin().face();

        // 计算与边相关的内角的余切值
        double cotAlpha = 0.0, cotBeta = 0.0;
        if (1) {
            std::vector<Vertex> verticesF1;
            for (Vertex vert : f1.adjacentVertices()) verticesF1.push_back(vert);
            int v1Index = -1;
            for (int i = 0; i < 3; ++i)
                if (verticesF1[i] == v) v1Index = i;
            Vector3 vPrevF1 =
                (v1Index == 0) ? inputVertexPositions[verticesF1[2]] : inputVertexPositions[verticesF1[v1Index - 1]];
            Vector3 vNextF1 =
                (v1Index == 2) ? inputVertexPositions[verticesF1[0]] : inputVertexPositions[verticesF1[v1Index + 1]];
            Vector3 edgePrevF1 = vPrevF1 - inputVertexPositions[v];
            Vector3 edgeNextF1 = vNextF1 - inputVertexPositions[v];
            double lenPrevF1 = sqrt(dot(edgePrevF1, edgePrevF1));
            double lenNextF1 = sqrt(dot(edgeNextF1, edgeNextF1));
            edgePrevF1 = {edgePrevF1.x / lenPrevF1, edgePrevF1.y / lenPrevF1, edgePrevF1.z / lenPrevF1};
            edgeNextF1 = {edgeNextF1.x / lenNextF1, edgeNextF1.y / lenNextF1, edgeNextF1.z / lenNextF1};
            double cosAlpha = dot(edgePrevF1, edgeNextF1);
            cosAlpha = std::max(-1.0, std::min(1.0, cosAlpha));
            double sinAlpha = sqrt(1.0 - cosAlpha * cosAlpha);
            cotAlpha = (sinAlpha < 1e-10) ? 0.0 : cosAlpha / sinAlpha;
        }
        if (1) {
            std::vector<Vertex> verticesF2;
            for (Vertex vert : f2.adjacentVertices()) verticesF2.push_back(vert);
            int v2Index = -1;
            for (int i = 0; i < 3; ++i)
                if (verticesF2[i] == v) v2Index = i;
            Vector3 vPrevF2 =
                (v2Index == 0) ? inputVertexPositions[verticesF2[2]] : inputVertexPositions[verticesF2[v2Index - 1]];
            Vector3 vNextF2 =
                (v2Index == 2) ? inputVertexPositions[verticesF2[0]] : inputVertexPositions[verticesF2[v2Index + 1]];
            Vector3 edgePrevF2 = vPrevF2 - inputVertexPositions[v];
            Vector3 edgeNextF2 = vNextF2 - inputVertexPositions[v];
            double lenPrevF2 = sqrt(dot(edgePrevF2, edgePrevF2));
            double lenNextF2 = sqrt(dot(edgeNextF2, edgeNextF2));
            edgePrevF2 = {edgePrevF2.x / lenPrevF2, edgePrevF2.y / lenPrevF2, edgePrevF2.z / lenPrevF2};
            edgeNextF2 = {edgeNextF2.x / lenNextF2, edgeNextF2.y / lenNextF2, edgeNextF2.z / lenNextF2};
            double cosBeta = dot(edgePrevF2, edgeNextF2);
            cosBeta = std::max(-1.0, std::min(1.0, cosBeta));
            double sinBeta = sqrt(1.0 - cosBeta * cosBeta);
            cotBeta = (sinBeta < 1e-10) ? 0.0 : cosBeta / sinBeta;
        }

        // 计算边贡献
        double weight = (cotAlpha + cotBeta) / 2.0;
        zeronormal = zeronormal + weight * edgeVec;

        // 移动到下一个半边
        he = he.twin().next();
    } while (he != startHe);

    // 归一化法向量
    double normalMag = sqrt(dot(zeronormal, zeronormal));
    if (normalMag < 1e-10) return zeronormal.zero(); // 避免除以零
    zeronormal = {zeronormal.x / normalMag, zeronormal.y / normalMag, zeronormal.z / normalMag};

    // 应用混合面积归一化
    zeronormal = zeronormal * (0.5 / dualArea);

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
    // 初始化对偶面积
    double area = 0.0;

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

        // 确定顶点 v 和相邻顶点
        Vector3 vPos = inputVertexPositions[v];
        Vector3 vPrev = (vIndex == 0) ? inputVertexPositions[vertices[2]] : inputVertexPositions[vertices[vIndex - 1]];
        Vector3 vNext = (vIndex == 2) ? inputVertexPositions[vertices[0]] : inputVertexPositions[vertices[vIndex + 1]];

        // 计算边向量和长度
        Vector3 edgePrev = vPrev - vPos;
        Vector3 edgeNext = vNext - vPos;
        double lenPrev = sqrt(dot(edgePrev, edgePrev));
        double lenNext = sqrt(dot(edgeNext, edgeNext));
        if (lenPrev < 1e-10 || lenNext < 1e-10) continue;

        // 归一化边向量
        Vector3 edgePrevNorm = {edgePrev.x / lenPrev, edgePrev.y / lenPrev, edgePrev.z / lenPrev};
        Vector3 edgeNextNorm = {edgeNext.x / lenNext, edgeNext.y / lenNext, edgeNext.z / lenNext};

        // 计算内角的余切值
        double cosTheta = dot(edgePrevNorm, edgeNextNorm);
        cosTheta = std::max(-1.0, std::min(1.0, cosTheta)); // 限制范围
        double sinTheta = sqrt(1.0 - cosTheta * cosTheta);
        double cotTheta = (sinTheta < 1e-10) ? 0.0 : cosTheta / sinTheta;

        // 计算对偶面积贡献：(l1^2 + l2^2) * cot θ / 8
        double contribution = (lenPrev * lenPrev + lenNext * lenNext) * cotTheta / 8.0;
        area += contribution;
    }

    // 确保面积非负
    if (area < 0) area = 0.0;
    return area;
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
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
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