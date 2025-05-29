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

    // ��ȡ�����еĶ�������
    int nVertices = mesh.getVertexIndices().size();

    // ����һ��ϡ����󣬳�ʼΪ�����
    SparseMatrix<double> hodgeStar0(nVertices, nVertices);

    // ���öԽ�Ԫ��Ϊ 1������Լ����ÿ����������Ϊ 1��
    for (int i = 0; i < nVertices; ++i) {
        hodgeStar0.insert(i, i) = barycentricDualArea(mesh.vertex(i));
    }

    // ��ɾ�����
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

    // ��ȡ����ı���
    size_t nEdges = mesh.nEdges();

    // ������Ԫ���б�洢����Ԫ�أ��Խ��ߣ�
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nEdges);

    // �������б�
    for (Edge e : mesh.edges()) {
        double cotSum = 0.0;
        int faceCount = 0;

        // ��ȡ��߹������������
        Halfedge he1 = e.halfedge();
        Halfedge he2 = he1.twin();

        // �����һ�������Σ�������ڣ�
        if (!he1.face().isBoundaryLoop()) {
            // ��ȡ�������������ԵĶ���
            Vertex v_opposite = he1.next().next().vertex();

            // ����ߵ���������
            Vector3 v0 = inputVertexPositions[he1.vertex()];
            Vector3 v1 = inputVertexPositions[he1.twin().vertex()];
            Vector3 v_opp = inputVertexPositions[v_opposite];

            // �������Զ��㵽�߶˵������
            Vector3 vec1 = v0 - v_opp;
            Vector3 vec2 = v1 - v_opp;

            // ��������ֵ��cos��/sin�� = dot(v1,v2)/|v1��v2|
            double dotProd = dot(vec1, vec2);
            double crossNorm = norm(cross(vec1, vec2));
            double cot = dotProd / crossNorm;

            cotSum += cot;
            faceCount++;
        }

        // ����ڶ��������Σ�������ڣ�
        if (!he2.face().isBoundaryLoop()) {
            // ��ȡ�������������ԵĶ���
            Vertex v_opposite = he2.next().next().vertex();

            // ����ߵ���������
            Vector3 v0 = inputVertexPositions[he2.vertex()];
            Vector3 v1 = inputVertexPositions[he2.twin().vertex()];
            Vector3 v_opp = inputVertexPositions[v_opposite];

            // �������Զ��㵽�߶˵������
            Vector3 vec1 = v0 - v_opp;
            Vector3 vec2 = v1 - v_opp;

            // ��������ֵ
            double dotProd = dot(vec1, vec2);
            double crossNorm = norm(cross(vec1, vec2));
            double cot = dotProd / crossNorm;

            cotSum += cot;
            faceCount++;
        }

        // ����ƽ��ֵ�����ֻ��һ�����������Σ���ֱ��ʹ�ø�ֵ
        double value = (faceCount > 0) ? cotSum / 2.0 : 0.0;

        // ��ӶԽ���Ԫ��
        triplets.emplace_back(e.getIndex(), e.getIndex(), value);
    }

    // ���������ϡ�����
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

      // ��ȡ�����е��棨�����Σ�����
    int nFaces =mesh.nFaces();

    // ����һ��ϡ����󣬳�ʼΪ�����
    SparseMatrix<double> hodgeStar2(nFaces, nFaces);

    // ����������
    for (Face f : mesh.faces()) {
        // ��ȡ�����������
        std::vector<Vertex> vertices;
        for (Vertex v : f.adjacentVertices()) {
            vertices.push_back(v);
        }

        // ȷ�����������Σ�Ӧ���� 3 �����㣩
        if (vertices.size() != 3) continue;

        // ��ȡ���������λ�� (A, B, C)
        Vector3 A = inputVertexPositions[vertices[0]];
        Vector3 B = inputVertexPositions[vertices[1]];
        Vector3 C = inputVertexPositions[vertices[2]];

        // �������� AB �� AC
        Vector3 AB = B - A;
        Vector3 AC = C - A;

        // ʹ�ò�����������ε����
        double area = 0.5 * (cross(AB,AC)).norm();

        // ��������㣻�������ӽ���������
        if (area < 1e-10) continue;

        // ����Խ�Ԫ��Ϊ 1/���
        int faceIndex = f.getIndex(); // ���� Face �ṩ getIndex() ����
        hodgeStar2.insert(faceIndex, faceIndex) = 1.0 / area;
    }

    // ��ɾ�����
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

    // ��ȡ���������ͱ�����
    int nVertices = mesh.nVertices();
    int nEdges = mesh.nEdges();

    // ����һ��ϡ���������Ϊ����������Ϊ������
    SparseMatrix<double> d0(nEdges, nVertices);

    // �������б�
    for (Edge e : mesh.edges()) {
        // ��ȡ�ߵ�����
        int edgeIndex = e.getIndex(); // ���� Edge �ṩ getIndex() ����

        // ��ȡ�ߵ������յ㣨ͨ�� halfedge��
        Vertex v0 = e.halfedge().vertex();        // ���
        Vertex v1 = e.halfedge().next().vertex(); // �յ�

        // ��ȡ��������
        int v0Index = v0.getIndex(); // ���� Vertex �ṩ getIndex() ����
        int v1Index = v1.getIndex(); // ���� Vertex �ṩ getIndex() ����

        // ��΢�־�����Ŀ����� +1���յ� -1
        d0.insert(edgeIndex, v0Index) = 1.0;  // ��㹱�� +1
        d0.insert(edgeIndex, v1Index) = -1.0; // �յ㹱�� -1
    }

    // ��ɾ�����
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

    // ��ȡ�������ͱ�����
    int nFaces = mesh.nFaces();
    int nEdges = mesh.nEdges();

    // ����һ��ϡ���������Ϊ����������Ϊ����
    SparseMatrix<double> d1(nFaces, nEdges);

    // ����������
    for (Face f : mesh.faces()) {
        // ��ȡ�������
        int faceIndex = f.getIndex(); // ���� Face �ṩ getIndex() ����

        // ��ȡ��İ�ߣ�halfedge�����Ա����߽�
        Halfedge he = f.halfedge();
        Halfedge heStart = he;

        // ���������εı߽磨�����ߣ�
        do {
            // ��ȡ��ǰ��߶�Ӧ�ı�
            Edge e = he.edge();
            int edgeIndex = e.getIndex(); // ���� Edge �ṩ getIndex() ����

            // ȷ���ߵķ����Ƿ�����߽緽��һ��
            // �����ߵķ�����ߵķ���һ�£�ϵ��Ϊ +1������Ϊ -1
            bool isEdgeDirectionSame = (he == e.halfedge());
            double coefficient = isEdgeDirectionSame ? 1.0 : -1.0;

            // ���þ�����Ŀ
            d1.insert(faceIndex, edgeIndex) = coefficient;

            // �ƶ�����һ�����
            he = he.next();
        } while (he != heStart); // ֱ���ص���ʼ���
    }

    // ��ɾ�����
    d1.makeCompressed();
    return d1;

    // TODO
    return identityMatrix<double>(1); // placeholder
}

} // namespace surface
} // namespace geometrycentral