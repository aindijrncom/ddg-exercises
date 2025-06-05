// Implement member functions for TreeCotree class.
#include "tree-cotree.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <set>

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
TreeCotree::TreeCotree(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build a primal spanning tree on a mesh without boundary. More specifically, populate the member variable
 * <vertexParent>, which is a std::map that maps each vertex of the input mesh to its parent in the primal spanning
 * tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildPrimalSpanningTree() {

    // TODO
    if (mesh->nVertices()==0) return; // Handle empty mesh case

    Vertex root = mesh->vertex(0); // Select the first vertex as the root
    std::queue<Vertex> queue;
    std::unordered_set<Vertex> visited;
    queue.push(root);
    visited.insert(root);

    while (!queue.empty()) {
        Vertex current = queue.front();
        queue.pop();

        // Iterate over all adjacent vertices of the current vertex
        for (Vertex neighbor : current.adjacentVertices()) {
            if (visited.find(neighbor) == visited.end()) {
                queue.push(neighbor);
                visited.insert(neighbor);
                vertexParentIndex[neighbor.getIndex()] = current.getIndex();
                vertexParent[neighbor] = current; // Set the parent of the neighbor
            }
        }
    }

}

/*
 * Check whether a halfedge is in the primal spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the primal spanning tree, false otherwise.
 */
bool TreeCotree::inPrimalSpanningTree(Halfedge he) {

    // TODO
    Vertex v1 = he.tailVertex(); // Source vertex of the halfedge
    Vertex v2 = he.tipVertex(); // Target vertex of the halfedge

    // Check if v2 is the parent of v1
    auto it1 = vertexParent.find(v1);
    if (it1 != vertexParent.end() && it1->second == v2) {
        return true;
    }

    // Check if v1 is the parent of v2
    auto it2 = vertexParent.find(v2);
    if (it2 != vertexParent.end() && it2->second == v1) {
        return true;
    }

    return false; // Halfedge is not in the spanning tree

    return false; // placeholder
}

/*
 * Build a dual spanning tree on a mesh without boundary. More specificially, populate the member variable <faceParent>,
 * which is a std::map that maps each face of the input mesh to its parent in the dual spanning tree.
 *
 * Input:
 * Returns:
/*
 * ������ż����������dual spanning cotree�����ޱ߽������ϡ�
 * ������˵������Ա���� <faceParent>������һ�� std::map�������������ÿ����ӳ�䵽���ڶ�ż���������еĸ��档
 *
 * ���룺�ޣ���������ĳ�Ա�������� mesh, vertexParent��
 * ���أ��ޣ�����洢�� faceParent �У�
 */
void TreeCotree::buildDualSpanningCoTree() {
    if (mesh->nFaces() == 0) return; // ������������

    Face root = mesh->face(0); // ѡ���һ������Ϊ��
    std::queue<Face> queue;
    std::unordered_set<Face> visited;

    queue.push(root);
    visited.insert(root);

    while (!queue.empty()) {
        Face current = queue.front();
        queue.pop();

        // ������ǰ�������������
        for (Halfedge he : current.adjacentHalfedges()) {
            Face neighbor = he.twin().face();
            if (visited.find(neighbor) == visited.end()) {
                // ������� current �� neighbor �İ���Ƿ���ԭʼ��������
                if (!inPrimalSpanningTree(he)) {
                    queue.push(neighbor);
                    visited.insert(neighbor);
                    faceParent[neighbor] = current; // �����ھӵĸ���
                    // �����Ҫ������ faceParentIndex
                    faceParentIndex[neighbor.getIndex()] = current.getIndex();
                }
            }
        }
    }
}

/*
 * Check whether a halfedge is in the dual spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the dual spanning tree, false otherwise.
 */
bool TreeCotree::inDualSpanningCotree(Halfedge he) {

    // TODO
    Face f1 = he.face();        // ��ȡ�����������
    Face f2 = he.twin().face(); // ��ȡ��ߵĶ�ż�����������

    // ��� f1 �Ƿ��� f2 �ĸ���
    auto it1 = faceParent.find(f2);
    if (it1 != faceParent.end() && it1->second == f1) {
        return true;
    }

    // ��� f2 �Ƿ��� f1 �ĸ���
    auto it2 = faceParent.find(f1);
    if (it2 != faceParent.end() && it2->second == f2) {
        return true;
    }

    return false; // ��߲��ڶ�ż��������
    return false; // placeholder
}

/*
 * Returns a halfedge lying on the shared edge between face f and g.
 *
 * Input: Two adjacent faces <f> and <g>.
 * Returns: A halfedge lying on the shared edge between face f and g.
 */
Halfedge TreeCotree::sharedHalfedge(Face f, Face g) const {

    for (Halfedge he : f.adjacentHalfedges()) {
        if (he.twin().face() == g) {
            return he;
        }
    }
    // Should never get here!
    std::cerr << "Oops, TreeCotree::sharedHalfedge() received bad input." << std::endl;
    return f.halfedge();
}

std::vector<Halfedge> TreeCotree::findPathInPrimalTree(Vertex v1, Vertex v2) {
    // ׷�� v1 ������·��
    std::vector<Vertex> path1;
    Vertex current = v1;
    while (vertexParent.find(current) != vertexParent.end()) {
        path1.push_back(current);
        current = vertexParent[current];
    }
    path1.push_back(current); // ���������

    // ׷�� v2 ������·��
    std::vector<Vertex> path2;
    current = v2;
    while (vertexParent.find(current) != vertexParent.end()) {
        path2.push_back(current);
        current = vertexParent[current];
    }
    path2.push_back(current); // ���������

    // �ҵ������ͬ���ȣ�LCA��
    Vertex lca;
    size_t i = path1.size() - 1;
    size_t j = path2.size() - 1;
    while (i >= 0 && j >= 0 && path1[i] == path2[j]) {
        lca = path1[i];
        if (i > 0) --i;
        if (j > 0) --j;
    }

    // ������ v1 �� LCA ��·��
    std::vector<Vertex> path;
    for (size_t k = 0; k <= i; ++k) {
        path.push_back(path1[k]);
    }
    path.push_back(lca);

    // ������ LCA �� v2 ��·������ת path2 �ĺ�벿�֣�
    for (size_t k = j; k > 0; --k) {
        path.push_back(path2[k]);
    }
    if (j == 0) {
        path.push_back(path2[0]);
    }

    // ������·��ת��Ϊ���·��
    std::vector<Halfedge> halfedgePath;
    for (size_t k = 0; k < path.size() - 1; ++k) {
        Vertex from = path[k];
        Vertex to = path[k + 1];
        // �����Ӷ��� from ���������а��
        for (const Halfedge& he : from.outgoingHalfedges()) {
            // ����ߵ�Ŀ�궥���Ƿ�Ϊ to
            if (he.tipVertex() == to) {
                halfedgePath.push_back(he); // �ҵ�ƥ��İ�ߣ�����
            }
        }
    }

    return halfedgePath;
}

void writeSetToFile(const std::set<int>& intSet, const std::string& filename, const std::string& delimiter = ",") {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    bool first = true;
    for (const auto& num : intSet) {
        if (!first) {
            outFile << delimiter;
        }
        outFile << num;
        first = false;
    }

    outFile.close();
}

void writeVectorToFile(const std::vector<int>& vec, const std::string& filename) {
    // ������ļ���
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }

    // ����vector��д���ļ�
    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i];

        // ����������һ��Ԫ�أ���Ӷ���
        if (i != vec.size() - 1) {
            outFile << ",";
        }
    }

    outFile.close();
    std::cout << "Successfully wrote " << vec.size() << " integers to " << filename << std::endl;
}

/*
 * ���������ͬ������Ԫ��
 *
 * ���룺�ޣ���������ĳ�Ա�������� mesh, vertexParent, faceParent��
 * ���أ��ޣ�����洢�� this->generators �У�
 */
void TreeCotree::buildGenerators() {
    // ����ԭ�����������洢�� vertexParent ��
    buildPrimalSpanningTree();
    // ������ż�����������洢�� faceParent ��
    buildDualSpanningCoTree();

    // ���������е����а�ߣ�Ѱ��ͬ������Ԫ
    std::map<int,bool> edgeIndexs; // �ؼ��ߵ�index
    for (const auto& e : mesh->halfedges()) {
        // ͬ������Ԫ��Ӧ�ڼȲ���ԭ���������С����ż��Ҳ���ڶ�ż�������еİ��
        if (!inPrimalSpanningTree(e) && !inDualSpanningCotree(e)) {
            // ��ȡ��� e ��������
            Face leftFace = e.face();
            Face rightFace = e.twin().face();

            if (edgeIndexs[e.getIndex()]) {
               continue;
            }
            edgeIndexs[e.getIndex()] = true;
            edgeIndexs[e.twin().getIndex()] = true;

            // ���������浽��ż����������·��
            std::vector<Halfedge> rightFaceCollections;
            Face currentRightFace = rightFace;
            while (faceParent.find(currentRightFace) != faceParent.end()) {
                Face parentFace = faceParent[currentRightFace];
                Halfedge sharedHe = sharedHalfedge(currentRightFace, parentFace);
                rightFaceCollections.push_back(sharedHe);
                currentRightFace = parentFace;
            }

            // ���������浽��ż����������·��
            std::vector<Halfedge> leftFaceCollections;
            Face currentLeftFace = leftFace;
            while (faceParent.find(currentLeftFace) != faceParent.end()) {
                Face parentFace = faceParent[currentLeftFace];
                Halfedge sharedHe = sharedHalfedge(currentLeftFace, parentFace);
                leftFaceCollections.push_back(sharedHe);
                currentLeftFace = parentFace;
            }

            // �����պ�ѭ��
            std::vector<Halfedge> cycle;
            cycle.push_back(e); // ��ӹؼ����

            // �ҵ������ͬ���ȣ�LCA�����ض�·��
            int m = rightFaceCollections.size() - 1;
            int n = leftFaceCollections.size() - 1;
            while (m >= 0 && n >= 0 && rightFaceCollections[m] == leftFaceCollections[n]) {
                m--;
                n--;
            }

            // �������·������ͷ�� m��
            for (int i = 0; i <= m; ++i) {
                cycle.push_back(rightFaceCollections[i]);
            }

            // �������·������ n ��ͷ������
            for (int i = n; i >= 0; --i) {
                cycle.push_back(leftFaceCollections[i].twin());
            }

            // ���պ�ѭ���洢������Ԫ�б�
            this->generators.push_back(cycle);
        }
    }

    std::set<int> indexs;
    for (auto loop : generators) {
        for (auto halfEdge : loop) {
            indexs.insert(halfEdge.tipVertex().getIndex()+1);
            indexs.insert(halfEdge.tailVertex().getIndex()+1);
        }
    }
    writeSetToFile(indexs, "data.txt", ","); // �����10;20;30
}