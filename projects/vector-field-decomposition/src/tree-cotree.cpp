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
 * 构建对偶生成余树（dual spanning cotree）在无边界网格上。
 * 具体来说，填充成员变量 <faceParent>，这是一个 std::map，将输入网格的每个面映射到其在对偶生成余树中的父面。
 *
 * 输入：无（依赖于类的成员变量，如 mesh, vertexParent）
 * 返回：无（结果存储在 faceParent 中）
 */
void TreeCotree::buildDualSpanningCoTree() {
    if (mesh->nFaces() == 0) return; // 处理空网格情况

    Face root = mesh->face(0); // 选择第一个面作为根
    std::queue<Face> queue;
    std::unordered_set<Face> visited;

    queue.push(root);
    visited.insert(root);

    while (!queue.empty()) {
        Face current = queue.front();
        queue.pop();

        // 遍历当前面的所有相邻面
        for (Halfedge he : current.adjacentHalfedges()) {
            Face neighbor = he.twin().face();
            if (visited.find(neighbor) == visited.end()) {
                // 检查连接 current 和 neighbor 的半边是否不在原始生成树中
                if (!inPrimalSpanningTree(he)) {
                    queue.push(neighbor);
                    visited.insert(neighbor);
                    faceParent[neighbor] = current; // 设置邻居的父面
                    // 如果需要，更新 faceParentIndex
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
    Face f1 = he.face();        // 获取半边所属的面
    Face f2 = he.twin().face(); // 获取半边的对偶半边所属的面

    // 检查 f1 是否是 f2 的父面
    auto it1 = faceParent.find(f2);
    if (it1 != faceParent.end() && it1->second == f1) {
        return true;
    }

    // 检查 f2 是否是 f1 的父面
    auto it2 = faceParent.find(f1);
    if (it2 != faceParent.end() && it2->second == f2) {
        return true;
    }

    return false; // 半边不在对偶生成树中
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
    // 追溯 v1 到根的路径
    std::vector<Vertex> path1;
    Vertex current = v1;
    while (vertexParent.find(current) != vertexParent.end()) {
        path1.push_back(current);
        current = vertexParent[current];
    }
    path1.push_back(current); // 加入根顶点

    // 追溯 v2 到根的路径
    std::vector<Vertex> path2;
    current = v2;
    while (vertexParent.find(current) != vertexParent.end()) {
        path2.push_back(current);
        current = vertexParent[current];
    }
    path2.push_back(current); // 加入根顶点

    // 找到最近共同祖先（LCA）
    Vertex lca;
    size_t i = path1.size() - 1;
    size_t j = path2.size() - 1;
    while (i >= 0 && j >= 0 && path1[i] == path2[j]) {
        lca = path1[i];
        if (i > 0) --i;
        if (j > 0) --j;
    }

    // 构建从 v1 到 LCA 的路径
    std::vector<Vertex> path;
    for (size_t k = 0; k <= i; ++k) {
        path.push_back(path1[k]);
    }
    path.push_back(lca);

    // 构建从 LCA 到 v2 的路径（反转 path2 的后半部分）
    for (size_t k = j; k > 0; --k) {
        path.push_back(path2[k]);
    }
    if (j == 0) {
        path.push_back(path2[0]);
    }

    // 将顶点路径转换为半边路径
    std::vector<Halfedge> halfedgePath;
    for (size_t k = 0; k < path.size() - 1; ++k) {
        Vertex from = path[k];
        Vertex to = path[k + 1];
        // 遍历从顶点 from 发出的所有半边
        for (const Halfedge& he : from.outgoingHalfedges()) {
            // 检查半边的目标顶点是否为 to
            if (he.tipVertex() == to) {
                halfedgePath.push_back(he); // 找到匹配的半边，返回
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
    // 打开输出文件流
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }

    // 遍历vector并写入文件
    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i];

        // 如果不是最后一个元素，添加逗号
        if (i != vec.size() - 1) {
            outFile << ",";
        }
    }

    outFile.close();
    std::cout << "Successfully wrote " << vec.size() << " integers to " << filename << std::endl;
}

/*
 * 计算网格的同调生成元。
 *
 * 输入：无（依赖于类的成员变量，如 mesh, vertexParent, faceParent）
 * 返回：无（结果存储在 this->generators 中）
 */
void TreeCotree::buildGenerators() {
    // 构建原生生成树，存储在 vertexParent 中
    buildPrimalSpanningTree();
    // 构建对偶生成余树，存储在 faceParent 中
    buildDualSpanningCoTree();

    // 遍历网格中的所有半边，寻找同调生成元
    std::map<int,bool> edgeIndexs; // 关键边的index
    for (const auto& e : mesh->halfedges()) {
        // 同调生成元对应于既不在原生生成树中、其对偶边也不在对偶生成树中的半边
        if (!inPrimalSpanningTree(e) && !inDualSpanningCotree(e)) {
            // 获取半边 e 的两个面
            Face leftFace = e.face();
            Face rightFace = e.twin().face();

            if (edgeIndexs[e.getIndex()]) {
               continue;
            }
            edgeIndexs[e.getIndex()] = true;
            edgeIndexs[e.twin().getIndex()] = true;

            // 构建从右面到对偶生成树根的路径
            std::vector<Halfedge> rightFaceCollections;
            Face currentRightFace = rightFace;
            while (faceParent.find(currentRightFace) != faceParent.end()) {
                Face parentFace = faceParent[currentRightFace];
                Halfedge sharedHe = sharedHalfedge(currentRightFace, parentFace);
                rightFaceCollections.push_back(sharedHe);
                currentRightFace = parentFace;
            }

            // 构建从左面到对偶生成树根的路径
            std::vector<Halfedge> leftFaceCollections;
            Face currentLeftFace = leftFace;
            while (faceParent.find(currentLeftFace) != faceParent.end()) {
                Face parentFace = faceParent[currentLeftFace];
                Halfedge sharedHe = sharedHalfedge(currentLeftFace, parentFace);
                leftFaceCollections.push_back(sharedHe);
                currentLeftFace = parentFace;
            }

            // 构建闭合循环
            std::vector<Halfedge> cycle;
            cycle.push_back(e); // 添加关键半边

            // 找到最近共同祖先（LCA）并截断路径
            int m = rightFaceCollections.size() - 1;
            int n = leftFaceCollections.size() - 1;
            while (m >= 0 && n >= 0 && rightFaceCollections[m] == leftFaceCollections[n]) {
                m--;
                n--;
            }

            // 添加右面路径（从头到 m）
            for (int i = 0; i <= m; ++i) {
                cycle.push_back(rightFaceCollections[i]);
            }

            // 添加左面路径（从 n 到头，反向）
            for (int i = n; i >= 0; --i) {
                cycle.push_back(leftFaceCollections[i].twin());
            }

            // 将闭合循环存储到生成元列表
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
    writeSetToFile(indexs, "data.txt", ","); // 输出：10;20;30
}