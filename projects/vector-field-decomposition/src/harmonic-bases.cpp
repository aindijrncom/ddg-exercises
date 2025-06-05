// Implement member functions for HarmonicBases class.
#include "harmonic-bases.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HarmonicBases::HarmonicBases(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build a closed, but not exact, primal 1-form ω.
 *
 * Input: A std::vector of Halfedges representing a homology generator of the mesh.
 * Returns: A vector representing a closed primal 1-form.
 */
Vector<double> HarmonicBases::buildClosedPrimalOneForm(const std::vector<Halfedge>& generator) const {

    // TODO
    // Initialize a vector of zeros with size equal to the number of halfedges in the mesh
    Vector<double> omega = Vector<double>::Zero(mesh->nEdges());

    // For each halfedge in the generator, assign 1 to it and -1 to its twin
    for (const auto& h : generator) {
        double value = 1.0f;
        if (h.edge().halfedge() != h) {
            value = -value;
        }
        omega[h.edge().getIndex()] = value; // Assign 1 to halfedges in the generator
    }

    return omega;
    return Vector<double>::Zero(1); // placeholder
}

/*
 * Compute the harmonic bases [γ1, γ2 ... γn] of the input mesh.
 *
 * Input: A std::vector of homology generators of the mesh (which are in turn represented as std::vectors of halfedges),
 * and a HodgeDecomposition object. Returns:
 */
std::vector<Vector<double>> HarmonicBases::compute(const std::vector<std::vector<Halfedge>>& generators,
                                                   const HodgeDecomposition& hodgeDecomposition) const {
    // Initialize an empty vector to store the harmonic bases
    std::vector<Vector<double>> gammas;

    // Process each homology generator
    for (const auto& generator : generators) {
        // Build a closed primal 1-form for the current generator
        Vector<double> omega = buildClosedPrimalOneForm(generator);

        // Compute its harmonic part using the HodgeDecomposition object
        Vector<double> harmonic_part = hodgeDecomposition.computeHarmonicComponent(
            omega, hodgeDecomposition.computeExactComponent(omega), hodgeDecomposition.computeCoExactComponent(omega));

        // Add the harmonic part to the list of harmonic bases
        gammas.push_back(harmonic_part);
    }

    return gammas;

    // TODO
    return gammas; // placeholder
}