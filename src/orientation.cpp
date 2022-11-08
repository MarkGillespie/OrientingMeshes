#include "orientation.h"

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

std::pair<std::unique_ptr<ManifoldSurfaceMesh>, OrientationCoverMapping>
constructOrientationCover(SurfaceMesh& mesh, const VertexData<Vector3>& pos) {

    size_t nH = mesh.nInteriorHalfedges();

    // Index halfedges
    // halfedge iH in the original mesh maps to halfedge iH (with the same
    // orientation) and halfedge iH+nH (with the opposite orientation) in the
    // orientation cover
    HalfedgeData<size_t> hIdx(mesh), hIdxOpp(mesh);
    size_t iH = 0;
    for (Halfedge he : mesh.interiorHalfedges()) {
        hIdx[he]    = iH;
        hIdxOpp[he] = iH + nH;
        iH++;
    }

    size_t nCoverH = 2 * nH;

    std::vector<size_t> next = std::vector<size_t>(nCoverH, INVALID_IND);
    std::vector<size_t> twin = std::vector<size_t>(nCoverH, INVALID_IND);

    // Build next and twin maps for cover halfedges
    for (Halfedge he : mesh.interiorHalfedges()) {
        next[hIdx[he]]           = hIdx[he.next()];
        next[hIdxOpp[he.next()]] = hIdxOpp[he];

        // Handle boundary halfedges
        if (he.sibling() == he || !he.sibling().isInterior()) continue;

        if (he.orientation() != he.sibling().orientation()) {
            // neighboring face orientation agrees
            twin[hIdx[he]]    = hIdx[he.sibling()];
            twin[hIdxOpp[he]] = hIdxOpp[he.sibling()];
        } else {
            // neighboring face orientation disagrees
            twin[hIdx[he]]    = hIdxOpp[he.sibling()];
            twin[hIdxOpp[he]] = hIdx[he.sibling()];
        }
    }

    std::unique_ptr<ManifoldSurfaceMesh> coverMesh;
    std::vector<Halfedge> meshHalfedges;
    std::tie(coverMesh, meshHalfedges) =
        makeManifoldSurfaceMeshAndHalfedgeIndices(next, twin);

    OrientationCoverMapping mapping;
    mapping.coverHalfedgeToBase =
        HalfedgeData<std::pair<Halfedge, bool>>(*coverMesh);

    mapping.baseHalfedgeToCover =
        HalfedgeData<std::pair<Halfedge, Halfedge>>(mesh);

    // record mapping between cover and base mesh
    for (Halfedge he : mesh.interiorHalfedges()) {
        Halfedge positiveLift = meshHalfedges[hIdx[he]];
        Halfedge negativeLift = meshHalfedges[hIdxOpp[he]];

        mapping.coverHalfedgeToBase[positiveLift] = std::make_pair(he, true);
        mapping.coverHalfedgeToBase[negativeLift] = std::make_pair(he, false);

        // halfedge with same orientation comes first
        mapping.baseHalfedgeToCover[he] =
            std::make_pair(positiveLift, negativeLift);
    }

    return std::make_pair(std::move(coverMesh), mapping);
}
