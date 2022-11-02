#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"

#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct OrientationCoverMapping {
    HalfedgeData<std::pair<Halfedge, bool>> coverHalfedgeToBase;

    // halfedge with same orientation comes first
    HalfedgeData<std::pair<Halfedge, Halfedge>> baseHalfedgeToCover;
};

std::pair<std::unique_ptr<ManifoldSurfaceMesh>, OrientationCoverMapping>
constructOrientationCover(SurfaceMesh& mesh, const VertexData<Vector3>& pos);


// Assumes next, twin describe the halfedges of a manifold, oriented surface
// mesh
// returns #vertices, #faces
std::pair<size_t, size_t> indexMeshElements(const std::vector<size_t>& next,
                                            const std::vector<size_t>& twin,
                                            std::vector<size_t>& heVertex,
                                            std::vector<size_t>& heFace);
