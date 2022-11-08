#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

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
