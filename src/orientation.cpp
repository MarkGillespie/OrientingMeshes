#include "orientation.h"

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

std::pair<std::unique_ptr<ManifoldSurfaceMesh>, OrientationCoverMapping>
constructOrientationCover(SurfaceMesh& mesh, const VertexData<Vector3>& pos) {

    size_t nV = mesh.nVertices();
    size_t nH = mesh.nInteriorHalfedges();
    size_t nF = mesh.nFaces();

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

    std::vector<size_t> heVertexMap, heFaceMap;

    size_t nCoverV, nCoverF;
    std::tie(nCoverV, nCoverF) =
        indexMeshElements(next, twin, heVertexMap, heFaceMap);

    std::vector<std::vector<size_t>> polygons(nCoverF, std::vector<size_t>{});
    std::vector<size_t> faceDegrees(nCoverF);
    std::vector<size_t> halfedgeIndexInFace(nCoverH);
    std::vector<std::vector<std::tuple<size_t, size_t>>> twinFS(
        nCoverF, std::vector<std::tuple<size_t, size_t>>{}); // twin face-side

    auto isValid = [&](size_t iH) -> bool {
        return iH != INVALID_IND && iH < nCoverH;
    };

    // build face vertex list and local halfedge indices
    // initialize face-side twin map, but don't fill until next loop when we
    // know all local indices
    for (size_t iH = 0; iH < nCoverH; iH++) {
        size_t iF = heFaceMap[iH];

        // If we've already built the face's arrays, continue
        if (!polygons[iF].empty()) continue;

        size_t currHe = iH;
        size_t iSide  = 0;
        do {
            size_t iV = heVertexMap[currHe];
            polygons[iF].push_back(heVertexMap[currHe]);
            halfedgeIndexInFace[currHe] = iSide;
            currHe                      = next[currHe];
            iSide++;

            verbose_assert(isValid(currHe),
                           "reached invalid halfedge using next map");
        } while (currHe != iH);

        twinFS[iF] = std::vector<std::tuple<size_t, size_t>>(iSide);
    }

    // build twin map in terms of face-sides
    for (size_t iH = 0; iH < nCoverH; iH++) {
        size_t iF    = heFaceMap[iH];
        size_t iSide = halfedgeIndexInFace[iH];

        if (twin[iH] == INVALID_IND) {
            twinFS[iF][iSide] = std::make_tuple(INVALID_IND, INVALID_IND);
        } else {
            size_t iFTwin     = heFaceMap[twin[iH]];
            size_t iSideTwin  = halfedgeIndexInFace[twin[iH]];
            twinFS[iF][iSide] = std::make_tuple(iFTwin, iSideTwin);
        }
    }

    std::unique_ptr<ManifoldSurfaceMesh> coverMesh;
    coverMesh.reset(new ManifoldSurfaceMesh(polygons, twinFS));

    // Work out new halfedge identities from face-sides
    // Assumes that geometry-central respects input face-sides
    std::vector<Halfedge> meshHalfedges;
    meshHalfedges.reserve(nCoverH);
    for (size_t iH = 0; iH < nCoverH; iH++) {
        Face f       = coverMesh->face(heFaceMap[iH]);
        size_t iSide = halfedgeIndexInFace[iH];
        Halfedge he  = f.halfedge();
        for (size_t iStep = 0; iStep < iSide; iStep++) {
            he = he.next();
        }
        meshHalfedges.push_back(he);
    }


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

std::pair<size_t, size_t> indexMeshElements(const std::vector<size_t>& next,
                                            const std::vector<size_t>& twin,
                                            std::vector<size_t>& heVertex,
                                            std::vector<size_t>& heFace) {

    size_t nH = next.size();
    verbose_assert(twin.size() == nH,
                   "next array has " + std::to_string(next.size()) +
                       " elements, but twin array has " +
                       std::to_string(twin.size()) + " elements");

    std::vector<size_t> prev(nH);
    for (size_t iH = 0; iH < nH; iH++) {
        prev[next[iH]] = iH;
    }

    heVertex = std::vector<size_t>(nH, INVALID_IND);
    heFace   = std::vector<size_t>(nH, INVALID_IND);

    auto isValid = [&](size_t iHe) -> bool {
        return iHe != INVALID_IND && iHe < nH;
    };

    // Loop over halfedges reachable by applying next to iHe and assign them
    // face index iF
    auto labelHeFace = [&](size_t iHe, size_t iF) {
        size_t currHe = iHe;
        do {
            verbose_assert(heFace[currHe] == INVALID_IND,
                           "halfedge '" + std::to_string(currHe) +
                               "' already had a face index?");

            heFace[currHe] = iF;
            currHe         = next[currHe];
        } while (currHe != iHe && isValid(currHe));

        verbose_assert(isValid(currHe),
                       "obtained invalid halfedge while iterating over face?");
    };

    auto labelHeVertex = [&](size_t iHe, size_t iV) {
        size_t currHe = iHe;
        do {
            verbose_assert(heVertex[currHe] == INVALID_IND,
                           "halfedge '" + std::to_string(currHe) +
                               "' already had a vertex index?");

            heVertex[currHe] = iV;
            if (twin[currHe] == INVALID_IND) {
                currHe = INVALID_IND;
            } else {
                currHe = next[twin[currHe]];
            }
        } while (currHe != iHe && isValid(currHe));

        if (!isValid(currHe)) {
            verbose_assert(currHe == INVALID_IND,
                           "obtained out-of-range halfedge " +
                               std::to_string(currHe) +
                               " which is not equal to INVALID_IND (" +
                               std::to_string(INVALID_IND) + ")?");
        }

        // If we hit a boundary, loop back the other way as well
        if (currHe == INVALID_IND) {
            currHe = twin[prev[iHe]];
            while (currHe != iHe && isValid(currHe)) {
                verbose_assert(heVertex[currHe] == INVALID_IND,
                               "halfedge '" + std::to_string(currHe) +
                                   "' already had a vertex index?");

                heVertex[currHe] = iV;
                currHe           = twin[prev[currHe]];
            }

            if (!isValid(currHe)) {
                verbose_assert(currHe == INVALID_IND,
                               "obtained out-of-range halfedge " +
                                   std::to_string(currHe) +
                                   " which is not equal to INVALID_IND (" +
                                   std::to_string(INVALID_IND) + ")?");
            }
        }
    };

    size_t iV = 0;
    size_t iF = 0;
    for (size_t iHe = 0; iHe < nH; iHe++) {
        if (heFace[iHe] == INVALID_IND) {
            labelHeFace(iHe, iF);
            iF++;
        }

        if (heVertex[iHe] == INVALID_IND) {
            labelHeVertex(iHe, iV);
            iV++;
        }
        if (heVertex[iHe] == INVALID_IND) {
            WATCH(twin[iHe]);
            verbose_assert(heVertex[iHe] != INVALID_IND, "???");
        }
    }

    return std::make_pair(iV, iF);
}
