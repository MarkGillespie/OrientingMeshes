#include "orientation.h"

#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

std::pair<std::unique_ptr<ManifoldSurfaceMesh>, OrientationCoverMapping>
constructOrientationCover(SurfaceMesh& mesh, const VertexData<Vector3>& pos) {

    // clang-format off
    /*
    // Build a halfedge mesh from connectivity information (0-indexed as always)
    // - `polygons` is the usual vertex indices for each face
    // - `twins` is indices for the halfedge twin pointers. For each halfedge, holds the index of the twin face and
    // halfedge within that face. In each face, the 0'th halfedge goes from vert 0-->1. Use INVALID_IND for boundary.

    ManifoldSurfaceMesh(const std::vector<std::vector<size_t>>& polygons,
                        const std::vector<std::vector<std::tuple<size_t, size_t>>>& twins);
    */
    // clang-format on

    size_t nV = mesh.nVertices();
    size_t nH = mesh.nInteriorHalfedges();
    size_t nF = mesh.nFaces();

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

    WATCH(nH);

    size_t nCoverH = 2 * nH;

    std::vector<size_t> coverNextMap =
        std::vector<size_t>(nCoverH, INVALID_IND);
    std::vector<size_t> coverTwinMap =
        std::vector<size_t>(nCoverH, INVALID_IND);

    // Hook up next and twin maps for cover halfedges
    for (Halfedge he : mesh.interiorHalfedges()) {
        coverNextMap[hIdx[he]]           = hIdx[he.next()];
        coverNextMap[hIdxOpp[he.next()]] = hIdxOpp[he];

        // Handle boundary halfedges
        if (he.sibling() == he || !he.sibling().isInterior()) continue;

        if (he.orientation() != he.sibling().orientation()) {
            // neighboring face orientation agrees
            coverTwinMap[hIdx[he]]    = hIdx[he.sibling()];
            coverTwinMap[hIdxOpp[he]] = hIdxOpp[he.sibling()];
        } else {
            // neighboring face orientation disagrees
            coverTwinMap[hIdx[he]]    = hIdxOpp[he.sibling()];
            coverTwinMap[hIdxOpp[he]] = hIdx[he.sibling()];
        }
    }

    //== Validate twin map
    for (size_t iH = 0; iH < nCoverH; iH++) {
        verbose_assert(coverTwinMap[iH] != iH, "twin has fixed point?");
        if (coverTwinMap[iH] != INVALID_IND) {
            verbose_assert(coverTwinMap[coverTwinMap[iH]] == iH,
                           "twin not involution?");
        }
    }

    std::vector<size_t> heVertexMap, heFaceMap;

    size_t nCoverV, nCoverF;
    std::tie(nCoverV, nCoverF) =
        indexMeshElements(coverNextMap, coverTwinMap, heVertexMap, heFaceMap);

    //== Validate vertex and face indices
    for (size_t iH = 0; iH < nCoverH; iH++) {

        verbose_assert(heFaceMap[iH] != INVALID_IND, "face map err");
        verbose_assert(heFaceMap[iH] == heFaceMap[coverNextMap[iH]],
                       "face map err");

        if (coverTwinMap[iH] == INVALID_IND) continue;

        verbose_assert(heVertexMap[iH] ==
                           heVertexMap[coverNextMap[coverTwinMap[iH]]],
                       "vertex map err");
    }

    std::vector<std::vector<size_t>> polygons(nCoverF, std::vector<size_t>{});
    std::vector<size_t> faceDegrees(nCoverF);
    std::vector<size_t> halfedgeIndexInFace(nCoverH);
    std::vector<std::vector<std::tuple<size_t, size_t>>> twins(
        nCoverF, std::vector<std::tuple<size_t, size_t>>{});

    auto isValid = [&](size_t iH) -> bool {
        return iH != INVALID_IND && iH < nCoverH;
    };

    // build face vertex list and local halfedge indices
    // initialize face-side twin map, but don't fill until next loop when we
    // know all local indices
    for (size_t iH = 0; iH < nCoverH; iH++) {
        size_t iF = heFaceMap[iH];

        verbose_assert(iF != INVALID_IND && iF < nCoverF,
                       "Invalid face index " + std::to_string(iF) +
                           " (there should only be " + std::to_string(nF) +
                           " faces)");

        // If we've already built the face's arrays, continue
        if (!polygons[iF].empty()) continue;

        size_t currHe = iH;
        size_t iSide  = 0;
        do {
            size_t iV = heVertexMap[currHe];
            if (std::find(polygons[iF].begin(), polygons[iF].end(), iV) !=
                polygons[iF].end()) {
                std::cout << "Vertex " << iV
                          << " appears multiple times in polygon " << iF << ": "
                          << polygons[iF] << vendl;
            }
            polygons[iF].push_back(heVertexMap[currHe]);
            halfedgeIndexInFace[currHe] = iSide;
            currHe                      = coverNextMap[currHe];
            iSide++;

            verbose_assert(isValid(currHe),
                           "reached invalid halfedge using next map");
        } while (currHe != iH);

        twins[iF] = std::vector<std::tuple<size_t, size_t>>(iSide);
    }

    // build twin map in terms of face-sides
    for (size_t iH = 0; iH < nCoverH; iH++) {
        size_t iF    = heFaceMap[iH];
        size_t iSide = halfedgeIndexInFace[iH];

        if (coverTwinMap[iH] == INVALID_IND) {
            twins[iF][iSide] = std::make_tuple(INVALID_IND, INVALID_IND);
        } else {
            size_t iFTwin    = heFaceMap[coverTwinMap[iH]];
            size_t iSideTwin = halfedgeIndexInFace[coverTwinMap[iH]];
            twins[iF][iSide] = std::make_tuple(iFTwin, iSideTwin);
        }
    }

    //== Validate polygon and twin lists
    /*
    auto heEndpoints = [&](size_t iF, size_t iS) -> std::pair<size_t, size_t> {
        const std::vector<size_t>& poly = polygons[iF];
        size_t degree                   = poly.size();

        size_t indTail = poly[iS];
        size_t indTip  = poly[(iS + 1) % degree];

        return std::make_pair(indTail, indTip);
    };
    for (size_t iF = 0; iF < twins.size(); iF++) {
        for (size_t iS = 0; iS < twins[iF].size(); iS++) {
            size_t jF = std::get<0>(twins[iF][iS]);
            size_t jS = std::get<1>(twins[iF][iS]);

            std::pair<size_t, size_t> iEndpoints = heEndpoints(iF, iS);
            std::pair<size_t, size_t> jEndpoints = heEndpoints(jF, jS);

            if (iEndpoints.first != jEndpoints.second ||
                iEndpoints.second != jEndpoints.first) {
                WATCH2(iEndpoints, jEndpoints);
                WATCH(polygons[iF]);
                WATCH(polygons[jF]);
                throw_verbose_runtime_error("incompatible twin map");
            }
        }
    }
    */

    std::unique_ptr<ManifoldSurfaceMesh> coverMesh;
    coverMesh.reset(new ManifoldSurfaceMesh(polygons, twins));

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

    // halfedge with same orientation comes first
    mapping.baseHalfedgeToCover =
        HalfedgeData<std::pair<Halfedge, Halfedge>>(mesh);

    for (Halfedge he : mesh.interiorHalfedges()) {
        Halfedge positiveLift = meshHalfedges[hIdx[he]];
        Halfedge negativeLift = meshHalfedges[hIdxOpp[he]];

        mapping.coverHalfedgeToBase[positiveLift] = std::make_pair(he, true);
        mapping.coverHalfedgeToBase[negativeLift] = std::make_pair(he, false);

        mapping.baseHalfedgeToCover[he] =
            std::make_pair(positiveLift, negativeLift);
    }

    if (true) {
        WATCH(coverMesh->nBoundaryLoops());

        HalfedgeData<int> heVertexData(*coverMesh), heFaceData(*coverMesh),
            heTwin(*coverMesh), heNext(*coverMesh);

        for (Halfedge he : mesh.interiorHalfedges()) {
            size_t posId     = hIdx[he];
            size_t negId     = hIdxOpp[he];
            Halfedge posLift = meshHalfedges[hIdx[he]];
            Halfedge negLift = meshHalfedges[hIdxOpp[he]];

            heVertexData[posLift] = heVertexMap[posId];
            heVertexData[negLift] = heVertexMap[negId];
            heFaceData[posLift]   = heFaceMap[posId];
            heFaceData[negLift]   = heFaceMap[negId];
            heTwin[posLift] = meshHalfedges[coverTwinMap[posId]].getIndex();
            heTwin[negLift] = meshHalfedges[coverTwinMap[negId]].getIndex();
            heNext[posLift] = meshHalfedges[coverNextMap[posId]].getIndex();
            heNext[negLift] = meshHalfedges[coverNextMap[negId]].getIndex();
        }

        VertexData<Vector3> positions(*coverMesh);
        for (Vertex v : coverMesh->vertices()) {
            Halfedge he        = v.halfedge();
            Halfedge heBase    = mapping.coverHalfedgeToBase[he].first;
            bool heOrientation = mapping.coverHalfedgeToBase[he].second;

            Vertex vBase =
                heOrientation ? heBase.tailVertex() : heBase.tipVertex();
            positions[v] = pos[vBase];
        }

        std::unique_ptr<VertexPositionGeometry> orientationCoverGeom;
        orientationCoverGeom.reset(
            new VertexPositionGeometry(*coverMesh, positions));

        auto psCoverMesh = polyscope::registerSurfaceMesh(
            "orientation cover", orientationCoverGeom->vertexPositions,
            coverMesh->getFaceVertexList(), polyscopePermutations(*coverMesh));

        orientationCoverGeom->requireVertexNormals();
        psCoverMesh->addVertexVectorQuantity(
            "normal", orientationCoverGeom->vertexNormals);

        VertexData<bool> isBoundary(*coverMesh);
        for (Vertex v : coverMesh->vertices()) {
            isBoundary[v] = v.isBoundary();
        }

        psCoverMesh->addVertexScalarQuantity("isBoundary", isBoundary);

        VertexData<Vector3> displacedPositions =
            positions + 0.1 * orientationCoverGeom->vertexNormals;

        auto psOffsetMesh = polyscope::registerSurfaceMesh(
            "displaced orientation cover", displacedPositions,
            coverMesh->getFaceVertexList(), polyscopePermutations(*coverMesh));
        psOffsetMesh->addVertexVectorQuantity(
            "normal", orientationCoverGeom->vertexNormals);
        psOffsetMesh->addHalfedgeScalarQuantity("v", heVertexData);
        psOffsetMesh->addHalfedgeScalarQuantity("f", heFaceData);
        psOffsetMesh->addHalfedgeScalarQuantity("twin", heTwin);
        psOffsetMesh->addHalfedgeScalarQuantity("next", heNext);
        orientationCoverGeom->unrequireVertexNormals();

        polyscope::show();
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
        // std::cout << "Starting from halfedge " << iHe << vendl;
        size_t currHe = iHe;
        do {
            // std::cout << "\tprocessing halfedge " << currHe << vendl;
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
                // std::cout << "\tprocessing halfedge " << currHe << "
                // (reversed)"
                //           << vendl;
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
