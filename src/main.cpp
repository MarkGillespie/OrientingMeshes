#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "orientation.h"
#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;

std::unique_ptr<ManifoldSurfaceMesh> orientationCoverMesh;
std::unique_ptr<VertexPositionGeometry> orientationCoverGeom;
OrientationCoverMapping mapping;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh, *psCoverMesh;

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    if (ImGui::Button("Orient")) {
        std::tie(orientationCoverMesh, mapping) =
            constructOrientationCover(*mesh, geom->vertexPositions);
        WATCH(orientationCoverMesh->nBoundaryLoops());

        VertexData<Vector3> positions(*orientationCoverMesh);
        for (Vertex v : orientationCoverMesh->vertices()) {
            Halfedge he        = v.halfedge();
            Halfedge heBase    = mapping.coverHalfedgeToBase[he].first;
            bool heOrientation = mapping.coverHalfedgeToBase[he].second;

            Vertex vBase =
                heOrientation ? heBase.tailVertex() : heBase.tipVertex();
            positions[v] = geom->vertexPositions[vBase];
        }

        orientationCoverGeom.reset(
            new VertexPositionGeometry(*orientationCoverMesh, positions));

        psCoverMesh = polyscope::registerSurfaceMesh(
            "orientation cover", orientationCoverGeom->vertexPositions,
            orientationCoverMesh->getFaceVertexList(),
            polyscopePermutations(*orientationCoverMesh));

        orientationCoverGeom->requireVertexNormals();
        psCoverMesh->addVertexVectorQuantity(
            "normal", orientationCoverGeom->vertexNormals);

        VertexData<bool> isBoundary(*orientationCoverMesh);
        for (Vertex v : orientationCoverMesh->vertices()) {
            isBoundary[v] = v.isBoundary();
        }

        psCoverMesh->addVertexScalarQuantity("isBoundary", isBoundary);

        VertexData<Vector3> displacedPositions =
            positions + 0.1 * orientationCoverGeom->vertexNormals;

        polyscope::registerSurfaceMesh(
            "displaced orientation cover", displacedPositions,
            orientationCoverMesh->getFaceVertexList(),
            polyscopePermutations(*orientationCoverMesh))
            ->addVertexVectorQuantity("normal",
                                      orientationCoverGeom->vertexNormals);
        orientationCoverGeom->unrequireVertexNormals();
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Geometry program");
    args::Positional<std::string> inputFilename(parser, "mesh",
                                                "Mesh to be processed.");

    std::cout << std::boolalpha;

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string filename = "../../meshes/bunny_small.obj";
    // Make sure a mesh name was given
    if (inputFilename) {
        filename = args::get(inputFilename);
    }

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    // Load mesh
    std::tie(mesh, geom) = readSurfaceMesh(filename);

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        polyscope::guessNiceNameFromPath(filename), geom->vertexPositions,
        mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    geom->requireVertexNormals();
    psMesh->addVertexVectorQuantity("normal", geom->vertexNormals);
    geom->unrequireVertexNormals();

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
