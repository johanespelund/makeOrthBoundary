#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "cellSet.H"
#include "pointField.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Adjust points on specified patches to make edges orthogonal to the internal mesh."
    );

    #include "addOverwriteOption.H"
    argList::addArgument
    (
        "patches",
        "List of patch names or regex - Eg, '(wall inlet outlet)'"
    );
    // Add list of excluded patches
    argList::addArgument
    (
        "excludePatches",
        "List of patch names or regex to exclude - Eg, '(wall inlet outlet)'"
    );

    argList::noFunctionObjects();  // Do not use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word oldInstance = mesh.pointsInstance();
    /* pointField newPoints(mesh.points().size()); */
    pointField newPoints(mesh.points());

    // Get the list of patches to process
    const wordRes patches(args.getList<wordRe>(1));
    const bool overwrite = args.found("overwrite");

    const labelHashSet patchSet(mesh.boundaryMesh().patchSet(patches));
    if (!patchSet.size())
    {
        FatalErrorInFunction
            << "Cannot find any patches in set " << patches << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    // Mark all mesh points on the specified patches
    bitSet vertOnPatch(mesh.nPoints());
    /* const labelList allPatches = mesh.boundaryMesh(); */

    for (const label patchi : mesh.boundaryMesh().patchID())
    {
        // Check if patch type is wedge, symmetry, empty or cyclic
        
        if (mesh.boundaryMesh()[patchi].type() == "wedge"
            || mesh.boundaryMesh()[patchi].type() == "symmetry"
            || mesh.boundaryMesh()[patchi].type() == "empty"
            || mesh.boundaryMesh()[patchi].type() == "cyclic")
        {
            continue;
        }

        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        const labelList& meshPoints = pp.meshPoints();

        vertOnPatch.set(meshPoints);

    }

    // Mark all mesh points on the excluded patches
         // Get the list of patches to process
    const wordRes excludePatches(args.getList<wordRe>(2));
    const labelHashSet excludePatchSet(mesh.boundaryMesh().patchSet(excludePatches));
    if (!excludePatchSet.size())
    {
        FatalErrorInFunction
            << "Cannot find any patches in set " << excludePatches << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    // Create list of points on the excluded patches
    bitSet vertOnExcludePatch(mesh.nPoints());
    for (const label excludePatchi : excludePatchSet)
    {
        const polyPatch& pp = mesh.boundaryMesh()[excludePatchi];
        const labelList& meshPoints = pp.meshPoints();

        vertOnExcludePatch.set(meshPoints);

    }



    // Adjust the points on the patches to be orthogonal to the internal points
    for (const label patchi : patchSet)
    {
        Info << "Adjusting points on patch " << mesh.boundaryMesh()[patchi].name() << endl;
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        const labelList& meshPoints = pp.meshPoints();

        for (const label meshPointi : meshPoints)
        {


            // Skip points on excluded patches
            if (vertOnExcludePatch.test(meshPointi))
            {
               Info << "Skipping point " << meshPointi << " on excluded patch" << endl;
               continue;
             }

            const point& boundaryPoint = mesh.points()[meshPointi];
            const label localPointI = pp.whichPoint(meshPointi);
            Info << "  Adjusting point " << meshPointi << ": " << boundaryPoint << endl;
            
            // Get the list of points connected to the boundary point
            const labelList& pPoints = mesh.pointPoints()[meshPointi];

            // Find the point that is not on the boundary if there is one
            label internalPointi = -1;
            for (const label pPointi : pPoints)
            {
                Info << "    Checking point " << pPointi << endl;
                if (!vertOnPatch.test(pPointi))
                {
                    Info << "    Found internal point " << pPointi << endl;
                    internalPointi = pPointi;
                    break;
                }
            }
            
            if (internalPointi < 0)
            {
                    Info << "No internal point found for boundary point " << meshPointi << endl;    
                    continue;
            }
            else
            {

                // Calculate the patch normal at the boundary point, by taking the average of the normals of the faces connected to the point
                vector boundaryNormal = vector::zero;
                const labelList& faces = pp.pointFaces()[localPointI];
                /* const auto faces = pp.pointFaces(); */
                /* const auto faces = pp.localFaces(); */
                Info << "    Faces: " << faces << endl;
                for (const label facei : faces)
                {
                    /* Info << "    Face " << facei << ": " << faces[facei] << endl; */
                    const vector& faceNormal = pp.faceNormals()[facei];
                    boundaryNormal += faceNormal;
                }
                boundaryNormal.normalise();


                // For now, assumes that the desired normal is (0,-1,0)
                /* vector boundaryNormal = vector(0, -1, 0); */
                vector internalPoint = mesh.points()[internalPointi];
                vector internalToBoundary = boundaryPoint - internalPoint;
                vector movement = internalToBoundary - dot(internalToBoundary, boundaryNormal) * boundaryNormal;
                /* point newPoint = internalPoint + internalToBoundary - dot(internalToBoundary, boundaryNormal) * boundaryNormal; */

                /* mesh.points()[meshPointi] = newPoint; */
                newPoints[meshPointi] = mesh.points()[meshPointi] - movement;
            }
            
            

        }
    }
    mesh.movePoints(newPoints);

    if (!overwrite)
    {
        ++runTime;
    }

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    // Write the modified mesh
    Info << "Writing adjusted mesh to time " << runTime.timeName() << endl;
    mesh.write();

    Info << "End\n" << endl;

    return 0;
}
