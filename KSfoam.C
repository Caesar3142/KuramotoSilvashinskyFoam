#include "fvCFD.H"

// Main program
int main(int argc, char *argv[])
{
    // Set up the time, mesh, and fields
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField u
    (
        IOobject
        (
            "u",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField v
    (
        IOobject
        (
            "v",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField w
    (
        IOobject
        (
            "w",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
	const dimensionedScalar A ("A", dimensionSet(0,-1, 1, 0, 0, 0,0) , 1);
	const dimensionedScalar B ("B", dimensionSet(0, 1, 0, 0, 0, 0,0) , 1);
	const dimensionedScalar C ("C", dimensionSet(0, 3, 0, 0, 0, 0,0) , 1);

    // Time-stepping loop
    while (runTime.loop())
    {
        // Update time
        // runTime++;

        // Calculate v = d2u/dx2
        v = fvc::laplacian(u);
        w = fvc::laplacian(v); //second order laplacian of u
        // Solve the Kuramoto-Sivashinsky system
        solve
        (
            A*fvm::ddt(u)
          + B*v
          + C*w
        );

        // Write results for each time step
        runTime.write();
    }

    return 0;
}

//*************************************** caesar wiratama ************************************ //
//*************************************** pttensor.com *************************************** //
