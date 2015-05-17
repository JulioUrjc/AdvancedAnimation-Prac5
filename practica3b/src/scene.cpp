
#include "scene.h"

#include <iostream>
#include <GL/glut.h>

using namespace std;

bool Scene::pauseFlag		= true;
unsigned int Scene::nCellsX = 200;
unsigned int Scene::nCellsY = 200;
float Scene::step			= 0.01f;
float Scene::kDensity       = 1.0f;
float Scene::kGravity       = -1.0f;
float Scene::kViscosity     = 0.001f;
unsigned int Scene::particlesPerCell = 4;

Scene::Scene( void )
	: fluid( NULL ), fluidDisplay( NULL )
{
}

Scene::~Scene( void )
{
	if( fluidDisplay ) delete fluidDisplay;
	if( fluid ) delete fluid;
}

void Scene::printSettings( void )
{
   cerr << endl << "Current Settings:" << endl;
   cerr << "\t-gridsize " << nCellsX << " " << nCellsY << endl;
   cerr << "\t-step " << step << endl;
}

void Scene::init( int argc, char* argv[] )
{
	int arg = 1;
	while( arg < argc )
	{
		// grid size
		if( !strcmp( argv[ arg ], "-gridsize" ) )
		{
			arg++;
		}
		// step
		else if( !strcmp( argv[ arg ], "-step" ) )
		{
			arg++;
		}
		else
		{
			// print instructions
			cerr << endl << "Unrecognized option " << argv[ arg ] << endl;
			cerr << "Usage: practica3.exe -[option1] [settings] -[option2] [settings] ..." << endl;
			cerr << "Options:" << endl;
			cerr << "\t-test [advection|smoke|smoke_with_forces|smoke_complete]" << endl;
			cerr << "\t-gridsize [gridcells x] [gridcells y]" << endl;
			cerr << "\t-step [step size in secs]" << endl;
			break;
		}
	}
   
	init();
}

void Scene::init( void )
{
	printSettings();
	
    const Index2 gridSize( nCellsX, nCellsY );
    const Bbox2  gridDomain( -2.0f, -2.0f, 2.0f,  2.0f );
    const Grid2  grid( gridDomain, gridSize );

    fluid = new Fluid2( grid );
	fluid->init();

    fluidDisplay = new FluidVisualizer2( *fluid );
	fluidDisplay->init();

    // initialize test
    initAnimation();
}

void Scene::initAnimation( void )
{
}

void Scene::pause( void )
{
   pauseFlag = !pauseFlag;
}

void Scene::update( void )
{
   if( pauseFlag ) return;

   animate();
}

void Scene::animate( void )
{
	fluid->advanceStep( step );
}

void Scene::display( void )
{
	// draw ink
	fluidDisplay->drawInkField();
    
    // draw grid
    fluidDisplay->drawGrid();

	// draw velocity field
	fluidDisplay->drawVelocityField();

    // draw partices
    fluidDisplay->drawParticles();
}

