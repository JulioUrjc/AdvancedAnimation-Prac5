
#ifndef _FLUID2_HAS_BEEN_INCLUDED_
#define _FLUID2_HAS_BEEN_INCLUDED_

#include "array2.h"
#include "grid2.h"
#include "particles2.h"

namespace
{
    inline float myrand( void )
    {
        return -1 + 2 * ( (float)rand() ) / RAND_MAX;
    }
}

class Fluid2
{
public:
    Fluid2( const Grid2& grid_ )
		: grid( grid_ ), flipEnabled( true ), particles( 0 )
    {
    }

    void init( void )
    {
		velocityX.resize( grid.getSizeFaces( 0 ) );
		velocityY.resize( grid.getSizeFaces( 1 ) );
		ink.resize( grid.getSize() );
		pressure.resize( grid.getSize() );

        if( flipEnabled )
        {
            oldVelocityX.resize( grid.getSizeFaces( 0 ) );
		    oldVelocityY.resize( grid.getSizeFaces( 1 ) );

            initParticles();
        }
    }

    const Grid2& getGrid( void ) const
    {
        return grid;
    }

	const Array2< float >& getVelocityX( void ) const 
	{
		return velocityX;
	}
	
	const Array2< float >& getVelocityY( void ) const 
	{
		return velocityY;
	}

	const Array2< float >& getPressure( void ) const 
	{
		return pressure;
	}
	
	const Array2< float >& getInk( void ) const 
	{
		return ink;
	}
	
    const Particles2& getParticles( void ) const
    {
        return particles;
    }

	void advanceStep( const float dt )
	{
		fluidEmission();
		fluidAdvection( dt );
		fluidVolumeForces( dt );
		fluidViscosity( dt );
		fluidPressureProjection( dt );
	}

    // init particles
    void initParticles( void );

	// emission
	void fluidEmission( void );

	// advection
	void fluidAdvection( const float dt );

	// external forces
	void fluidVolumeForces( const float dt );

	// viscosity
	void fluidViscosity( const float dt );

	// pressure
	void fluidPressureProjection( const float dt );

public:
    Grid2 grid;

    Array2< float > velocityX;
    Array2< float > velocityY;
    Array2< float > pressure;
    Array2< float > ink;

    bool flipEnabled;

    Array2< float > oldVelocityX;
    Array2< float > oldVelocityY;

    Particles2 particles;
};

#endif
