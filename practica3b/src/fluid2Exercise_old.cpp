
#include "scene.h"
#include "pcg_solver.h"

namespace
{
	//////////////////////////////////////////////
	// Add any custom classes or functions here //
	//////////////////////////////////////////////
    
    inline int clamp( int x, int a, int b )
    {
        return x < a ? a : ( x > b ? b : x );
    }

	struct BaseSampler
	{
		BaseSampler( const Grid2& grid, const Array2< float >& data )
			: grid_( grid ), data_( data )
		{}

        virtual Vec2 getIndex( const Vec2& pos ) = 0;

		float getValue( const Vec2& pos )
		{
            const Vec2 ispos( getIndex( pos ) );
            const Vec2 isposmin( floorf( ispos.x ), floorf( ispos.y ) );
            const Vec2 isposmax( ceilf( ispos.x ), ceilf( ispos.y ) );
            const Vec2 t( ispos - isposmin );
            const Index2 mincorner( (int)isposmin.x, (int)isposmin.y );
            const Index2 maxcorner( (int)isposmax.x, (int)isposmax.y );
            
            const Index2& size = data_.getSize();
            const Index2 id1( clamp( mincorner.x, 0, size.x-1 ), clamp( mincorner.y, 0, size.y-1 ) );
            const Index2 id2( clamp( maxcorner.x, 0, size.x-1 ), clamp( mincorner.y, 0, size.y-1 ) );
            const Index2 id3( clamp( mincorner.x, 0, size.x-1 ), clamp( maxcorner.y, 0, size.y-1 ) );
            const Index2 id4( clamp( maxcorner.x, 0, size.x-1 ), clamp( maxcorner.y, 0, size.y-1 ) );

            const float value = bilerp( data_[ id1 ], data_[ id2 ], data_[ id3 ], data_[ id4 ], t.x, t.y );
			return value;
		}

        float bilerp( const float aa, const float ba, const float ab, const float bb, const float tx, const float ty )
        {
            const float y1 = aa * ( 1.0f - tx ) + ba * tx;
            const float y2 = ab * ( 1.0f - tx ) + bb * tx;
            return y1 * ( 1.0f - ty ) + y2 * ty;
        }

    protected:
        const Grid2& grid_;
		const Array2< float >& data_;
	};
    
	struct CellSampler : BaseSampler
	{
		CellSampler( const Grid2& grid, const Array2< float >& data )
			: BaseSampler( grid, data )
		{}

        Vec2 getIndex( const Vec2& pos )
        {
            return grid_.getCellIndex( pos );
        }
    };

	struct FaceSampler : BaseSampler
	{
		FaceSampler( const Grid2& grid, const Array2< float >& data, const unsigned int axis )
			: BaseSampler( grid, data ), axis_( axis )
		{}
        
        Vec2 getIndex( const Vec2& pos )
        {
            return grid_.getFaceIndex( pos, axis_ );
        }

        const unsigned int axis_;
	};
}

// advection
void Fluid2::fluidAdvection( const float dt )
{
    // ink
    {
	    Array2<float> inkcopy( ink );
	    CellSampler inksampler( grid, inkcopy );

        const Index2& size = ink.getSize();
	    for( unsigned int i = 0; i < size.x; ++i )
	    for( unsigned int j = 0; j < size.y; ++j )
	    {
		    const Index2 id( i, j );

		    const Vec2 pos( grid.getCellPos( id ) );
		    const Vec2 vel( ( velocityX[ id ] + velocityX[ Index2( i+1, j ) ] ) * 0.5f,
                            ( velocityY[ id ] + velocityY[ Index2( i, j+1 ) ] ) * 0.5f );
		    const Vec2 endpos( pos - dt * vel );

		    ink[ id ] = inksampler.getValue( endpos );;
	    }
    }

    // velocity
    {
	    Array2< float > ucopy( velocityX );
        Array2< float > vcopy( velocityY );
	    FaceSampler usampler( grid, ucopy, 0 );
	    FaceSampler vsampler( grid, vcopy, 1 );
        const Index2& sizeu = velocityX.getSize();
        const Index2& sizev = velocityY.getSize();

	    for( unsigned int i = 0; i < sizeu.x; ++i )
	    for( unsigned int j = 0; j < sizeu.y; ++j )
	    {
		    const Index2 id( i, j );
            const Index2 idv1( clamp( i-1, 0, sizev.x-1 ), clamp( j  , 0, sizev.y-1 ) );
            const Index2 idv2( clamp( i  , 0, sizev.x-1 ), clamp( j  , 0, sizev.y-1 ) );
            const Index2 idv3( clamp( i-1, 0, sizev.x-1 ), clamp( j+1, 0, sizev.y-1 ) );
            const Index2 idv4( clamp( i  , 0, sizev.x-1 ), clamp( j+1, 0, sizev.y-1 ) );

		    const Vec2 pos( grid.getFaceXPos( id ) );
		    const Vec2 vel( ucopy[ id ], ( vcopy[ idv1 ] + vcopy[ idv2 ] + vcopy[ idv3 ] + vcopy[ idv4 ] ) * 0.25f );
		    const Vec2 endpos( pos - dt * vel );

		    velocityX[ id ] = usampler.getValue( endpos );
	    }

	    for( unsigned int i = 0; i < sizev.x; ++i )
	    for( unsigned int j = 0; j < sizev.y; ++j )
	    {
		    const Index2 id( i, j );
            const Index2 idu1( clamp( i  , 0, sizeu.x-1 ), clamp( j-1, 0, sizeu.y-1 ) );
            const Index2 idu2( clamp( i  , 0, sizeu.x-1 ), clamp( j  , 0, sizeu.y-1 ) );
            const Index2 idu3( clamp( i+1, 0, sizeu.x-1 ), clamp( j-1, 0, sizeu.y-1 ) );
            const Index2 idu4( clamp( i+1, 0, sizeu.x-1 ), clamp( j  , 0, sizeu.y-1 ) );

		    const Vec2 pos( grid.getFaceYPos( id ) );
		    const Vec2 vel( ( ucopy[ idu1 ] + ucopy[ idu2 ] + ucopy[ idu3 ] + ucopy[ idu4 ] ) * 0.25f, vcopy[ id ] );
		    const Vec2 endpos( pos - dt * vel );

		    velocityY[ id ] = vsampler.getValue( endpos );
	    }
    }
}

// emission
void Fluid2::fluidEmission()
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        Bbox2 source( -0.18f, -1.9f, 0.18f, -1.7f );
        Vec2 ismin = grid.getCellIndex( source.minPosition );
        Vec2 ismax = grid.getCellIndex( source.maxPosition );
        Index2 cornermin( (int) floor( ismin.x ), (int) floor( ismin.y ) );
        Index2 cornermax( (int) ceil( ismax.x ), (int) ceil(ismax.y ) );

		for( unsigned int i = cornermin.x; i <= cornermax.x; ++i )
		for( unsigned int j = cornermin.y; j <= cornermax.y; ++j )
			ink[ Index2( i, j ) ] = 1.0f;
        
        for( unsigned int i = cornermin.x; i <= cornermax.x+1; ++i )
		for( unsigned int j = cornermin.y; j <= cornermax.y; ++j )
            velocityX[ Index2( i, j ) ] = 0.0f;

		for( unsigned int i = cornermin.x; i <= cornermax.x; ++i )
		for( unsigned int j = cornermin.y; j <= cornermax.y+1; ++j )
            velocityY[ Index2( i, j ) ] = 8.0f;
	}
}

// volume forces
void Fluid2::fluidVolumeForces( const float dt )
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        const float dtGravity = dt * Scene::kGravity;
        const Index2& sizev = velocityY.getSize();
		for( unsigned int i = 0, n = sizev.x * sizev.y; i < n; ++i )
			velocityY[ i ] += dtGravity;
	}
}

// viscosity
void Fluid2::fluidViscosity( const float dt )
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        Array2< float > ucopy( velocityX );
        Array2< float > vcopy( velocityY );
        const Index2& sizeu = velocityX.getSize();
        const Index2& sizev = velocityY.getSize();

        const Vec2 dx = grid.getCellDx();
        const Vec2 invDxSq( 1.0f / ( dx.x * dx.x ), 1.0f / ( dx.y * dx.y ) );
        const float dtMuOverRho = dt * Scene::kViscosity / Scene::kDensity;

        for( unsigned int i = 0; i < sizeu.x; ++i )
		for( unsigned int j = 0; j < sizeu.y; ++j )
        {
            const Index2 id( i, j );
            const Index2 id1( clamp( i-1, 0, sizeu.x-1 ), j );
            const Index2 id2( clamp( i+1, 0, sizeu.x-1 ), j );
            const Index2 id3( i , clamp( j-1, 0, sizeu.y-1 ) );
            const Index2 id4( i , clamp( j+1, 0, sizeu.y-1 ) );
            velocityX[ id ] += dtMuOverRho * (
                ( ucopy[ id1 ] - 2.0f * ucopy[ id ] + ucopy[ id2 ] ) * invDxSq.x +
                ( ucopy[ id3 ] - 2.0f * ucopy[ id ] + ucopy[ id4 ] ) * invDxSq.y );
        }
        
		for( unsigned int i = 0; i < sizev.x; ++i )
		for( unsigned int j = 0; j < sizev.y; ++j )
        {
            const Index2 id( i, j );
            const Index2 id1( clamp( i-1, 0, sizev.x-1 ), j );
            const Index2 id2( clamp( i+1, 0, sizev.x-1 ), j );
            const Index2 id3( i , clamp( j-1, 0, sizev.y-1 ) );
            const Index2 id4( i , clamp( j+1, 0, sizev.y-1 ) );
            velocityY[ id ] += dtMuOverRho * (
                ( vcopy[ id1 ] - 2.0f * vcopy[ id ] + vcopy[ id2 ] ) * invDxSq.x +
                ( vcopy[ id3 ] - 2.0f * vcopy[ id ] + vcopy[ id4 ] ) * invDxSq.y );
        }
	}
}

// pressure
void Fluid2::fluidPressureProjection( const float dt )
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        const Vec2 dx = grid.getCellDx();
        const Vec2 invDx( 1.0f / dx.x, 1.0f / dx.y );
        const Vec2 invDxSq( 1.0f / ( dx.x * dx.x ), 1.0f / ( dx.y * dx.y ) );
        
        const Index2& size = pressure.getSize();
        const Index2& sizeu = velocityX.getSize();
        const Index2& sizev = velocityY.getSize();

        // wall boundary conditions
        for( unsigned int j = 0; j < sizeu.y; ++j )
        {
            velocityX[ Index2( 0, j ) ] = 0.0f;
            velocityX[ Index2( sizeu.x-1, j ) ] = 0.0f;
        }
        for( unsigned int i = 0; i < sizev.x; ++i )
        {
            velocityY[ Index2( i, 0 ) ] = 0.0f;
            // TOP as solid
            //velocityY[ Index2( i, sizev.y-1 ) ] = 0.0f;
        }
        
        // rhs
        const float rhoOverDt = Scene::kDensity / dt;
        std::vector< double > rhs( size.x * size.y );
		for( unsigned int i = 0; i < size.x; ++i )
		for( unsigned int j = 0; j < size.y; ++j )
        {
            const Index2 id( i, j );
            rhs[ pressure.getLinearIndex( i, j ) ] = - rhoOverDt * 
                ( ( velocityX[ Index2( i+1, j ) ] - velocityX[ id ] ) * invDx.x + 
                  ( velocityY[ Index2( i, j+1 ) ] - velocityY[ id ] ) * invDx.y );
        }
        
        // A
        SparseMatrix< double > A( size.x * size.y, 5 );
        for( unsigned int i = 0; i < size.x; ++i )
		for( unsigned int j = 0; j < size.y; ++j )
        {
            const unsigned int id = pressure.getLinearIndex( i, j );
            if( i > 0 ) {
                const unsigned int id1 = pressure.getLinearIndex( i-1, j );
                A.add_to_element( id, id, 1. * invDxSq.x );
                A.add_to_element( id, id1, -1. * invDxSq.x ); }
            if( i < size.x-1 ) {
                const unsigned int id1 = pressure.getLinearIndex( i+1, j );
                A.add_to_element( id, id, 1. * invDxSq.x );
                A.add_to_element( id, id1, -1. * invDxSq.x ); }
            if( j > 0 ) {
                const unsigned int id1 = pressure.getLinearIndex( i, j-1 );
                A.add_to_element( id, id, 1. * invDxSq.y );
                A.add_to_element( id, id1, -1. * invDxSq.y ); }
            // TOP as air
            A.add_to_element( id, id, 1. * invDxSq.y );
            if( j < size.y-1 ) {
                const unsigned int id1 = pressure.getLinearIndex( i, j+1 );
                A.add_to_element( id, id1, -1. * invDxSq.y ); }
            //// TOP as wall
            //if( j < size.y-1 ) {
            //    const unsigned int id1 = pressure.getLinearIndex( i, j+1 );
            //    A.add_to_element( id, id, 1. * invDxSq.y );
            //    A.add_to_element( id, id1, -1. * invDxSq.y ); }
        }
        
        // pcg solver
        PCGSolver< double > solver;
        solver.set_solver_parameters( 1e-6, 10000 );
        
        double residual_out;
        int iterations_out;
        std::vector< double > p( size.x * size.y );
        solver.solve( A, rhs, p, residual_out, iterations_out );
        std::cout << "Pressure system result: res=" << residual_out << ", iter=" << iterations_out << std::endl;

        // set pressure
        for( unsigned int i = 0, n = size.x * size.y; i < n; ++i )
            pressure[ i ] = (float) p[ i ];

        // apply pressure gradient
        const float dtOverRho = dt / Scene::kDensity;
		for( unsigned int i = 1; i < sizeu.x - 1; ++i )
		for( unsigned int j = 0; j < sizeu.y; ++j )
        {
            const Index2 id( i, j );
            const float gradp = ( pressure[ id ] - pressure[ Index2( i-1, j ) ] ) * invDx.x;
            velocityX[ id ] -= dtOverRho * gradp;
        }
		for( unsigned int i = 0; i < sizev.x; ++i )
		for( unsigned int j = 1; j < sizev.y - 1; ++j )
        {
            const Index2 id( i, j );
            const float gradp = ( pressure[ id ] - pressure[ Index2( i, j-1 ) ] ) * invDx.y;
            velocityY[ id ] -= dtOverRho * gradp;
        }
        // apply pressure gradient: TOP as air
        for( unsigned int i = 0; i < sizev.x; ++i )
        {
            const Index2 id( i, sizev.y-1 );
            const float gradp = ( 0.0f - pressure[ Index2( i, size.y-1 ) ] ) * invDx.y;
            velocityY[ id ] -= dtOverRho * gradp;
        }
	}
}
