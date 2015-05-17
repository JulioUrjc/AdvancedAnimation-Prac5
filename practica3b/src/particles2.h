
#ifndef _PARTICLES2_HAS_BEEN_INCLUDED_
#define _PARTICLES2_HAS_BEEN_INCLUDED_

#include "vec2.h"
#include <vector>

class Particles2
{
public:
    Particles2( unsigned int size=0 )
    {
        size_ = size;
        position_.resize( size_ );
        velocity_.resize( size_ );
        ink_.resize( size_ );
    }
    
    unsigned int getSize() const
    {
        return size_;
    }

    void resize( const unsigned int size )
    {
        if( size_ > size )
        {
            size_ = size;
            position_.resize( size_ );
            velocity_.resize( size_ );
            ink_.resize( size_ );
        }
        else if( size_ < size )
        {
            position_.resize( size );
            velocity_.resize( size );
            ink_.resize( size );
            
            const unsigned int diff = size - size_;
            std::memset( &position_[ size_ ], 0, sizeof(Vec2) * diff );
            std::memset( &velocity_[ size_ ], 0, sizeof(Vec2) * diff );
            std::memset( &ink_[ size_ ], 0, sizeof(float) * diff );
            size_ = size;
        }
    }

    const Vec2& getPosition( const unsigned int i ) const
    {
        return position_[ i ];
    }

    const Vec2& getVelocity( const unsigned int i ) const
    {
        return velocity_[ i ];
    }

    const float getInk( const unsigned int i ) const
    {
        return ink_[ i ];
    }

    void setPosition( const unsigned int i, const Vec2& pos )
    {
        position_[ i ] = pos;
    }
    
    void setVelocity( const unsigned int i, const Vec2& vel )
    {
        velocity_[ i ] = vel;
    }

    void setInk( const unsigned int i, const float s )
    {
        ink_[ i ] = s;
    }

    unsigned int addParticle(
        const Vec2& pos, const Vec2& vel = Vec2(), const float s = 0 )
    {
        unsigned int id = size_++;
        
        position_.push_back( pos );
        velocity_.push_back( vel );
        ink_.push_back( s );

        return id;
    }

private:
    unsigned int size_;
    std::vector<Vec2> position_;
    std::vector<Vec2> velocity_;
    std::vector<float> ink_;
};

#endif
