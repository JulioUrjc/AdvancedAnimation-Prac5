
#ifndef _ARRAY2_HAS_BEEN_INCLUDED_
#define _ARRAY2_HAS_BEEN_INCLUDED_

#include "index2.h"

template< class T >
class Array2
{
public:
	Array2()
        : ptr_( NULL )
    {}
    
	Array2( const int i, const int j )
        : ptr_( NULL )
    {
        resize( Index2( i, j ) )
	}

	Array2( const Index2& size )
        : ptr_( NULL )
    {
		resize( size );
	}

	Array2( const Array2< T >& array )
        : ptr_( NULL )
    {
		copy( array );
	}

	~Array2()
    {
		if( ptr_ ) delete[] ptr_;
	}
	
    const Index2& getSize() const
    {
		return size_;
	}

	const T* getData() const
    {
		return ptr_;
	}
	
    unsigned int getLinearIndex( const unsigned int i, const unsigned int j ) const
    {
        return j * size_.x + i;
    }

	bool resize( const Index2& size )
    {
		if( size_ != size )
        {
			if( ptr_ ) delete[] ptr_;

			size_ = size;
			ptr_ = new T[ size.x * size.y ];

            clear();

			return true;
		}
		return false;
	}

	void copy( const Array2< T >& src )
    {
		resize( src.size_ );
        
        for( unsigned int i = 0, n = size_.x * size_.y; i < n; ++i )
			ptr_[ i ] = src[ i ];
	}

	void operator=( const Array2< T >& src )
    {
		copy( src );
	}
    
    const T& getValue( const Index2& i ) const
    {
        const unsigned int id = getLinearIndex( i.x, i.y );
		return ptr_[ id ];
    }
    const T& getValue( const unsigned int i, const unsigned int j ) const
    {
        const unsigned int id = getLinearIndex( i, j );
		return ptr_[ id ];
    }
    
    void setValue( const Index2& i, const T& value )
    {
        const unsigned int id = getLinearIndex( i.x, i.y );
		ptr_[ id ] = value;
    }
    void setValue( const unsigned int i, const unsigned int j, const T& value )
    {
        const unsigned int id = getLinearIndex( i, j );
		ptr_[ id ] = value;
    }

	const T& operator[]( unsigned int i ) const
    {
		return ptr_[ i ];
	}
	T& operator[]( unsigned int i )
    {
		return ptr_[ i ];
	}

    const T& operator[]( const Index2& id ) const
    {
        const unsigned int i = getLinearIndex( id.x, id.y );
		return ptr_[ i ];
	}
	T& operator[]( const Index2& id )
    {
        const unsigned int i = getLinearIndex( id.x, id.y );
		return ptr_[ i ];
	}

    void clear()
    {
        for( int i = 0, n = size_.x * size_.y; i < n; ++i )
			ptr_[ i ] = T();
	}

private:
	Index2 size_;
	T* ptr_;
};

#endif